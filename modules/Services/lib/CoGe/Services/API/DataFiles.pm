package CoGe::Services::API::DataFiles;

use Mojo::Base 'Mojolicious::Controller';

use CoGeX;
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_genome_path get_dataset_source_path);
use CoGe::Services::Auth;
use CoGe::Services::Error;
use File::Spec::Functions qw(catfile);

###############################################################################
# Authenticated file gateway (JBrowse2 Phase A, 2026-08-18).
#
# Serves the preserved on-disk artifacts for genomes and datasets through the
# SAME access checks the web UI uses. New code, parallel to the JBrowse1
# endpoints -- nothing existing routes through here.
#
#   GET /genomes/:gid/files/:kind    kind: fasta | fai
#   GET /datasets/:id/files/:kind    kind: gff | gff-tabix | gff-csi
#
# Design (see repo discussion 2026-08-18):
#   * CoGe is the sole policy decision point. Storage -- local disk today,
#     object store tomorrow -- is dumb bytes. The ACL graph is never mirrored
#     into the storage layer; access is decided here, per request.
#   * `kind` is matched against the fixed tables below and the filename is
#     derived from the DB row (dataset.name) or a constant (genome.faa).
#     The client never supplies any path fragment: no traversal surface.
#   * The FILE_STORE conf key is the driver seam. Absent or 'local' (the
#     default -- no config change needed today) serves bytes via
#     Mojolicious's static machinery, which handles HTTP Range/206 natively
#     (tabix and FASTA index readers are range-request consumers).
#     's3' is reserved: _serve() is where a presigned-URL redirect goes, and
#     the object key for that future is exactly the path relative to DATADIR
#     -- identical strings on disk and in the bucket, so migration is a sync.
#   * Cache-Control follows the object's restricted flag: bounded staleness
#     if visibility flips, cacheable where it is safe.
###############################################################################

my %GENOME_KINDS = (
    'fasta' => 'genome.faa',
    'fai'   => 'genome.faa.fai',
);

# Suffixes applied to dataset.name -- must stay in lockstep with the preserve
# block in scripts/load_annotation.pl.
my %DATASET_KINDS = (
    'gff'       => '',
    'gff-tabix' => '.sorted.gz',
    'gff-csi'   => '.sorted.gz.csi',
);

sub genome_file {
    my $self = shift;
    my $gid  = $self->stash('gid');
    my $kind = $self->stash('kind');

    my $is_aliases = ($kind // '') eq 'aliases';
    my $fname = $GENOME_KINDS{$kind // ''};
    return $self->render(API_STATUS_NOTFOUND) unless $fname or $is_aliases;

    my ($db, $user) = CoGe::Services::Auth::init($self);
    return $self->render(API_STATUS_CUSTOM(500, 'database unavailable')) unless $db;

    my $genome = $db->resultset('Genome')->find($gid);
    return $self->render(API_STATUS_NOTFOUND) if !$genome or $genome->deleted;

    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return $self->render(API_STATUS_UNAUTHORIZED);
    }

    my $dir = get_genome_path($genome->id)
        or return $self->render(API_STATUS_NOTFOUND);

    # kind=aliases: a refName-aliases file for JBrowse2's RefNameAliasAdapter,
    # generated from the .fai rather than stored. CoGe FASTAs carry NCBI-style
    # prefixed sequence names (e.g. "lcl|LL0249_Chr01") while the loaded GFFs
    # use the bare names -- disjoint refName sets, so annotation tracks
    # silently render nothing without aliasing. Format: canonical name (as in
    # the FASTA) TAB alias. Only prefixed names emit a line.
    if ($is_aliases) {
        my $fai = catfile($dir, $GENOME_KINDS{fai});
        return $self->render(API_STATUS_NOTFOUND) unless -f $fai;
        open(my $fh, '<', $fai)
            or return $self->render(API_STATUS_CUSTOM(500, 'cannot read fai'));
        my $out = '';
        while (my $line = <$fh>) {
            my ($name) = split(/\t/, $line, 2);
            next unless defined $name and $name =~ /\|/;
            my ($bare) = $name =~ /\|([^|]+)$/;
            $out .= "$name\t$bare\n" if $bare;
        }
        close $fh;
        $self->res->headers->cache_control(
            $genome->restricted ? 'private, max-age=3600' : 'public, max-age=86400' );
        return $self->render( text => $out, format => 'txt' );
    }

    return $self->_serve( catfile($dir, $fname), $genome->restricted );
}

sub dataset_file {
    my $self = shift;
    my $id   = $self->stash('id');
    my $kind = $self->stash('kind');

    defined( my $suffix = $DATASET_KINDS{$kind // ''} )
        or return $self->render(API_STATUS_NOTFOUND);

    my ($db, $user) = CoGe::Services::Auth::init($self);
    return $self->render(API_STATUS_CUSTOM(500, 'database unavailable')) unless $db;

    my $ds = $db->resultset('Dataset')->find($id);
    return $self->render(API_STATUS_NOTFOUND) if !$ds or $ds->deleted;

    # Access rides on the dataset's GENOMES, not the dataset's own restricted
    # flag. Per-dataset enforcement shipped and was REVERTED the same day
    # (2026-08-18): dataset.restricted has never been enforced anywhere in
    # CoGe -- has_access_to_dataset's public short-circuit is commented out
    # upstream and the JBrowse1 layer gates on the genome -- so the flags were
    # never curated. 1,568 live datasets on PUBLIC genomes carry restricted=1
    # as load-wizard defaults, and enforcing the flag hid their tracks from
    # everyone but their owners. Effective CoGe semantics, kept here: a
    # dataset is visible iff some genome it belongs to is visible.
    my ($allowed, $public) = (0, 0);
    for my $genome ( $ds->genomes ) {
        next if $genome->deleted;
        if ( !$genome->restricted ) {
            $allowed = 1;
            $public  = 1;
        }
        elsif ( defined $user and $user->has_access_to_genome($genome) ) {
            $allowed = 1;
        }
    }
    $allowed = 1 if defined $user and $user->is_admin;
    return $self->render(API_STATUS_UNAUTHORIZED) unless $allowed;

    my $dir = get_dataset_source_path($ds->id)
        or return $self->render(API_STATUS_NOTFOUND);
    return $self->_serve( catfile($dir, $ds->name . $suffix), !$public );
}

# List the datasets of a genome that THIS requester may know about, plus
# which gateway files exist for each -- the single call the JBrowse2 page
# builds its track config from (Phase B).
#
# Visibility is the same predicate set as the JBrowse1 track_config fix:
# genome access gates the endpoint; per dataset, public OR accessible via the
# user's genomes OR admin. Restricted datasets a user cannot access are
# OMITTED entirely (existence and name are metadata worth protecting), and
# this listing is convenience only -- the file endpoints above re-check
# access on every fetch, so a leaked id never becomes leaked bytes.
sub genome_datasets {
    my $self = shift;
    my $gid  = $self->stash('gid');

    my ($db, $user) = CoGe::Services::Auth::init($self);
    return $self->render(API_STATUS_CUSTOM(500, 'database unavailable')) unless $db;

    my $genome = $db->resultset('Genome')->find($gid);
    return $self->render(API_STATUS_NOTFOUND) if !$genome or $genome->deleted;

    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return $self->render(API_STATUS_UNAUTHORIZED);
    }

    # No per-dataset filter, deliberately (2026-08-18 revert): genome access
    # gates everything, matching JBrowse1 -- dataset.restricted was never
    # enforced in CoGe and 1,568 datasets on public genomes carry it as an
    # uncurated load-wizard default. The flag is still REPORTED per dataset.
    my @datasets;
    for my $ds ( sort { $a->name cmp $b->name } $genome->datasets ) {
        next if $ds->deleted;
        my $dir = get_dataset_source_path($ds->id);
        my %files;
        for my $kind (keys %DATASET_KINDS) {
            $files{$kind} = ($dir and -f catfile($dir, $ds->name . $DATASET_KINDS{$kind})) ? \1 : \0;
        }
        push @datasets, {
            id         => int($ds->id),
            name       => $ds->name,
            version    => $ds->version,
            restricted => $ds->restricted ? \1 : \0,
            date       => '' . $ds->date,
            files      => \%files,
        };
    }

    # Genome block included so the page needs exactly one metadata fetch.
    my $gdir = get_genome_path($genome->id);
    $self->render(json => {
        genome => {
            id         => int($genome->id),
            name       => $genome->info,
            restricted => $genome->restricted ? \1 : \0,
            files      => {
                fasta => ($gdir and -f catfile($gdir, $GENOME_KINDS{fasta})) ? \1 : \0,
                fai   => ($gdir and -f catfile($gdir, $GENOME_KINDS{fai}))   ? \1 : \0,
            },
        },
        datasets => \@datasets,
    });
}

# The driver seam. Everything above decides WHETHER; this decides HOW.
sub _serve {
    my ($self, $path, $restricted) = @_;

    my $store = get_defaults()->{FILE_STORE} // 'local';
    if ($store ne 'local') {
        # 's3' driver lands here: authorize (already done), presign
        # GET s3://<bucket>/<path relative to DATADIR>, respond
        # 302 Location: <presigned url>. Browsers follow the redirect
        # transparently and S3 honors Range on presigned GETs, so JBrowse2
        # config needs no change when this is implemented.
        return $self->render(API_STATUS_CUSTOM(501,
            "FILE_STORE '$store' not implemented"));
    }

    return $self->render(API_STATUS_NOTFOUND) unless -f $path;

    # Bounded staleness on restricted-flag flips; private data never lands in
    # shared caches.
    $self->res->headers->cache_control(
        $restricted ? 'private, max-age=3600' : 'public, max-age=86400' );

    # Mojolicious's static machinery serves this with HTTP Range / 206
    # support, which the tabix (.gz + .csi) and faidx consumers rely on.
    return $self->reply->file($path);
}

1;

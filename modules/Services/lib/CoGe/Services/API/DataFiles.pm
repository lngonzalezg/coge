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

    my $fname = $GENOME_KINDS{$kind // ''}
        or return $self->render(API_STATUS_NOTFOUND);

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

    # NB: has_access_to_dataset does NOT grant public datasets (its
    # !restricted short-circuit is commented out upstream) -- callers must
    # test the flag themselves, which is exactly what the web UI does.
    if ( $ds->restricted
        and ( not defined $user or not $user->has_access_to_dataset($ds) ) )
    {
        return $self->render(API_STATUS_UNAUTHORIZED);
    }

    my $dir = get_dataset_source_path($ds->id)
        or return $self->render(API_STATUS_NOTFOUND);
    return $self->_serve( catfile($dir, $ds->name . $suffix), $ds->restricted );
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

    # Accessible-dataset set computed once (same reasoning as track_config:
    # has_access_to_dataset is uncached and quadratic in a loop).
    my %ds_visible;
    if ($user) {
        if ($user->is_admin) {
            $ds_visible{$_->id} = 1 for $genome->datasets;
        }
        else {
            for my $g ($user->genomes(include_deleted => 1)) {
                $ds_visible{$_->id} = 1 for $g->datasets;
            }
        }
    }

    my @datasets;
    for my $ds ( sort { $a->name cmp $b->name } $genome->datasets ) {
        next if $ds->deleted;
        next if $ds->restricted and not $ds_visible{$ds->id};
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

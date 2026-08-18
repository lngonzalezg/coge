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

package CoGe::Builder::IndexAnnotation;

use Moose;
extends 'CoGe::Builder::Buildable';

use File::Spec::Functions qw(catdir catfile);
use CoGe::Core::Storage qw(get_dataset_source_path);

###############################################################################
# Backfill the preserved-annotation triplet for a dataset loaded before file
# preservation existed, so JBrowse2 can render it (job type index_annotation,
# submitted by GenomeView2's "Visualize older annotations" button).
#
# One task: scripts/backfill_dataset_gff.pl, which exports the dataset from
# the database (Dataset->gff), sanitizes, sorts, bgzips and CSI-indexes into
# DATADIR/annotation/<tiered dsid>/. The script is idempotent and race-safe;
# outputs declared here let JEX short-circuit if the files already exist.
###############################################################################

sub get_name {
    my $self = shift;
    my $ds = $self->request->dataset;
    return 'Index annotation "' . $ds->name . '" (dataset ' . $ds->id . ') for the genome browser';
}

sub build {
    my $self = shift;
    my $ds   = $self->request->dataset;

    my $dir = get_dataset_source_path($ds->id);
    $self->add({
        cmd  => catfile($self->conf->{SCRIPTDIR}, 'backfill_dataset_gff.pl'),
        args => [
            ['-dsid',   $ds->id, 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0],
        ],
        inputs  => [],
        outputs => [
            catfile($dir, $ds->name . '.sorted.gz'),
            catfile($dir, $ds->name . '.sorted.gz.csi'),
        ],
        description => 'Exporting and indexing annotation for the genome browser',
    });

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;

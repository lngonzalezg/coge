package CoGe::Request::Dataset;

use Moose;
extends 'CoGe::Request::Request';

use CoGe::Exception::AccessDenied;
use CoGe::Exception::ItemNotFound;
use CoGe::Exception::MissingField;

has dataset => (is => 'rw', isa => 'CoGeX::Result::Dataset');

sub is_valid { # called first
    my $self = shift;

    my $dsid = $self->parameters->{dsid} || $self->parameters->{dataset_id};
    unless ($dsid) {
        CoGe::Exception::MissingField->throw(message => "Missing dataset_id");
    }

    my $ds = $self->db->resultset("Dataset")->find($dsid);
    if (!$ds or $ds->deleted) {
        CoGe::Exception::ItemNotFound->throw(type => 'dataset', id => $dsid);
    }
    $self->dataset($ds);

    return 1;
}

sub has_access { # called second
    my $self = shift;

    # Genome-gated, matching CoGe's effective dataset semantics everywhere
    # (JBrowse1 track_config, the DataFiles gateway): a dataset is accessible
    # iff some genome it belongs to is. dataset.restricted itself is
    # deliberately NOT enforced -- the flags were never enforced historically
    # and are uncurated (see the 2026-08-18 revert in DataFiles.pm).
    return 1 if $self->user and $self->user->is_admin;
    for my $genome ($self->dataset->genomes) {
        next if $genome->deleted;
        return 1 if !$genome->restricted;
        return 1 if $self->user and $self->user->has_access_to_genome($genome);
    }
    CoGe::Exception::AccessDenied->throw(type => 'dataset', id => $self->dataset->id);
}

__PACKAGE__->meta->make_immutable;

1;

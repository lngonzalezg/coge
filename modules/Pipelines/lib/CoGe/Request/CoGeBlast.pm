package CoGe::Request::CoGeBlast;

use Moose;
extends 'CoGe::Request::Request';

use CoGe::Builder::Tools::CoGeBlast qw( get_genomes );

sub is_valid {
    my $self = shift;
    return unless $self->parameters->{query_seq};
    return unless scalar get_genomes($self->parameters->{genomes}, $self->parameters->{notebooks}, $self->db);
    return 1;
}

sub has_access {
    my $self = shift;
    my @gids = get_genomes($self->parameters->{genomes}, $self->parameters->{notebooks}, $self->db);
    for (@gids) {
        if ($self->user) {
           return unless $self->user->has_access_to_genome($self->db->resultset("Genome")->find($_));
        }
        else {
            # Security pass 1 (S3): was raw SQL concatenating the genome id into the
            # anonymous-user restricted-genome check -- an injection here also defeated
            # the access gate. Genome ids are now integer-validated in get_genomes, and
            # this uses the same ORM lookup as the authenticated branch above (bound).
            my $genome = $self->db->resultset("Genome")->find($_);
            return if $genome && $genome->restricted;
        }
    }
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;

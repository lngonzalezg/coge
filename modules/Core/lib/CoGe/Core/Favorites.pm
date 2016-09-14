package CoGe::Core::Favorites;
use strict;
use warnings;

use Moose;

# Attributes -------------------------------------------------------------------
has 'user' => (
    is       => 'ro',
    required => 1
    #isa      =>
);

# Private methods --------------------------------------------------------------
sub BUILD { # called after constructor
    my $self = shift;
    
    if ($self->user->is_public) {
        # error
        $self->user(undef);
        return;
    }
}

# Public methods ---------------------------------------------------------------
sub is_favorite {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return unless $item;
    
    my ($favorite) = $self->user->favorites({ child_id => $item->id, child_type => $item->item_type });
    return defined($favorite);
}

sub get {
    my $self = shift;
    my @favorites = map { $_->child } $self->user->favorites();
    return wantarray ? @favorites : \@favorites;
}

sub add {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return if $self->is_favorite($item);
    
    my $favorite = $self->user->add_to_favorites({ child_id => $item->id, child_type => $item->item_type });
    return $favorite;
}

sub remove {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return unless $self->is_favorite($item);
    
    my ($favorite) = $self->user->favorites({ child_id => $item->id, child_type => $item->item_type });
    $favorite->delete;
    return 1;
}

sub toggle {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    
    if ($self->is_favorite($item)) {
        $self->remove($item);
        return 0;
    }
    else {
        $self->add($item);
        return 1;
    }
}

1;

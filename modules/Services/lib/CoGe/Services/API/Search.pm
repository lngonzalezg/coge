package CoGe::Services::API::Search;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::API::Job;
use CoGe::Core::Search;
use CoGe::Services::Error;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(API_STATUS_SEARCHTERM);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $show_users = 0;
    if ($user && $user->is_admin) {
		$show_users = 1;
    }

    # Optional paging (per result category): the UI requests one page at a
    # time instead of everything up to the server-side cap.
    my $limit  = $self->param('limit');
    my $offset = $self->param('offset');
    $limit  = (defined $limit  && $limit  =~ /^\d+$/) ? int($limit)  : undef;
    $offset = (defined $offset && $offset =~ /^\d+$/) ? int($offset) : undef;

    my %fetched; # per-category pre-filter window fill, drives UI paging
    my @results = CoGe::Core::Search::search(
        db => $db,
        user => $user,
        search_term => $search_term,
        show_users => $show_users,
        limit => $limit,
        offset => $offset,
        fetched => \%fetched
    );

    $self->render(json => { results => \@results, fetched => \%fetched });
}

1;

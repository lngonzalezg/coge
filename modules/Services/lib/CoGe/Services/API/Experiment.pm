package CoGe::Services::API::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(decode_json);
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::Utils;
use CoGe::Core::Experiment qw( delete_experiment );
use CoGe::Core::Favorites;
use CoGe::Core::Metadata qw( create_annotation delete_annotation get_annotation get_annotations );
use CoGe::Services::Auth;
use CoGe::Services::API::Job;
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

    # Search experiments
    my $search_term2 = '%' . $search_term . '%';
    my @experiments = $db->resultset("Experiment")->search(
        \[
            'experiment_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'experiment_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
         !$_->restricted || (defined $user && $user->has_access_to_experiment($_))
    } @experiments;

    # Get user's favorites
    my $favorites = CoGe::Core::Favorites->new(user => $user);

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
        restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        favorited  => $favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false
      }
    } @filtered;

    $self->render(json => { experiments => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 0, $db, $user);
    return unless $experiment;

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group ? { name => $_->type->group->name, description => $_->type->group->description } : undef
        }
    } $experiment->annotations;

    # Format types
    my @types = map {
        {
            name => $_->name
        }
    } $experiment->types;

    $self->render(json => {
        id => int($experiment->id),
        name => $experiment->name,
        description => $experiment->description,
        version => $experiment->version,
        genome_id  => int($experiment->genome->id),
        source => {
            name => $experiment->source->name,
            description => $experiment->source->description,
            link => $experiment->source->link
        },
        types => \@types,
        additional_metadata => \@metadata,
        restricted => $experiment->restricted ? Mojo::JSON->true : Mojo::JSON->false,
		num_items => $experiment->row_count
    });
}

sub fetch_annotations {
    my $self = shift;
    my $id = int($self->stash('id'));
    my ($db) = CoGe::Services::Auth::init($self);
    #TODO add error checking on ID parameter

    $self->render(json => get_annotations($id, 'Experiment', $db, 1));
}

sub fetch_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));
    #TODO add error checking on ID parameters

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 0, $db, $user);
    return unless $experiment;

    my $annotation = get_annotation($aid, 'Experiment', $db);
    $self->render(json => $annotation) if $annotation;
}

sub add {
    my $self = shift;
    my $data = $self->req->body; #$self->req->json; # mdb replaced 11/22/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    unless ($data) {
        $self->render(API_STATUS_MISSING_BODY);
        return;
    }
    $data = decode_json($data);

    # Valid data items
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(API_STATUS_MISSING_DATA);
        warn Dumper $data->{source_data};
        return;
    }
    
    # Marshall incoming payload into format expected by Job Submit.
    # Note: This is kind of a kludge -- is there a better way to do this using
    # Mojolicious routing?
    my $request = {
        type => 'load_experiment',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

sub add_annotation {
    my $self = shift;
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    #TODO add error checking on parameters
    create_annotation(
        conf => $conf,
        db => $db,
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => int($self->stash('id')),
        target_type => 'experiment',
        text => $self->param('annotation'),
        type_name => $self->param('type_name'),
        user => $user
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub delete_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $error = CoGe::Core::Metadata::delete_annotation($aid, $id, 'Experiment', $db, $user);
    if ($error) {
        $self->render(status => 400, json => { error => { message => $error} });
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub remove {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);

    my $error = delete_experiment($id, $db, $user);
    if ($error) {
        $self->render(status => 400, json => { error => { message => $error} });
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub update {
	my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $experiment = $self->_get_experiment($id, 1, $db, $user);
    return unless $experiment;

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$experiment->update($data->{metadata});
	$self->render(json => { success => Mojo::JSON->true });
}

sub update_annotation {
    my $self = shift;
    my ($db, $user) = CoGe::Services::Auth::init($self);
    CoGe::Core::Metadata::update_annotation(
        annotation_id => int($self->stash('aid')),
        db => $db,
        delete_bisque_image => $self->param('delete_bisque_image'),
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => $self->stash('id'),
        target_type => 'experiment',
        text => $self->param('annotation'),
        type_name => $self->param('type_name'),
        user => $user
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub _get_experiment {
    my ($self, $id, $own_or_edit, $db, $user) = @_;
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }
    if ($own_or_edit) {
        unless ($user) {
            $self->render(status => 401, json => { error => { message => "User not logged in"} });
            return;
        }
        return $experiment if $user->is_admin;
        unless ($user->is_owner_editor(experiment => $id)) {
            $self->render(API_STATUS_UNAUTHORIZED);
            return;
        }
    }
    unless ( !$experiment->restricted || (defined $user && $user->has_access_to_experiment($experiment)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    return $experiment;
}

1;

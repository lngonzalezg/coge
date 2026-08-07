package CoGe::Services::API::Log;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use Data::Dumper;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub fetch {
    my $self = shift;
    my $id   = $self->stash('id');
    my $type = $self->stash('type');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    # Fetch log entries
    my @entries;
    my $node_types = CoGeX::node_types();
    my $type_code = $node_types->{$type};

    # §7.8 require access to the target object before returning its activity log
    # (entries leak descriptions and full-URL links). Content objects are gated
    # by object access (public objects remain viewable to any logged-in user);
    # any other node type requires admin.
    if ($type_code) {
        my %rsname = (genome => 'Genome', experiment => 'Experiment', notebook => 'List', list => 'List');
        my $rs = $rsname{lc($type)};
        if ($rs) {
            my $object = $db->resultset($rs)->find($id);
            unless ($object && ($user->is_admin || $user->has_access_to($object))) {
                $self->render(API_STATUS_UNAUTHORIZED);
                return;
            }
        }
        elsif (!$user->is_admin) {
            $self->render(API_STATUS_UNAUTHORIZED);
            return;
        }
    }

    @entries = $db->resultset('Log')->search({ parent_id => $id, parent_type => $type_code }) if $type_code;
    
    my @formatted = map { 
        {   id => $_->id, 
            description => $_->description,
            time => $_->time,
            link => $_->link,
            user => {
                id => $_->user_id,
                name => $_->user->display_name,
                image_id => $_->user->image_id
            },
            parent_id => $_->id,
            parent_type => CoGeX->node_type_name($_->parent_type)
        } 
    } @entries;

    $self->render(json => {
        id => int($id),
        type => $type,
        entries => \@formatted
    });
}

1;

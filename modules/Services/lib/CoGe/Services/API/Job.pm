package CoGe::Services::API::Job;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::Asset::File;
use Mojo::JSON qw(decode_json);
#use IO::Compress::Gzip 'gzip';
use Data::Dumper;
use File::Spec::Functions qw( catfile );
use CoGe::Services::Auth qw( init );
use CoGe::Services::Error;
use CoGe::Core::Storage qw( get_workflow_paths get_workflow_results );
use CoGe::Factory::RequestFactory;
use CoGe::Factory::PipelineFactory;

sub add {
    my $self = shift;
    my $payload = shift; # allow special payload to be passed in from other controllers
    my $json = $self->req->body; #$self->req->json; # mdb replaced 11/30/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    warn "CoGe::Services::API::Job::add\n", Dumper $payload, Dumper $json;
    unless ($payload || $json) {
        $self->render(API_STATUS_MISSING_BODY);
        return;
    }
    $payload = decode_json($json) if !$payload;

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    
    # Create request
    my $request = CoGe::Factory::RequestFactory->new(db => $db, conf => $conf, user => $user)->get($payload);

    # Validate the request's parameters
    unless ($request and $request->is_valid) {
        return $self->render(API_STATUS_BAD_REQUEST("Invalid request"));
    }

    # Check if authentication is required
    if ($request->authRequired && !$user) {
        return $self->render(API_STATUS_UNAUTHORIZED);
    }

    # Check user's permission to execute the request
    unless ($request->has_access) {
        return $self->render(API_STATUS_UNAUTHORIZED);
    }

    # Create pipeline to execute job
    my $pipeline = CoGe::Factory::PipelineFactory->new()->get($request);
    unless ($pipeline && $pipeline->workflow) {
        return $self->render(API_STATUS_CUSTOM(200, "Failed to generate pipeline")); # this will probably never happen, will be preempted by an exception
    }

    # Submit pipeline
    my $response = $pipeline->submit();
    unless ($response->{success} && $response->{id}) {
        return $self->render(API_STATUS_CUSTOM(200, "Failed to start workflow"));
    }
    
    # Log submission
    CoGe::Accessory::Web::log_history(
        db          => $db,
        parent_id   => $response->{id},
        parent_type => 7, #FIXME magic number
        user_id     => ($user ? $user->id : 0),
        page        => ($pipeline->page ? $pipeline->page : 'API'), # will be 'API' for external API requests
        description => $pipeline->workflow->name,
        link        => ($response->{site_url} ? $response->{site_url} : '')
    );
    print STDERR "CoGe::Services::API::Job::add submitted workflow ", $response->{id}, "\n";
    
    # Convert 'success' to boolean
    $response->{success} = Mojo::JSON->true;

    return $self->render(status => 201, json => $response);
}

sub fetch {
    my $self = shift;
    my $id = $self->stash('id');
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required -- removed 5/26/2016 because synmap3d needs public access
    # unless (defined $user) {
    #     $self->render(status => 401, json => {
    #         error => { Auth => "Access denied" }
    #     });
    #     return;
    # }

    # Get job status from JEX
    my $jex = CoGe::JEX::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    my $job_status = $jex->get_job($id);
    unless ($job_status) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

#    unless ($job_status->{status} =~ /completed|running/i) {
#        $self->render(json => {
#            id => int($id),
#            status => $job_status->{status}
#        });
#        return;
#    }

    # TODO COGE-472: add read from "debug.log" in results path if JEX no longer has the log in memory

    # Add tasks (if any)
    my @tasks;
    foreach my $task (@{$job_status->{jobs}}) {
        my $t = {
            started => $task->{started},
            ended => $task->{ended},
            elapsed => $task->{elapsed},
            description => $task->{description},
            status => lc($task->{status}),
            log => undef
        };

        if (defined $task->{output}) {
            foreach (split(/\\n/, $task->{output})) {
                #print STDERR $_, "\n";
                next unless ($_ =~ /^'?log\: /);
                $_ =~ s/^'?log\: //;
                $t->{log} .= $_ . "\n";
            }
        }

        push @tasks, $t;
    }

    # Add results
    my $user_name = $user ? ($user->is_admin ? undef : $user->name) : 'public'; # mdb added 8/12/15 - enables admins to see all workflow results
    my $results = get_workflow_results($user_name, $id);

    $self->render(json => {
        id => int($id),
        status => lc($job_status->{status}),
        tasks => \@tasks,
        results => $results
    });
}

sub results { # legacy for Genome Export via HTTP
    my $self = shift;
    my $id = $self->stash('id');
    my $name = $self->stash('name');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $id);
    my $result_file = catfile($result_dir, $name);
    warn $result_file;

    unless (-r $result_file) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Download the file
    $self->res->headers->content_disposition("attachment; filename=$name;");
    $self->res->content->asset(Mojo::Asset::File->new(path => $result_file));
    $self->rendered(200);
}

1;

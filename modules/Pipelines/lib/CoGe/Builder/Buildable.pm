package CoGe::Builder::Buildable;

use Moose::Role;

use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);
use JSON qw(encode_json);
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;

use CoGe::Accessory::IRODS qw(irods_set_env irods_iput);
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage qw(get_workflow_paths get_workflow_results_file);

requires qw(build);

# Public attributes
has 'workflow'      => ( is => 'rw', isa  => 'CoGe::JEX::Workflow' );
has 'params'        => ( is => 'ro', required => 1 );
has 'site_url'      => ( is => 'rw' );
has 'page'          => ( is => 'rw' );
has 'db'            => ( is => 'ro', required => 1, isa  => 'CoGeX' );
has 'user'          => ( is => 'ro', required => 1, isa  => 'CoGeX::Result::User' );
has 'conf'          => ( is => 'ro', required => 1 );
has 'outputs'       => ( is => 'rw', default => sub { [] } );
has 'assets'        => ( is => 'rw', default => sub { [] } );
has 'errors'        => ( is => 'rw', default => sub { [] } );

# Private attributes
has 'staging_dir'   => ( is => 'rw');#, traits => ['Private'] );
has 'result_dir'    => ( is => 'rw');#, traits => ['Private'] );

my $previous_outputs = [];

sub get_site_url { # override this or use default in pre_build
    return '';
}

sub pre_build { # Default method, SynMap & SynMap3D override this
    my ($self, %params) = @_;

    # Initialize workflow -- NOTE: init => 1 means that a new workflow will be created right away
    $self->workflow( $params{jex}->create_workflow(name => $self->get_name, init => 1 ) );
    return unless ($self->workflow && $self->workflow->id);

    # Setup workflow paths
    my ($staging_dir, $result_dir) = get_workflow_paths(($self->user ? $self->user->name : 'public'), $self->workflow->id);
    make_path($staging_dir);
    make_path($result_dir);
    $self->staging_dir($staging_dir);
    $self->result_dir($result_dir);
    $self->workflow->logfile( catfile($result_dir, "debug.log") );

    # Set "page" and "site_url" attributes
    my $site_url = $self->get_site_url();
    if ($params{requester}) { # request is from internal web page - external API requests will not have a 'requester' field
        my $page = $params{requester}->{page}; # page name used for logging
        my $url  = $params{requester}->{url};  # used to set site_url with special query params
        $self->page($page) if $page;
        unless ($site_url) {
            $url = $page unless $url;
            $site_url = url_for($url, wid => $self->workflow->id);
        }
    }
    $self->site_url($site_url) if $site_url;
}

sub post_build {
    my $self = shift;

    # Capture outputs from legacy pipelines
    $self->sync_outputs();

    # Add task to add results to notebook
    if ($self->params->{notebook} || $self->params->{notebook_id}) {
        #TODO add_items_to_notebook_job and create_notebook_job and their respective scripts can be consolidated
        if ($self->params->{notebook_id}) { # Use existing notebook
            $self->add_task_chain_all(
                add_items_to_notebook( notebook_id => $self->params->{notebook_id} )
            );
        }
        else { # Create new notebook
            $self->add_task_chain_all(
                create_notebook( metadata => $metadata )
            );
        }
    }

    # Add task to send notification email
    if ( $self->params->{email} ) {
        $self->add_task_chain_all(
            $self->send_email(
                to => $self->params->{email},
                subject => 'CoGe Workflow Notification',
                body => 'Your workflow finished: ' . $self->get_name() .
                        ($self->site_url ? "\nLink: " . $self->site_url : '') .
                        "\n\nNote: you received this email because you submitted a job on " .
                        "CoGe (http://genomevolution.org) and selected the option to be notified."
            )
        );
    }

    # Add task to send notification to callback url
    if ( $self->params->{callback_url} ) {
        $self->add_task_chain_all(
            $self->curl_get( url => $self->params->{callback_url} )
        );
    }
}

# Add independent task
sub add_task {
    my ($self, $task) = @_;
    return unless $task;
    
    push @{$self->outputs}, @{$task->{outputs}};
    $previous_outputs = $task->{outputs};
    
    return $self->workflow->add_job($task);    
}

# Add independent tasks
#sub add_tasks {
#    my ($self, $tasks) = @_;
#    return unless $tasks;
#
#    $previous_outputs = [];
#    foreach my $task (@$tasks) {
#        push @{$self->outputs}, @{$task->{outputs}};
#        push @$previous_outputs, @{$task->{outputs}};
#        unless ($self->workflow->add_job($task)) {
#            return;
#        }
#    }
#
#    return 1;
#}

# Chain this task to the previous task assuming add_task*() routines were used
sub add_task_chain {
    my ($self, $task) = @_;

    push @{$task->{inputs}}, @{$previous_outputs} if $previous_outputs;
    return $self->add_task($task);
}

# Chain this task to all previous tasks
sub add_task_chain_all {
    my ($self, $task) = @_;
    
    # Chain this task to all previous tasks
    push @{$task->{inputs}}, @{$self->outputs} if $self->outputs;
    return $self->add_task($task);
}

sub add_output {
    my ($self, $output) = @_;
    push @{$self->outputs}, $output;
    $previous_outputs = [$output];
}

sub previous_output {
    my ($self, $index) = @_;
 
    $index = 0 unless defined $index;
    if ($previous_outputs && @$previous_outputs >= $index) {
        return $previous_outputs->[$index];    
    }
    
    return;
}

# Copies unknown ouputs from workflow into outputs array.  For legacy pipelines that don't use
# the add_task*() routines in this module.
sub sync_outputs {
    my $self = shift;

    my %seen;
    my @merged = grep( !$seen{$_}++, $self->workflow->get_outputs(), @{$self->outputs});
    $self->outputs(\@merged);
}

# Add a pipeline asset.  An asset is a file produced by the pipeline for use by downstream pipelines.
# The words input and output were purposely avoided.
sub add_asset {
    my ($self, $name, $value) = @_;
    return unless $name;
    push @{$self->assets}, [ $name, $value ];
}

sub get_asset {
    my ($self, $match_name) =  @_;
    
    foreach (@{$self->assets}) {
        my ($name, $value) = @{$_};
        if ($name eq $match_name) {
            return $value;
        }
    }
    
    return;
}

sub get_assets {
    my ($self, $match_name) =  @_;
    
    my @outputs;
    foreach (@{$self->assets}) {
        my ($name, $value) = @{$_};
        if ($name eq $match_name) {
            push @outputs, $value;
        }
    }
    
    return wantarray ? @outputs : \@outputs;
}


###############################################################################
# Task Library (replaces CommonTasks.pm)
###############################################################################

# Generate GFF file of genome annotations
sub create_gff {
    my ($self, %params) = @_;
    
    # Build argument list
    my $args = [
        ['-gid',     $params{gid},                0],
        ['-f',       $params{output_file},        1],
        ['-config',  $self->conf->{_CONFIG_PATH}, 0],
        ['-cds',     $params{cds} // 0,           0],
        ['-annos',   $params{annos} // 0,         0],
        ['-nu',      $params{nu} // 0,            0],
        ['-id_type', $params{id_type} // 0,       0],
        ['-upa',     $params{upa} // 0,           0]
    ];
    push @$args, ['-chr',     $params{chr},     0] if (defined $params{chr});
    push @$args, ['-add_chr', $params{add_chr}, 0] if (defined $params{add_chr});
    
    return {
        cmd         => catfile($self->conf->{SCRIPTDIR}, "coge_gff.pl"),
        args        => $args,
        outputs     => [ $params{output_file} ],
        description => "Generating GFF..."
    };
}

# Generate BED file of genome annotations
sub create_bed {
    my ($self, %params) = @_;

    return {
        cmd  => catfile($self->conf->{SCRIPTDIR}, "coge2bed.pl"), #FIXME script name makes no sense
        args => [
            ['-gid',    $params{gid},                0],
            ['-f',      $params{output_file},        1],
            ['-config', $self->conf->{_CONFIG_PATH}, 0]
        ],
        outputs => [ $params{output_file} ]
    };
}

# Generate TBL file of genome annotations
sub create_tbl {
    my ($self, %params) = @_;

    return {
        cmd     => catfile($self->conf->{SCRIPTDIR}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid',    $params{gid},                0],
            ['-f',      $params{output_file},        1],
            ["-config", $self->conf->{_CONFIG_PATH}, 0]
        ],
        outputs => [ $params{output_file} ]
    };
}

sub export_to_irods {
    my ($self, %params) = @_;

    irods_set_env(catfile($self->conf->{_HOME_PATH}, 'irodsEnv')); # mdb added 2/9/16 -- for hypnotoad, use www-data's irodsEnvFile
    my $cmd = irods_iput($params{src_file}, $params{dest_file}, { no_execute => 1, overwrite => ($params{overwrite} // 0) });

    my $done_file = catfile($self->staging_dir, basename($params{src_file})) . '.iput.done';

    return {
        cmd => qq[$cmd && touch $done_file],
        description => "Exporting file to IRODS " . $params{dest_file},
        args => [],
        inputs => [ $params{src_file} ],
        outputs => [ $done_file ]
    };
}

sub create_irods_imeta {
    my ($self, %params) = @_;
    
    my $done_file = catfile($self->staging_dir, basename($params{dest_file})) . '.imeta.done';
    
    my $cmd = catdir($self->conf->{SCRIPTDIR}, 'irods.pl') .
        " -cmd metadata" .
        " -dest " . $params{dest_file} . 
        " -metafile " . $params{metadata_file} .
        " -env " . catfile($self->conf->{_HOME_PATH}, 'irodsEnv') .
        " && touch $done_file";
    
    return {
        cmd => $cmd,
        description => "Setting IRODS metadata for " . $params{dest_file},
        args => [],
        inputs => [],
        outputs => [ $done_file ]
    };
}

sub untar {
    my ($self, %params) = @_;
    
    my $input_file = $params{input_file};
    my $output_path = $params{output_path};
    my $done_file = "$input_file.untarred";

    my $cmd = get_command_path('TAR');

    return {
        cmd => "mkdir -p $output_path && $cmd -xf $input_file --directory $output_path && touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file,
#            $input_file . '.done' # ensure file is done transferring
        ],
        outputs => [
            [$output_path, '1'],
            $done_file
        ],
        description => "Unarchiving " . basename($input_file) . "..."
    };
}

sub gunzip {
    my ($self, %params) = @_;
    my $input_file = $params{input_file};
    my $output_file = $input_file;
    $output_file =~ s/\.gz$//;

    my $cmd = get_command_path('GUNZIP');
    
    return {
        cmd => "$cmd -c $input_file > $output_file && touch $output_file.decompressed",
        script => undef,
        args => [],
        inputs => [
            $input_file,
#            $input_file . '.done' # ensure file is done transferring
        ],
        outputs => [
            $output_file,
            "$output_file.decompressed"
        ],
        description => "Decompressing " . basename($input_file) . "..."
    };
}

sub join_files {
    my ($self, %params) = @_;
    my $input_files = $params{input_files};
    my $input_dir   = $params{input_dir};
    my $output_file = $params{output_file};
    
    my @files;
    push @files, @$input_files   if $input_files;
    push @files, $input_dir.'/*' if $input_dir;
    
    my $cmd = "mkdir -p \$(dirname $output_file) && cat " . join(' ', @files) . ' > ' . $output_file;
    
    return {
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [
            @$input_files,
#            @$done_files
        ],
        outputs => [
            $output_file
        ],
        description => 'Joining files...'
    };
}

sub send_email {
    my ($self, %params) = @_;
    my $from    = 'CoGe Support <coge.genome@gmail.com>';
    my $to      = $params{to};
    my $subject = $params{subject};
    my $body    = $params{body};
    
    my $done_file = catfile($self->staging_dir, "send_email.done");
    
    my $args = [
        ['-from',      '"'.escape($from).'"',    0],
        ['-to',        '"'.escape($to).'"',      0],
        ['-subject',   '"'.escape($subject).'"', 0],
        ['-body',      '"'.escape($body).'"',    0],
    ];

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "send_email.pl"),
        script => undef,
        args => $args,
        inputs => [],
        outputs => [],
        description => "Sending notification email..."
    };
}

sub curl_get {
    my ($self, %params) = @_;
    my $url = $params{url};

    my $cmd = get_command_path('CURL');
    my $output_file = catfile($self->staging_dir, 'curl_output.log');

    return {
        cmd => "$cmd -s -o $output_file $url",
        script => undef,
        args => [],
        inputs => [],
        outputs => [ $output_file ],
        description => "Sending GET request to $url..."
    };
}

sub add_result {
    my ($self, %params) = @_;
    my $username = $self->user->name;
    my $wid      = $self->workflow->id;
    my $result   = encode_json($params{result});
    
    my $result_file = get_workflow_results_file($username, $wid);

    return {
        cmd  => catfile($self->conf->{SCRIPTDIR}, "add_workflow_result.pl"),
        args => [
            ['-user_name', $username,       0],
            ['-wid',       $wid,            0],
            ['-result',    "'".$result."'", 0] #TODO pass via temp file instead
        ],
        inputs  => [],
        outputs => [
#            $result_file,  # force this to run (for case of multiple results)
        ],
        description => "Adding workflow result..."
    };
}

sub add_items_to_notebook {
    my ($self, %params) = @_;
    my $notebook_id = $params{notebook_id};

    my $result_file = get_workflow_results_file($user->name, $wid);

    my $log_file = catfile($self->staging_dir, 'add_items_to_notebook', 'log.txt');

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'add_items_to_notebook.pl'),
        script => undef,
        args => [
            ['-uid', $self->user->id, 0],
            ['-wid', $self->workflow->id, 0],
            ['-notebook_id', $notebook_id, 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0],
            ['-log', $log_file, 0]
        ],
        inputs => [],
        outputs => [
            $result_file,
            $log_file
        ],
        description => "Adding experiment to notebook..."
    };
}

sub create_notebook {
    my ($self, %params) = @_;
    my $metadata = $params{metadata};
    my $annotations = $params{annotations}; # array ref

    my $result_file = get_workflow_results_file($self->user->name, $self->workflow->id);

    my $log_file = catfile($self->staging_dir, 'create_notebook', 'log.txt');

    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'create_notebook.pl'),
        script => undef,
        args => [
            ['-uid', $self->user->id, 0],
            ['-wid', $self->workflow->id, 0],
            ['-name', shell_quote($metadata->{name}), 0],
            ['-desc', shell_quote($metadata->{description}), 0],
            ['-type', 2, 0],
            ['-restricted', $metadata->{restricted}, 0],
            ['-annotations', qq{"$annotations_str"}, 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0],
            ['-log', $log_file, 0]
        ],
        inputs => [],
        outputs => [
            $result_file,
            $log_file
        ],
        description => "Creating notebook of results..."
    };
}

1;

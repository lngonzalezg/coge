package CoGe::Services::API::Download;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;

use File::Spec;
use File::Slurp;
use File::Basename;
use Data::Dumper;

use CoGe::Services::Auth qw(init);
use CoGe::Services::Error;
use CoGe::Core::Storage qw(get_download_path get_genome_cache_path get_experiment_cache_path get_workflow_log_file);
use CoGe::Accessory::Validate qw(valid_filename);

sub get {
    my $self = shift;
    my $filename = $self->param('filename');
    my $gid = $self->param('gid');
    my $eid = $self->param('eid');
    my $wid = $self->param('wid');
    my $attachment = $self->param('attachment') // 1;
    
    # Validate inputs
    unless ($gid || $eid || $wid) {
        print STDERR "CoGe::Services::Download invalid request\n";
        return $self->render(API_STATUS_MISSING_ID);
    }
    
    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Security pass 1 (F1): $filename was joined into a cache path with no containment
    # check, so ?filename=../../../../opt/apache2/coge/coge.conf read arbitrary files
    # (DB creds, JWT secret, /etc/passwd) with only a public-genome id as the gate.
    # Reduce to a single safe path component up front; every branch below uses it, so
    # catdir(<cache dir>, $filename) can no longer escape the cache directory.
    if (defined $filename && length $filename) {
        my $safe = valid_filename($filename);
        return $self->render(API_STATUS_BAD_REQUEST) unless defined $safe;
        $filename = $safe;
    }

    # Determine path to file
    my $file_path = '';
    if ($gid) { # genome FASTA/GFF export
        my $genome = $db->resultset('Genome')->find($gid);
        if ( $genome->restricted
            and ( not defined $user or not $user->has_access_to_genome($genome) ) )
        {
            print STDERR "CoGe::Services::Download access denied to genome $gid\n";
            return $self->render(API_STATUS_UNAUTHORIZED);
        }

        if ($filename) { # GFF file
            my $dl_path = get_genome_cache_path($gid);
            $file_path = File::Spec->catdir( $dl_path, $filename );
        }
        else { # FASTA file
            $file_path = get_genome_file($gid);
        }
    } 
    elsif ($eid) { # experiment tarball export
        my $exp = $db->resultset('Experiment')->find($eid);
        if ($exp->restricted and
            (not defined $user or not $user->has_access_to_experiment($exp))) {
            print STDERR "CoGe::Services::Download access denied to experiment $eid\n";
            return $self->render(API_STATUS_UNAUTHORIZED);
        }

        my $dl_path = get_experiment_cache_path($eid);
        $file_path = File::Spec->catdir($dl_path, $filename);
    }
    elsif ($wid) { # workflow debug log file
        unless ($user) {
            print STDERR "CoGe::Services::Download ERROR: not logged in\n";
            return $self->render(API_STATUS_UNAUTHORIZED);
        }
        # Default to the authenticated user. download_url_for has never emitted a "user"
        # param (its username lines are commented out, and would have spelled it
        # "username" anyway), so every workflow result URL ever written into .results
        # arrives without one. That left $user_name undef, and because
        # get_download_path/get_workflow_log_file build the path with catfile, the undef
        # collapsed to nothing instead of erroring: the lookup went to
        # .../downloads/jobs/<wid>/<file> with the username level missing entirely, and
        # 404'd. Non-admins never got that far -- undef ne their name tripped the check
        # below and returned 401.
        #
        # Defaulting here rather than adding the param to new URLs also repairs every
        # already-stored result link, and keeps usernames out of URLs. An admin fetching
        # another user's result can still pass user= explicitly.
        my $user_name = $self->param('user');
        $user_name = $user->name unless defined $user_name && length $user_name;
        if ($user_name ne $user->name && !$user->is_admin) {
            warn 'CoGe::Services::Download ERROR: different user';
            return $self->render(API_STATUS_UNAUTHORIZED);
        }
        # get workflow log by default
        unless ($filename) {
            $file_path = get_workflow_log_file($user_name, $wid);
            unless (-r $file_path) {
                print STDERR "CoGe::Services::Download ERROR: workflow log not found for wid $wid in $file_path\n";
                return $self->render(API_STATUS_NOTFOUND);
            }
            $filename = "workflow_$wid.log";
        }
        else {
            my $dl_path = get_download_path('jobs', $user_name, $wid);
            $file_path = File::Spec->catdir($dl_path, $filename);
        }
    }
    
    warn "CoGe::Services::Download file=$file_path";
    return $self->render(API_STATUS_NOTFOUND) unless $file_path && -e $file_path;

    # Send file
    $self->res->headers->content_disposition("attachment; filename=$filename;") if $attachment; # tell browser to download file
    my $content;
    eval {
        $content = read_file($file_path) if -r $file_path;
    };

    $self->render(data => $content);
}

1;

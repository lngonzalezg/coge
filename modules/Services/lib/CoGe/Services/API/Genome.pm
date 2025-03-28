package CoGe::Services::API::Genome;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(decode_json);
use Data::Dumper;

use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Favorites;
use CoGe::Core::Genome qw(genomecmp search_genomes);
use CoGe::Core::Metadata qw( create_annotation delete_annotation get_annotation get_annotations );
use CoGe::Core::Storage qw(get_genome_file get_genome_seq);
use CoGe::Services::API::Job;
use CoGe::Services::Auth qw(init);
use CoGe::Services::Error;
use CoGeDBI qw(get_feature_counts get_feature_count_summary get_features);

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $fast = $self->param('fast');
    $fast = (defined $fast && ($fast eq '1' || $fast eq 'true')); # default to false
    my $sort = $self->param('sort');
    $sort = (defined $sort && ($sort eq '1' || $sort eq 'true')); # default to false

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(API_STATUS_SEARCHTERM);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);
    
    # Search notebooks and filter based on user permissions
    my $filtered = search_genomes(db => $db, search_term => $search_term, user => $user, sort => $sort);

    # Get user's favorites
    my $favorites = CoGe::Core::Favorites->new(user => $user);

    # Format response
    my @result;
    if ($fast) {
        @result = map {
          {
            id => int($_->id),
            info => $_->info,
            certified  => $_->certified ? Mojo::JSON->true : Mojo::JSON->false,
            restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
            deleted    => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
            favorited  => $favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false
          }
        } @$filtered;
    }
    else {
        @result = map {
          {
            id          => int($_->id),
            name        => $_->name,
            description => $_->description,
            link        => $_->link,
            version     => $_->version,
            info        => $_->info,
            organism_id   => int($_->organism->id),
            sequence_type => {
                id          => ($_->type ? $_->type->id : 0),
                name        => ($_->type ? $_->type->name : ''),
                description => ($_->type ? $_->type->description : ''),
            },
            restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
            certified  => $_->certified ? Mojo::JSON->true : Mojo::JSON->false,
            deleted    => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
            chromosome_count => int($_->chromosome_count),
            organism => {
                id          => int($_->organism->id),
                name        => $_->organism->name,
                description => $_->organism->description
            }
          }
        } @$filtered;
    }
    
    $self->render(json => { genomes => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $genome = $self->_get_genome($id, 0, $db, $user);
    return unless $genome;

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $genome->annotations;
    
    # Build chromosome list
    my $chromosomes = $genome->chromosomes_all;
    my $feature_counts = get_feature_counts($db->storage->dbh, $genome->id);
    foreach (@$chromosomes) {
        my $name = $_->{name};
        $_->{gene_count} = int($feature_counts->{$name}{1}{count} // 0);
        $_->{CDS_count} = int($feature_counts->{$name}{3}{count} // 0);
    }

    # Build feature types
    my $ftypes = get_feature_count_summary($db->storage->dbh, $genome->id);
    foreach (@$ftypes) {
        $_->{count}   = int($_->{count});
        $_->{type_id} = int($_->{type_id});
    }

    # Generate response
    $self->render(json => {
        id => int($genome->id),
        name => $genome->name,
        description => $genome->description,
        link => $genome->link,
        version => $genome->version,
        restricted => $genome->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        certified => $genome->certified ? Mojo::JSON->true : Mojo::JSON->false,
        organism => {
            id => int($genome->organism->id),
            name => $genome->organism->name,
            description => $genome->organism->description
        },
        sequence_type => {
            id => int($genome->type ? $genome->type->id : 0),
            name => ($genome->type ? $genome->type->name : ''),
            description => ($genome->type ? $genome->type->description : ''),
        },
        chromosome_count => int($genome->chromosome_count),
        chromosomes => $chromosomes,
        feature_types => $ftypes,
        experiments => [ map { int($_->id) } $genome->experiments ],
        additional_metadata => \@metadata
    });
}

sub fetch_annotations {
    my $self = shift;
    my $id = int($self->stash('id'));
    my ($db) = CoGe::Services::Auth::init($self);
    #TODO add error checking on ID param

    $self->render(json => get_annotations($id, 'Genome', $db, 1));
}

sub fetch_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));
    #TODO add error checking on ID params

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $genome = $self->_get_genome($id, 0, $db, $user);
    return unless $genome;

    my $annotation = get_annotation($aid, 'Genome', $db);
    $self->render(json => $annotation) if $annotation;
}

sub sequence {
    my $self   = shift;
    my $gid    = $self->stash('id');
    return unless $gid;
    my $chr    = $self->stash('chr');    # optional
    my $start  = $self->param('start');  # optional
    my $stop   = $self->param('stop') || $self->param('end'); # optional
    my $strand = $self->param('strand'); # optional
    print STDERR "API::Genome::sequence gid=$gid ",
        (defined $chr ? "chr=$chr " : ''),
        (defined $start ? "start=$start " : ''),
        (defined $stop ? "stop=$stop " : ''), "\n";

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $genome = $self->_get_genome($gid, 0, $db, $user);
    return unless $genome;

    # Force browser to download whole genome as attachment
    my $format;
    my $genome_name = sanitize_name($genome->organism->name);
    $genome_name = 'genome_'.$gid unless $genome_name;
    if ( (!defined($chr) || $chr eq '') ) {
        $self->res->headers->content_disposition("attachment; filename=$genome_name.faa;");
    }
    elsif (defined($chr) && !defined($start) && !defined($stop)) {
        $genome_name .= '_' . $chr;
        $stop = $genome->get_chromosome_length($chr);
        $format = 'fasta';
        $self->res->headers->content_disposition("attachment; filename=$genome_name.faa;");
    }

    # Get sequence from file
    unless (defined $chr and $chr ne '') {
        $self->res->headers->content_type('text/plain');
        $self->res->content->asset(Mojo::Asset::File->new(path => get_genome_file($gid)));
        $self->rendered(200);
    }
    else {
        $self->render(text => get_genome_seq(
            gid   => $gid,
            chr   => $chr,
            start => $start,
            stop  => $stop,
            strand => $strand,
            format => $format
        ));
    }
}

sub features {
    my $self = shift;
    my $gid    = $self->stash('id');
    return unless $gid;
    my $type    = $self->stash('type'); # optional
    print STDERR "API::Genome::features gid=$gid ", (defined $type ? "type=$type " : '');

    # Connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    unless (defined $genome) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Check permissions
    unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    my $type_id;
    if ($type) {
        my $feature_type = $db->resultset('FeatureType')->find({ name => $type });
        unless ($feature_type) {
            $self->render(API_STATUS_BAD_REQUEST("Feature type not present"));
            return;
        }
        $type_id = $feature_type->id;
    }

    my $features = get_features($db->storage->dbh, $gid, undef, $type_id, 1);
    my @features = map {
        {   id         => int($_->{id}),
            type       => $_->{type_name},
            chromosome => $_->{chromosome},
            start      => int($_->{start}),
            stop       => int($_->{stop}),
            strand     => int($_->{strand}),
            name       => $_->{name}
        }
    } @$features;

    # Generate response
    $self->render(json => {
        id => int($genome->id),
        features => \@features
    });
}

sub export {
    my $self = shift;
    my $gid  = $self->stash('id');
    return unless $gid;
    my $data = $self->req->json;
    $data->{gid} = $gid;

    # Alias to Job Submit -- is there a better way to do this using Mojolicious routing?
    my $request = {
        type => 'export_genome',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

sub add {
    my $self = shift;
    my $data = $self->req->body; #$self->req->json; # mdb replaced 11/30/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    unless ($data) {
        $self->render(API_STATUS_MISSING_BODY);
        return;
    }
    $data = decode_json($data);
    #print STDERR "CoGe::Services::Data::Genome::add\n", Dumper $data, "\n";

    # Valid data items # TODO move into request validation
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(API_STATUS_MISSING_DATA);
        return;
    }
    
    # Alias to Job Submit -- is there a better way to do this using Mojolicious routing?
    my $request = {
        type => 'load_genome',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

sub add_annotation {
    my $self = shift;
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    #TODO add error checking on params
    create_annotation(
        conf => $conf,
        db => $db,
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => int($self->stash('id')),
        target_type => 'genome',
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
    my $error = CoGe::Core::Metadata::delete_annotation($aid, $id, 'Genome', $db, $user);
    if ($error) {
        $self->render(API_STATUS_BAD_REQUEST($error));
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub update {
	my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $genome = $self->_get_genome($id, 1, $db, $user);
    return unless $genome;

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$genome->update($data->{metadata});
	$self->render(json => { success => Mojo::JSON->true });
}

sub update_annotation {
    my $self = shift;
    my ($db, $user) = CoGe::Services::Auth::init($self);
    #TODO add error checking on params
    CoGe::Core::Metadata::update_annotation(
        annotation_id => int($self->stash('aid')),
        db => $db,
        delete_bisque_image => $self->param('delete_bisque_image'),
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => $self->stash('id'),
        target_type => 'genome',
        text => $self->param('annotation'),
        type_name => $self->param('type_name'),
        user => $user
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub _get_genome {
    my ($self, $id, $own_or_edit, $db, $user) = @_;
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }
    my $genome = $db->resultset("Genome")->find($id);
    unless (defined $genome) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }
    if ($own_or_edit) {
        unless ($user) {
            $self->render(API_STATUS_CUSTOM(401, "User not logged in"));
            return;
        }
        return $genome if $user->is_admin;
        unless ($user->is_owner_editor(genome => $id)) {
            $self->render(API_STATUS_UNAUTHORIZED);
            return;
        }
    }
    unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    return $genome;
}

1;

package CoGe::Builder::Tools::SynMap3D;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::JEX::Jex;
use CoGe::Builder::Tools::SynMap qw( defaults gen_org_name );
use CoGe::Accessory::Validate qw(valid_id valid_filename);
use CoGe::Exception::Generic;
use File::Spec::Functions;

sub pre_build { # override superclass method
	my $self = shift;

	# Connect to JEX
    my $jex = CoGe::JEX::Jex->new( host => $self->conf->{JOBSERVER}, port => $self->conf->{JOBPORT} );
    unless ($jex) {
        CoGe::Exception::Generic->throw(message => "Couldn't connect to JEX");
    }
    $self->jex($jex);

	# Initialize workflow -- NOTE: init => 0 means that a previous identical workflow will be reused when submitted
    $self->workflow( $jex->create_workflow(name => $self->get_name, init => 0 ) );
    return unless $self->workflow;

	# Set site_url attribute
	my %opts = ( %{ defaults() }, %{ $self->params } );
	$self->site_url( $opts{tinylink} || get_query_link( $self->conf, $self->db, %opts ) );
	my $requester = $self->request->requester;
    if ($requester) { # request is from internal web page - external API requests will not have a 'requester' field
        my $page = $requester->{page}; # page name used for logging
        $self->page($page) if $page;
    }

	return;
}

sub build {
	my $self = shift;
	# §4.3.4 all merge params below are spliced into the synmerge_3.py shell
	# command line; validate before use so nothing can inject shell tokens.
	_sanitize_3d_opts($self->params);
	my $xid = $self->params->{genome_id1};
	my $yid = $self->params->{genome_id2};
	my $zid = $self->params->{genome_id3};
	my ($dir1, $dir2) = sort($xid, $yid);
	my ($dir3, $dir4) = sort($xid, $zid);
	my ($dir5, $dir6) = sort($yid, $zid);

	my $SYN3DIR   = $self->conf->{SYN3DIR};
	my $SCRIPTDIR = catdir( $self->conf->{SCRIPTDIR}, 'synmap' );
	my $PYTHON    = $self->conf->{PYTHON} // 'python';
	my $MERGER    = 'nice ' . $PYTHON . ' ' . catfile($SCRIPTDIR, 'synmerge_3.py');

	#########################################################################
	# Add tasks for all pairwise SynMap comparisons
	#########################################################################

    my $synmap = CoGe::Builder::Tools::SynMap->new($self);
    $synmap->build();
    $self->add_to_all($synmap);

	#########################################################################
	# Add task to merge results
	#########################################################################

	# Input datafiles. TODO: Make this path specification based on the '$final_dagchainer_file' from SynMap.pm
	# AKB (2016-8-15) changed to reflect new naming conventions.
	#my $dot_xy_path = catfile($self->conf->{DIAGSDIR}, $dir1, $dir2, $dir1 . '_' . $dir2 . "_synteny.json");
    #my $dot_xz_path = catfile($self->conf->{DIAGSDIR}, $dir3, $dir4, $dir3 . '_' . $dir4 . "_synteny.json");
	#my $dot_yz_path = catfile($self->conf->{DIAGSDIR}, $dir5, $dir6, $dir5 . '_' . $dir6 . "_synteny.json");

	my $opts_name = '.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.Dm0.ma1.gcoords.dotplot_dots_synteny.json';  # check that this is right. Maybe there needs to be '.gcoords' after 'ma1'? or no '.Dm0.ma1'? Depends on state of final_dagchainer_file at point of assignment...
	my $dot_xy_path = catfile( $self->conf->{DIAGSDIR}, $dir1, $dir2, $dir1 . '_' . $dir2 . $opts_name );
	my $dot_xz_path = catfile( $self->conf->{DIAGSDIR}, $dir3, $dir4, $dir3 . '_' . $dir4 . $opts_name );
	my $dot_yz_path = catfile( $self->conf->{DIAGSDIR}, $dir5, $dir6, $dir5 . '_' . $dir6 . $opts_name );

	# Options.
	my $sort = $self->params->{sort};
    my $min_length = $self->params->{min_length};
    my $min_synteny = $self->params->{min_synteny};
    my $cluster = $self->params->{cluster};
    my $c_eps;
    my $c_min;
    if ($cluster eq 'true') {
        $c_eps = $self->params->{c_eps};
        $c_min = $self->params->{c_min};
    }
    my $ratio = $self->params->{ratio};
    my $r_by;
    my $r_min;
    my $r_max;
    if ($ratio ne 'false') {
        $r_by = $self->params->{r_by};
        $r_min = $self->params->{r_min};
        $r_max = $self->params->{r_max};
    }
    my $graph_out = $self->params->{graph_out};
    my $log_out = $self->params->{log_out};
    my $data_out = $self->params->{download};

	# Build command line arguments.
	my $merge_ids = ' -xid ' . $xid . ' -yid ' . $yid . ' -zid ' . $zid;
    my $merge_ins = ' -i1 ' . $dot_xy_path . ' -i2 ' . $dot_xz_path . ' -i3 ' . $dot_yz_path;
    my $merge_otp = ' -o ' . $SYN3DIR;
    my $merge_opt = ' -S ' . $sort . ' -ml ' . $min_length . ' -ms ' . $min_synteny;
    if ($cluster eq 'true') {
        $merge_opt .= ' -C -C_eps ' . $c_eps . ' -C_ms ' . $c_min;
    }
    if ($ratio ne 'false') {
        $merge_opt .= ' -R ' . $ratio . ' -Rby ' . $r_by . ' -Rmin ' . $r_min . ' -Rmax ' . $r_max;
    }

	$self->add_to_all({
		description => "Identifying common points & building graph object",
        cmd => $MERGER . $merge_ids . $merge_ins . $merge_opt . $merge_otp,
        inputs => [
			$dot_xy_path,
			$dot_xz_path,
			$dot_yz_path
		],
        outputs => [
			catfile($SYN3DIR, $graph_out),
			catfile($SYN3DIR, $log_out),
			catfile($SYN3DIR, $data_out)
		]
    });
}

sub get_name {
    my $self = shift;
    
    my $description;
    foreach my $key (keys %{$self->params}) {
        next unless $key =~ /^genome_id/;
        
        my ($genome_id) = $key =~ /(\d+)$/;
        my ( $org_name ) = gen_org_name(
            db        => $self->db,
            genome_id => $genome_id
        );
        
        $description .= ($description ? 'v. ' : '') . $org_name . ' ';
    }
    
    $description .= 'Ks' if $self->params->{ks_type};
    
	return $description;
}

# §4.3.4 Validate every SynMap3D merge parameter that is spliced into the
# synmerge_3.py command line or into an output filename. Fails closed (throws)
# on any value carrying shell metacharacters.
sub _sanitize_3d_opts {
	my $params = shift;
	my $float_re = qr/^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$/;

	for my $k (qw(genome_id1 genome_id2 genome_id3)) {
		CoGe::Exception::Generic->throw(message => "Invalid $k")
			unless defined valid_id( $params->{$k} );
	}
	# Numeric merge options ( -ml, -ms, -C_eps, -C_ms, -R, -Rmin, -Rmax ).
	for my $k (qw(min_length min_synteny c_eps c_min ratio r_min r_max)) {
		next unless defined $params->{$k} && length $params->{$k};
		next if $k eq 'ratio' && $params->{$k} eq 'false';
		CoGe::Exception::Generic->throw(message => "Invalid $k")
			unless $params->{$k} =~ $float_re;
	}
	# Token options ( -S sort key, -Rby ratio-by column ).
	for my $k (qw(sort r_by)) {
		next unless defined $params->{$k} && length $params->{$k};
		CoGe::Exception::Generic->throw(message => "Invalid $k")
			unless $params->{$k} =~ /^[\w.\-]+$/;
	}
	# Output filenames ( become catfile() args -> JEX outputs ). The raw value
	# is used as-is, so require it to already be a clean single-component name
	# (reject path separators and '..' traversal, not just sanitise them).
	for my $k (qw(graph_out log_out download)) {
		next unless defined $params->{$k} && length $params->{$k};
		my $clean = valid_filename( $params->{$k} );
		CoGe::Exception::Generic->throw(message => "Invalid $k")
			unless defined $clean && $clean eq $params->{$k};
	}
	return;
}

__PACKAGE__->meta->make_immutable;

1;

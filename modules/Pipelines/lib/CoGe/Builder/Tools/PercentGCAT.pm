package CoGe::Builder::Tools::PercentGCAT;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Validate qw(valid_id);
use CoGe::Exception::Generic;
use File::Spec::Functions;

sub get_name {
	my $self = shift;
    return 'PercentGCAT | ' . $self->params->{'gid'} . ' | ' . $self->params->{'chr'};
}

sub build {
	my $self = shift;

    # §4.3.4 validate inputs that reach shell command strings below.
    my $gid = valid_id( $self->params->{'gid'} );
    CoGe::Exception::Generic->throw(message => "Invalid gid")
        unless defined $gid;
    my $chr = $self->params->{'chr'};
    CoGe::Exception::Generic->throw(message => "Invalid chr")
        unless defined $chr && $chr =~ /^[\w.\-]+$/;
    my $wsize = $self->params->{'wsize'}; # sliding window size
    my $wstep = $self->params->{'wstep'}; # sliding window spacing
    CoGe::Exception::Generic->throw(message => "Invalid wsize")
        unless defined $wsize && $wsize =~ /^\d+$/;
    CoGe::Exception::Generic->throw(message => "Invalid wstep")
        unless defined $wstep && $wstep =~ /^\d+$/;
    my $dir = catfile($self->conf->{SECTEMPDIR}, "downloads/genome", $gid);

    my $fasta = catfile($dir, $gid . '_' . $chr . '.faa');
    $self->add({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'generate_chr_fasta.pl'),
        args        => [[ 'gid', $gid, 0 ], [ 'chr', $chr, 0 ]],
        outputs     => [$fasta],
        description => "Generating chromosome sequence",        
    });

    my $filename = $gid . '_' . $chr . '_' . $wsize . '_' . $wstep . '_out.txt';
    my $output = catfile($dir, $filename);
    $self->add({
        cmd         => catfile($self->conf->{SCRIPTDIR}, 'percent_gc_at.py') . ' ' . $fasta . ' ' . $wsize . ' ' . $wstep,
        inputs      => [$fasta],
        outputs     => [$output],
        description => "Generating nucleotide sliding window percentages",        
    });

    if ($self->params->{'irods'}) {
        my $irods_base = irods_get_base_path($self->user->name);
        $self->add(
            $self->export_to_irods(
                src_file  => $output,
                dest_file => catfile($irods_base, $filename)
            )
        );
    }
}

__PACKAGE__->meta->make_immutable;

1;

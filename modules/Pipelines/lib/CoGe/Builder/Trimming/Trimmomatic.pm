package CoGe::Builder::Trimming::Trimmomatic;

use Moose;
extends 'CoGe::Builder::Trimming::Trimmer';

use Data::Dumper;
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;
    my ($fastq1, $fastq2) =  @{$opts{data_files}}; # fastq2 is undef for single-ended
    unless ($fastq1 && @$fastq1) {
        CoGe::Exception::Generic->throw(message => 'Missing fastq');
    }

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    if ($read_type eq 'single') { # single-ended
        # Create trimmomatic task for each file
        foreach (@$fastq1) {
            $self->add(
                $self->trimmomatic($_),
                qq[$_.done] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, $self->previous_output;
        }
    }
    else { # paired-end
        # Create trimmomatic task for each file pair
        for (my $i = 0;  $i < @$fastq1;  $i++) {
            my ($f1, $f2) = ($fastq1->[$i], $fastq2->[$i]);
            $self->add(
                $self->trimmomatic([ $fastq1->[$i], $fastq2->[$i] ]),
                [ qq{$f1.done}, qq{$f2.done} ] # done file dependencies are created in Extractor
            );
            push @{$self->fastq}, @{$self->previous_outputs};
        }
    }
}

sub trimmomatic {
    my $self = shift;
    my $fastq = shift;
    $fastq = [ $fastq ] unless (ref($fastq) eq 'ARRAY');

    my $trimming_params = $self->params->{trimming_params} // {};
    my $read_params = $self->params->{read_params} // {};
    my $encoding = $read_params->{encoding} // 33;
    my $read_type = $read_params->{read_type} // 'single';

    my $cmd = get_command_path('TRIMMOMATIC');
    $cmd = 'nice java -jar ' . $cmd;

    # Initialize args properly
    my $args = [];

    if ($read_type eq 'paired') {
        push @$args, [ 'PE', '', 0 ];
	}
    else {
        # Single-ended mode
        push @$args, [ 'SE', '', 0 ];
	}
    # Add phred encoding first (correct order per documentation)
    push @$args, [($encoding == 64 ? '-phred64' : '-phred33'), '', 0];
    push @$args, [ '-threads', 8, 0 ];
    push @$args, [ '-trimlog', 'trimmomatic.log', 0 ];

    my @outputs;

    if ($read_type eq 'paired') {
        push @$args, map { [ '', $_, 1 ] } @$fastq;
        
        # Create proper output files for paired-end mode
        my $base = to_filename($fastq->[0]);
        my $ext = '.fastq' . to_compressed_ext($fastq->[0]);
        
        my $paired_out1 = catfile($self->staging_dir, $base . '_1P' . $ext);
        my $unpaired_out1 = catfile($self->staging_dir, $base . '_1U' . $ext);
        my $paired_out2 = catfile($self->staging_dir, $base . '_2P' . $ext);
        my $unpaired_out2 = catfile($self->staging_dir, $base . '_2U' . $ext);
        
        push @$args, [ '', $paired_out1, 1 ];
        push @$args, [ '', $unpaired_out1, 1 ];
        push @$args, [ '', $paired_out2, 1 ];
        push @$args, [ '', $unpaired_out2, 1 ];
        
        # Update outputs to what Trimmomatic will actually produce
        @outputs = ($paired_out1, $paired_out2, $unpaired_out1, $unpaired_out2);
    }
    else {
        # Single-ended mode
        my $output = catfile($self->staging_dir, 
                             basename(remove_fastq_ext($fastq->[0]) . 
                             '.trimmed.fastq' . to_compressed_ext($fastq->[0])));
        
        push @$args, [ '', $fastq->[0], 1 ];
        push @$args, [ '', $output, 1 ];
        
        @outputs = ($output);
    }

    # Add all trimming steps
    foreach ('ILLUMINACLIP', 'SLIDINGWINDOW', 'MAXINFO', 'LEADING', 
            'TRAILING', 'CROP', 'HEADCROP', 'MINLEN', 'AVGQUAL') {
        my $value = $trimming_params->{$_};
        push @$args, ["$_:$value", '', 0] if ($value);
    }

    return {
        cmd => $cmd,
        args => $args,
        inputs => [ @$fastq ],
        outputs => \@outputs,
        description => 'Trimming (Trimmomatic) ' . fastq_description($fastq, $read_type)
    };
}
1;

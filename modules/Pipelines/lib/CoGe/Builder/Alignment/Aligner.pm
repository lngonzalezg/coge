package CoGe::Builder::Alignment::Aligner;

use Moose;
extends 'CoGe::Builder::Buildable';

use Switch;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use Clone qw(clone);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);

use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata qw(to_annotations);
use CoGe::Builder::Trimming::Trimmer;
use CoGe::Builder::Alignment::HISAT2;
use CoGe::Builder::Alignment::GSNAP;
use CoGe::Builder::Alignment::BWA;
use CoGe::Builder::Alignment::Bowtie;
use CoGe::Builder::Alignment::Tophat;
use CoGe::Builder::Alignment::Bismark;
use CoGe::Builder::Alignment::BWAmeth;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

# Settings
has NUM_CPUS       => (is => 'ro', isa => 'Int', default => 16); # number of CPUs to use for alignment tasks
has VALIDATE_FASTQ => (is => 'ro', isa => 'Int', default => 0);  # flag to indicate FASTQ input files should be checked for correct format

# Outputs
has index   => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # index files
has bam     => (is => 'rw', isa => 'ArrayRef', default => sub { [] }); # processed bam files (sorted and indexed)

sub build {
    my $self = shift;
    my %opts = @_;
    my $fastq = $opts{data_files}; # array ref of FASTQ files
    unless ($fastq && @$fastq) {
        CoGe::Exception::MissingField->throw(message => 'Missing fastq');
    }

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    my $gid = $self->request->genome->id;

# mdb removed 11/6/15 COGE-673
#    # Check multiple files (if more than one file then all should be FASTQ)
#    my $numFastq = 0;
#    foreach (@$input_files) {
#        $numFastq++ if (is_fastq_file($_));
#    }
#    if ($numFastq > 0 and $numFastq != @$input_files) {
#        my $error = 'Unsupported combination of file types';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }
#    if ($numFastq == 0 and @$input_files > 1) {
#        my $error = 'Too many files';
#        print STDERR 'CoGe::Builder::Common::Alignment ERROR: ', $error, "\n";
#        return { error => $error };
#    }

    #
    # Build workflow
    #

    # Validate the input files
    if ($self->VALIDATE_FASTQ) {
        foreach my $input_file (@$fastq) {
            $self->add_to_previous( # previous task is file trasfer or decompression
                $self->validate_fastq($input_file)
            );
        }
    }

    # Reheader the fasta file - needs to happen before trimming for BBDuk
    my ($reheader_fasta) = $self->add(
        $self->reheader_fasta($gid)
    );

    # Index the fasta file
    $self->add_to_previous(
        $self->index_fasta($reheader_fasta)
    );

    # Trim the fastq input files #TODO consider moving this and validation into Load/Experiment
    my @trimmed;
    if ($self->params->{trimming_params}) {
        my $trimmer = CoGe::Builder::Trimming::Trimmer->new($self);
        $trimmer->build(data_files => $fastq);
        $self->add_to_all($trimmer);
        @trimmed = @{$trimmer->fastq};
    }
    else { # no trimming
        @trimmed = @$fastq;
    }

    my $aligner;
    switch( $self->_aligner ) {
        case 'hisat2'  { $aligner = CoGe::Builder::Alignment::HISAT2->new($self)  }
        case 'bowtie2' { $aligner = CoGe::Builder::Alignment::Bowtie->new($self)  }
        case 'tophat'  { $aligner = CoGe::Builder::Alignment::Tophat->new($self)  }
        case 'bismark' { $aligner = CoGe::Builder::Alignment::Bismark->new($self) }
        case 'bwameth' { $aligner = CoGe::Builder::Alignment::BWAmeth->new($self) }
        case 'bwa'     { $aligner = CoGe::Builder::Alignment::BWA->new($self)     }
        case 'gsnap'   { $aligner = CoGe::Builder::Alignment::GSNAP->new($self)   }
        default {
            CoGe::Exception::Generic->throw(message => 'Invalid aligner');
        }
    }
    $aligner->build(fasta_file => $reheader_fasta, data_files => \@trimmed);
    $self->add_to_all($aligner);
    push @{$self->bam}, @{$aligner->bam};

    # Process and load each alignment output
    foreach my $bam_file (@{$self->bam}) {
        # commented out 4/12/2017 at Jeff Grover's suggestion. bismark is only for methylation
        # Sort BAM file(s) -- bismark only (because methylation analysis requires unsorted BAM file), other aligner modules perform sort internally
        # if ($self->_aligner eq 'bismark') { #) && !$self->params->{methylation_params}) {
        #     ($bam_file) = $self->add(
        #         $self->sort_bam($bam_file)
        #     );
        # }

        # Index BAM file -- requires sorted BAM
        unless ($self->_aligner eq 'bismark') {
            $self->add(
                $self->index_bam($bam_file)
            );
        }

        # Load alignment
        if ($self->params->{alignment_params} && $self->params->{alignment_params}->{load_bam}) {
            # Get custom metadata to add to experiment
            my $annotations = $self->generate_additional_metadata();

            # Add bam filename to experiment name for ChIP-seq pipeline
            my $md = clone($metadata);
            if (@{$self->bam} > 1) {
                $md->{name} .= ' (' . to_filename_base($bam_file) . ')';
            }

            $self->add(
                $self->load_bam(
                    metadata    => $md,
                    annotations => $annotations,
                    bam_file    => $bam_file
                )
            );
        }
    }
}

sub validate_fastq {
    my $self = shift;
    my $fastq = shift;

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "validate_fastq.pl"),
        args => [
            ["", $fastq, 0] # mdb changed 3/1/16 from 1 to 0, COGE-707
        ],
        inputs => [
            $fastq
        ],
        outputs => [
            "$fastq.validated"
        ],
        description => "Validating " . basename($fastq)
    };
}

sub bowtie2_index { # shared between Bowtie and Tophat
    my $self = shift;
    my $fasta = shift;

    my $gid = $self->request->genome->id;
    my $cache_dir = catdir(get_genome_cache_path($gid), "bowtie_index");
	make_path($cache_dir) unless (-d $cache_dir);
    my $name = catfile($cache_dir, 'genome.reheader');

    my $cmd = get_command_path('BOWTIE_BUILD', 'bowtie2-build');
    $cmd = $cmd . " --large-index  ";

    $self->index([
        $name . ".1.bt2l",
        $name . ".2.bt2l",
        $name . ".3.bt2l",
        $name . ".4.bt2l",
        $name . ".rev.1.bt2l",
        $name . ".rev.2.bt2l"
    ]);

    return {
        cmd => $cmd,
        args => [
            ["", $fasta, 1],
            ["", $name, 0],
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            @{$self->index}
        ],
        description => "Indexing genome sequence with Bowtie"
    };
}

sub generate_additional_metadata { #TODO redo arg capture in a more automated fashion
    my $self = shift;
    my $read_params      = $self->params->{read_params};
    my $trimming_params  = $self->params->{trimming_params};
    my $alignment_params = $self->params->{alignment_params};

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php?title=LoadExp%2B||note|Generated by CoGe's NGS Analysis Pipeline};

    switch( lc($trimming_params->{trimmer}) ) {
        case 'cutadapt' {
            push @annotations, 'note|cutadapt '.join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '-m'));
        }
        case 'trimgalore' {
            push @annotations, 'note|trimgalore '.join(' ', map { $_.' '.$trimming_params->{$_} } ('-q', '--length', '-a'));
        }
        case 'trimmomatic' {
            push @annotations, 'note|trimmomatic '.join(' ', map { $_.':'.$trimming_params->{$_} } grep { $trimming_params->{$_} } ('ILLUMINACLIP', 'SLIDINGWINDOW', 'MAXINFO', 'LEADING', 'TRAILING', 'CROP', 'HEADCROP', 'MINLEN'));
        }
        case 'bbduk' {
            my @args = map { $_.'='.$trimming_params->{$_} } ('k', 'mink', 'hdist', 'tpe', 'tbo', 'qtrim', 'trimq', 'minlength');

            my $adapters = $trimming_params->{'adapters'};
            if ($adapters && $adapters ne 'none') {
                if ($adapters eq 'r' || $adapters eq 'both') {
                    push @args, 'rref=adapters.fa';
                }
                if ($adapters eq 'l' || $adapters eq 'both') {
                    push @args, 'lref=adapters.fa';
                }
            }

            push @annotations, 'note|bbduk '.join(' ', @args);
        }
    }

    switch( $self->_aligner() ) {
        case 'hisat2'  {
            push @annotations, qq{note|hisat2_build};
            my $params = ($read_params->{encoding} eq '64' ? '--phred64' : '--phred33');
            push @annotations, 'note|hisat2 ' . $params;
        }
        case 'bowtie2' {
            my $rg = $alignment_params->{'--rg-id'};
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|bowtie2 ' . $alignment_params->{'presets'} . ($rg ? " --rg-id $rg" : '');
        }
        case 'tophat'  {
            push @annotations, qq{note|bowtie2_build};
            push @annotations, 'note|tophat ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-g'));
        }
        case 'bismark' {
            push @annotations, qq{note|bismark_genome_preparation};
            push @annotations, 'note|bismark ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-L'));
        }
        case 'bwameth' {
            push @annotations, qq{note|bwameth index};
            push @annotations, 'note|bwameth (default options)';
        }
        case 'bwa'     {
            my $M = $alignment_params->{'-M'};
            my $R = $alignment_params->{'-R'};
            my $args_str = ($M ? '-M' : '') . ($R ? " -R $R" : '');
            push @annotations, qq{note|bwa index};
            push @annotations, 'note|bwa mem ' . ($args_str ? $args_str : ' (default options)');
        }
        case 'gsnap'   {
            push @annotations, qq{note|gmap_build};
            push @annotations, 'note|gsnap ' . join(' ', map { $_.' '.$alignment_params->{$_} } ('-N', '-n', '-Q', '--gap-mode', '--nofails'));
        }
    }

    return \@annotations;
}

sub _aligner {
    my $self = shift;
    my $alignment_params = $self->params->{alignment_params};

    if ($alignment_params && $alignment_params->{tool}) {
        return lc($alignment_params->{tool});
    }

    return 'gsnap'; # default aligner if not specified
}

__PACKAGE__->meta->make_immutable;

1;

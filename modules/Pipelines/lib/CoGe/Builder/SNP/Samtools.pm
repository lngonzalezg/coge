package CoGe::Builder::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;
use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename);

use CoGe::Accessory::Jex;
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths);
use CoGe::Builder::CommonTasks;

our $CONF = CoGe::Accessory::Web::get_defaults();
our $JEX = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION   = 0.1;
    @ISA       = qw(Exporter);
    @EXPORT    = qw(run build);
}

sub build {
    my $opts = shift;

    # Required arguments
    my $genome = $opts->{genome};
    my $input_file = $opts->{input_file}; # path to bam file
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $metadata = $opts->{metadata};
    my $params = $opts->{params};
    my $skipAnnotations = $opts->{skipAnnotations};

    # Setup paths
    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.faa";
    
    my $annotations = generate_additional_metadata($params);

    my $conf = {
        staging_dir => $staging_dir,
        result_dir  => $result_dir,

        bam         => $input_file,
        fasta       => catfile($FASTA_CACHE_DIR, $reheader_fasta),
        bcf         => catfile($staging_dir, qq[snps.raw.bcf]),
        vcf         => catfile($staging_dir, qq[snps.flt.vcf]),

        username    => $user->name,
        metadata    => $metadata,
        wid         => $wid,
        gid         => $gid,
        
        params      => $params,
        
        annotations => ($skipAnnotations ? '' : join(';', @$annotations))
    };

    # Build the workflow's tasks
    my @tasks;
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta_file,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );

    push @tasks, create_find_snps_job($conf);
    
    push @tasks, create_filter_snps_job($conf);
    
    my $load_vcf_task = create_load_vcf_job($conf);
    push @tasks, $load_vcf_task;

    return {
        tasks => \@tasks,
        metadata => $annotations,
        done_files => [ $load_vcf_task->{outputs}->[1] ]
    };
}

sub create_find_snps_job {
    my $opts = shift;

    # Required arguments
    my $reference = $opts->{fasta};
    my $alignment = $opts->{bam};
    my $snps = $opts->{bcf};

    # Pipe commands together
    my $sam_command = $CONF->{SAMTOOLS} || "samtools";
    $sam_command .= " mpileup -u -f " . basename($reference) . ' ' . basename($alignment);
    my $bcf_command = $CONF->{BCFTOOLS} || "bcftools";
    $bcf_command .= " view -b -v -c -g";

    # Get the output filename
    my $output = basename($snps);

    return {
        cmd => qq[$sam_command | $bcf_command - > $output],
        inputs => [
            $reference,
            $reference . '.fai',
            $alignment
        ],
        outputs => [ $snps ],
        description => "Identifying SNPs using SAMtools method ..."
    };
}

sub create_filter_snps_job {
    my $opts = shift;

    # Required arguments
    my $snps = $opts->{bcf};
    my $filtered_snps = $opts->{vcf};

    # Optional arguments
    my $params = $opts->{params};
    my $min_read_depth = $params->{'min-read-depth'} || 6;
    my $max_read_depth = $params->{'max-read-depth'} || 10;

    # Pipe commands together
    my $bcf_command = $CONF->{BCFTOOLS} || "bcftools";
    $bcf_command .= " view " . basename($snps);
    my $vcf_command = $CONF->{VCFTOOLS} || "vcfutils.pl";
    $vcf_command .= " varFilter -d $min_read_depth -D $max_read_depth";

    # Get the output filename
    my $output = basename($filtered_snps);

    return {
        cmd => qq[$bcf_command | $vcf_command > $output],
        inputs  => [ $snps ],
        outputs => [ $filtered_snps ],
        description => "Filtering SNPs ..."
    };
}

sub generate_additional_metadata {
    my $params = shift;
    my $min_read_depth = $params->{'min-read-depth'} || 6;
    my $max_read_depth = $params->{'max-read-depth'} || 10;
    return [ qq{note|SNPs generated using SAMtools method, min read depth $min_read_depth, max read depth $max_read_depth} ];
}

1;

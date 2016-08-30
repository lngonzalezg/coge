#!/usr/bin/perl -w

use strict;

use Data::Dumper;
use CoGeX;
use POSIX;

use Getopt::Long;

my ($dsgid, $delete_seqs, $db, $user, $pass, $space_min, $space_max, $length_min, $length_max, $use_strand);

GetOptions (
    "dsgid=i"=>\$dsgid,
    "database|db=s"=>\$db,
    "user|u=s"=>\$user,
    "password|pw=s"=>\$pass,
    "space_min|smin=i"=>\$space_min,
    "space_max|smax=i"=>\$space_max,
    "length_min|lmin=i"=>\$length_min,
    "length_max|lmax=i"=>\$length_max,
    "strand|s=i"=> \$use_strand,
);

unless ($dsgid && $db && $user && $pass) {
    print qq{
welcome to $0

Usage:  $0 -dsgid <database id for dataset group> -db <database name> -u <database user name> -pw <database password>

This program will generate fake quant data for a genome in CoGe.  Output is to STDOUT as a tab delimited file

 Options:

  -space_min | smin       min space (nt) between quantitative measurements (default 1)

  -space_max | smax       max space (nt) between quantitative measurements (default 2)

  -length_min | lmin       min length (nt) of quantitative measurements (default 0)  Note:  min of length 0 means the measurement will span one nt.  Nucleotide math is hard.

  -length_max | lmax       max length (nt) of quantitative measurements (default 2)

  -strand                  use this strand instead of random
};
    exit;
}

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3306";
my $coge = CoGeX->connect($connstr, $user, $pass );

my $dsg = $coge->resultset('Genome')->find($dsgid);
unless ($dsg) {
    print "unable to find entry for $dsgid\n";
    exit;
}
print join ("\t",qw(#CHR START STOP STRAND VALUE)), "\n";
gen_fake_data(dsg=>$dsg, space_min=>$space_min, space_max=>$space_max, length_min=>$length_min, length_max=>$length_max);

sub gen_fake_data {
    my %opts = @_;
    my $dsg = $opts{dsg};
    my $spacing_min = $opts{space_min}; #min space (nt) between quants
    my $spacing_max = $opts{space_max}; #max space (nt) between quants
    my $length_min = $opts{length_min}; #min length (nt) of the quant
    my $length_max = $opts{length_max}; #max length (nt) of the quant

    $spacing_min = 1 unless $spacing_min;
    $spacing_min = 1 if $spacing_min < 1; #can't be less than 1nt away
    $spacing_max = 2 unless $spacing_max;
    $spacing_max = $spacing_min if $spacing_max < $spacing_min; #value validation

    $length_min = 0 if not defined $length_min; #unless $length_min;
    $length_min = 0 if $length_min < 0; #can't be less than 0nt in size
    $length_max = 2 if not defined $length_max; #unless $length_max;
    $length_max = $length_min if $length_max < $length_min; #value validation

    print STDERR "smin=$spacing_min smax=$spacing_max lmin=$length_min lmax=$length_max\n";

    foreach my $gs ($dsg->genomic_sequences) {
    	my $pos = 1; #starting position
    	my $last = $gs->sequence_length; #get last position on chromosome;
    	my $chr = $gs->chromosome;
    	while ($pos <=$last) {
    	    my $spacing;
    	    $spacing= $spacing_min if $spacing_min == $spacing_max;
    	    $spacing = int(rand($spacing_max-$spacing_min+1))+$spacing_min unless $spacing;
    	    my $length;
    	    $length = $length_min if $length_min == $length_max;
    	    $length = int(rand($length_max-$length_min+1))+$length_min unless defined $length;
    	    my $val1 = rand();
    	    my $val2 = rand(10);
    	    my $strand = int(rand(2));
                $strand = $use_strand if defined $use_strand;
    	    $strand = -1 unless $strand; #if == 0, set to minus strand
    	    my $start = $pos+$spacing;
                last if $start > $last;
    	    my $stop = $start+$length;
                $stop = $last if ($stop > $last);
    	    $pos += $spacing+$length;
    
    	    print join (",", $chr, $start, $stop, $strand, $val1, $val2),"\n";
    	}
    }
}

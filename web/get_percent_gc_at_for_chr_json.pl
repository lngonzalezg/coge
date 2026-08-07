#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use CoGe::Accessory::Validate qw(valid_id valid_filename contained_path);
use Data::Dumper;
use File::Spec::Functions;
use CGI;

my $q = CGI->new;

# §6.3 Validate inputs and enforce genome access (see get_seq_for_chr.pl).
my $gid   = valid_id( $q->param('gid') );
my $chr   = valid_filename( $q->param('chr') );
my $wsize = valid_id( $q->param('wsize') );
my $wstep = valid_id( $q->param('wstep') );
unless ( defined $gid && defined $chr && defined $wsize && defined $wstep ) {
    print "Content-Type: application/json\n\n{}";
    exit;
}

my ( $db, $user, $conf ) = CoGe::Accessory::Web->init( cgi => $q );
my $genome = $db->resultset('Genome')->find($gid);
if (
    !$genome
    || ( $genome->restricted
        && ( !$user || $user->is_public || !$user->has_access_to_genome($genome) ) )
  )
{
    print "Content-Type: application/json\n\n{}";
    exit;
}

my $filename = $gid . "_" . $chr . "_" . $wsize . "_" . $wstep . "_out.txt";

print "Content-Type: application/json\n\n";
my $path = catfile( $conf->{SECTEMPDIR}, "downloads/genome", $gid );
my $file = contained_path( $path, $filename );
unless ( defined $file && -r $file ) {
    print "{}";
    exit;
}
open( my $fh, '<', $file ) or do { print "{}"; exit; };
my $l = <$fh>;    # skip header line
my ( @at, @gc, @n, @x );
while ( $l = <$fh> ) {
    chomp $l;
    my @tokens = split /\t/, $l;
    push @at, $tokens[2];
    push @gc, $tokens[3];
    push @n,  $tokens[4];
    push @x,  $tokens[5];
}
print '{"at":[';
print join( ',', @at );
print '],"gc":[';
print join( ',', @gc );
print '],"n":[';
print join( ',', @n );
print '],"x":[';
print join( ',', @x );
print ']}';
close $fh;

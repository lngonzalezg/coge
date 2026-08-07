#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use CoGe::Accessory::Validate qw(valid_id valid_filename contained_path);
use Data::Dumper;
use File::Spec::Functions;
use CGI;

my $q = CGI->new;

# §6.3 Validate inputs and enforce genome access (see get_seq_for_chr.pl).
my $gid = valid_id( $q->param('gid') );
my $chr = valid_filename( $q->param('chr') );
my $ws  = valid_id( $q->param('ws') );
unless ( defined $gid && defined $chr && defined $ws ) {
    print "Content-Type: text/plain\n\nInvalid request\n";
    exit;
}

my ( $db, $user, $conf ) = CoGe::Accessory::Web->init( cgi => $q );
my $genome = $db->resultset('Genome')->find($gid);
unless ($genome) {
    print "Content-Type: text/plain\n\nNot found\n";
    exit;
}
if ( $genome->restricted
    && ( !$user || $user->is_public || !$user->has_access_to_genome($genome) ) )
{
    print "Content-Type: text/plain\n\nAccess denied\n";
    exit;
}

my $filename = $gid . "_" . $chr . "_" . $ws . "_out.txt";

my $path = catfile( $conf->{SECTEMPDIR}, "downloads/genome", $gid );
my $file = contained_path( $path, $filename );
unless ( defined $file && -r $file ) {
    print "Content-Type: text/plain\n\nNot found\n";
    exit;
}

print "Content-Type: application/force-download\n";
print "Content-disposition: attachement; filename=chromosome_";
print $filename;
print "\n\n";

open( my $fh, '<', $file ) or exit;
while ( my $l = <$fh> ) {
    print $l;
}
close $fh;

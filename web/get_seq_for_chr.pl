#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use CoGe::Accessory::Validate qw(valid_id valid_filename contained_path);
use Data::Dumper;
use File::Spec::Functions;
use CGI;

my $q = CGI->new;

# §6.3 gid/chr were interpolated into a filesystem path with no validation,
# auth, or access check: chr could traverse (suffix-limited) and any restricted
# genome's chromosome sequence was downloadable anonymously (cross-genome IDOR).
my $gid = valid_id( $q->param('gid') );
my $chr = valid_filename( $q->param('chr') );
unless ( defined $gid && defined $chr ) {
    print "Content-Type: text/plain\n\nInvalid request\n";
    exit;
}

# Resolve db/user/conf from the session and enforce access to restricted genomes.
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

my $filename = $gid . "_" . $chr . ".faa";

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

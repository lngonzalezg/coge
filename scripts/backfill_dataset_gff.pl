#!/usr/bin/perl -w
#
# Backfill the preserved-annotation triplet for a dataset loaded BEFORE file
# preservation existed (pre 2026-08-18), so JBrowse2 can render it:
#
#   DATADIR/annotation/<tiered dsid>/<dataset.name>              GFF (DB export)
#   DATADIR/annotation/<tiered dsid>/<dataset.name>.sorted.gz    bgzip, sorted
#   DATADIR/annotation/<tiered dsid>/<dataset.name>.sorted.gz.csi
#   DATADIR/annotation/<tiered dsid>/.provenance                 how it got here
#
# PROVENANCE MATTERS: unlike files written by load_annotation.pl, the GFF here
# is a fresh export from the database (Dataset->gff), NOT the file the user
# originally uploaded -- byte-different, semantically equivalent. The
# .provenance file records that.
#
# Idempotent and race-safe: exits successfully if the triplet already exists;
# all writes go through a tmp name in the same directory and rename into
# place, so concurrent runs converge on identical files.
#
# Run by the JEX task submitted from GenomeView2's "Visualize older
# annotations" button (job type index_annotation); also runnable by hand:
#   perl backfill_dataset_gff.pl -dsid 12345 -config /opt/apache2/coge/coge.conf

use strict;
use Getopt::Long;
use File::Path qw(mkpath);
use File::Spec::Functions qw(catfile);
use String::ShellQuote qw(shell_quote);

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw(get_dataset_source_path);

our ($dsid, $config, $overwrite);
GetOptions(
    "dsid=i"      => \$dsid,
    "config=s"    => \$config,
    "overwrite=i" => \$overwrite,
);

unless ($dsid) {
    print STDOUT "log: error: dataset not specified, use -dsid\n";
    exit(-1);
}

my $P = CoGe::Accessory::Web::get_defaults($config);
my $connstr = "dbi:mysql:dbname=$P->{DBNAME};host=$P->{DBHOST};port=$P->{DBPORT};";
my $coge = CoGeX->connect( $connstr, $P->{DBUSER}, $P->{DBPASS} );
unless ($coge) {
    print STDOUT "log: error: couldn't connect to database\n";
    exit(-1);
}

my $ds = $coge->resultset('Dataset')->find($dsid);
if ( !$ds or $ds->deleted ) {
    print STDOUT "log: error: dataset $dsid not found\n";
    exit(-1);
}

my $dir = get_dataset_source_path($dsid);
unless ($dir) {
    print STDOUT "log: error: cannot determine storage path (DATADIR unset?)\n";
    exit(-1);
}

my $gff_file    = catfile( $dir, $ds->name );
my $sorted_file = $gff_file . '.sorted.gz';
my $csi_file    = $sorted_file . '.csi';

if ( -s $sorted_file and -s $csi_file and !$overwrite ) {
    print STDOUT "log: Annotation files already exist, nothing to do\n";
    exit;
}

print STDOUT 'log: Exporting dataset "' . $ds->name . "\" (id $dsid) from the database\n";
mkpath($dir) unless -d $dir;

# Export params mirror scripts/coge_gff.pl's defaults for a plain full export
# (all features, annotations included, IDs uniquified numerically).
my $tmp_gff = "$gff_file.tmp.$$";
open( my $fh, '>', $tmp_gff )
    or do { print STDOUT "log: error: cannot write '$tmp_gff': $!\n"; exit(-1); };
print $fh $ds->gff(
    print       => 0,
    annos       => 1,
    cds         => 0,
    name_unique => 0,
    id_type     => 0,
    base_url    => $P->{SERVER},
);
close($fh);
unless ( -s $tmp_gff ) {
    unlink $tmp_gff;
    print STDOUT "log: error: export produced an empty GFF (dataset has no features?)\n";
    exit(-1);
}

# Same recipe as load_annotation.pl's preserve block: coordinate-sorted bgzip
# + CSI index (CSI, not TBI: TBI cannot index sequences >2^29-1 bp and CoGe
# hosts plant chromosomes bigger than that; JBrowse2 reads CSI).
print STDOUT "log: Indexing (bgzip + CSI) for JBrowse2\n";
my $tmp_sorted = "$sorted_file.tmp.$$";
my $q_src = shell_quote($tmp_gff);
my $q_out = shell_quote($tmp_sorted);
my $q_dir = shell_quote($dir);
my $rc = system( '/bin/bash', '-c',
      "set -o pipefail; "
    # The awk filter drops malformed fragments: legacy Dataset->gff can emit
    # annotation values containing raw newlines, splitting a record into a
    # valid-prefix line plus a garbage fragment with <9 tab fields. tabix
    # otherwise warns ("Failed to parse TBX_GENERIC") and mis-parses.
    . "(grep '^#' $q_src; grep -v '^#' $q_src | awk -F'\\t' 'NF>=9' | LC_ALL=C sort -T $q_dir -t \$'\\t' -k1,1 -k4,4n) "
    . "| bgzip -c > $q_out && tabix -C -p gff $q_out" );
unless ( $rc == 0 and -s $tmp_sorted and -s "$tmp_sorted.csi" ) {
    unlink $tmp_gff, $tmp_sorted, "$tmp_sorted.csi";
    print STDOUT "log: error: bgzip/index failed (rc=$rc)\n";
    exit(-1);
}

# Atomic placement, GFF last: the listing endpoint keys renderability on the
# tabix pair, so a reader can never see a partial triplet that it would
# consider complete.
rename( "$tmp_sorted.csi", $csi_file )   or do { print STDOUT "log: error: rename csi: $!\n";    exit(-1); };
rename( $tmp_sorted,       $sorted_file ) or do { print STDOUT "log: error: rename sorted: $!\n"; exit(-1); };
rename( $tmp_gff,          $gff_file )    or do { print STDOUT "log: error: rename gff: $!\n";    exit(-1); };

open( my $pfh, '>', catfile( $dir, '.provenance' ) );
if ($pfh) {
    print $pfh "Backfilled from the CoGe database by backfill_dataset_gff.pl on "
        . localtime() . ".\n"
        . "The GFF is a database export (Dataset->gff), NOT the originally uploaded file.\n";
    close($pfh);
}

print STDOUT "log: Backfill complete\n";
exit;

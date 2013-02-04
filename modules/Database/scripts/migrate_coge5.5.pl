#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./migrate_coge5.5.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
# NOTE: CoGeX code must be updated to v5 BEFORE running this script!
#
# Migrate sharing-related database tables from v5 to v5.5 (makes things work
# with new user_connector table).
#
#-------------------------------------------------------------------------------

use DBI;
use strict;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge);

my ($db, $user, $pass, $go);
GetOptions (
	"debug=s" => \$DEBUG,
	"database|db=s"=>\$db,
	"user|u=s"=>\$user,
	"password|pw=s"=>\$pass,
	"go=i"=>\$go
);

$| = 1;
$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on

print STDERR "Running $0\n";

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dbh = $coge->storage->dbh;

#-------------------------------------------------------------------------------
# Create new tables
#-------------------------------------------------------------------------------

# Create user_connector table
#sql('drop table if exists user_connector');
#sql(<<SQL);
#create table user_connector
#(
#	user_connector_id INT(11) NOT NULL AUTO_INCREMENT,
#	parent_id INT(11) NOT NULL,
#	parent_type TINYINT(1) NOT NULL,
#	child_id INT(11) NOT NULL,
#	child_type TINYINT(1) NOT NULL,
#	role_id INT(11) NOT NULL DEFAULT 4,
#	PRIMARY KEY (user_connector_id)
#);
#SQL

#-------------------------------------------------------------------------------
# Migrate data:  remap group/list connectors to user_connector
#-------------------------------------------------------------------------------

foreach my $group ($coge->resultset('UserGroup')->all) {
	print STDERR "Group id" . $group->id . " " . $group->name . ":" . $group->description . "\n";
	if ($group->is_owner) {
		# Map owner list contents to users in user_connector
		my $owner_list = $group->owner_list;
		my $children = $owner_list->children_by_type;

		if (not keys %$children) {
			print STDERR "   empty, do nothing\n";
		}
		else {
			foreach my $type (keys %$children) {
				foreach my $child (@{$children->{$type}}) {
					print STDERR "   owner conn parent_id=" . $group->creator_user_id . " child_id=" . $child->id . "\n";
					if ($go) {
						my $conn = $coge->resultset('UserConnector')->find_or_create( { 
							parent_id => $group->creator_user_id,
							parent_type => 5, #FIXME hardcoded to "user"
							child_id => $child->id,
							child_type => $type,
							role_id => 2 #FIXME hardcoded to "owner"
						} );
						die unless $conn;
					}
				}
			}
		}

		# Delete owner group/list (auto-cascades down owner list and associated list_connectors)
		#$group->delete;
	}
	else {
		# Map non-owner lists to parent groups in user_connector
		foreach my $list ($group->lists) {
			print STDERR "   non-owner conn parent_id=" . $group->id . " child_id=" . $list->id . "\n";
			if ($go) {
				my $conn = $coge->resultset('UserConnector')->find_or_create( { 
					parent_id => $group->id,
					parent_type => 6, #FIXME hardcoded to "group"
					child_id => $list->id,
					child_type => 1, # FIXME hardcoded to "list"
					role_id => $group->role_id
				} );
				die unless $conn;
			}
		}
		
		# Map group member users (user_group_connector) in user_connector;
		foreach my $user ($group->users) {
			print STDERR "   group parent_id=" . $group->id . " child_id=" . $user->id . "\n";
			if ($go) {
				my $conn = $coge->resultset('UserConnector')->find_or_create( { 
					parent_id => $group->id,
					parent_type => 6, #FIXME hardcoded to "group"
					child_id => $user->id,
					child_type => 5, #FIXME hardcoded to "user"
					role_id => ($user->id == $group->creator_user_id ? 2 : $group->role_id) # FIXME hardcoded to "owner"
				} );
				die unless $conn;
			}
		}
	}
}

# Remove "locked" column from list table
#drop_column("list", "locked");
# Remove "user_group_id" column from list table
#drop_column("list", "user_group_id");
# Remove "locked" column from user_group
#drop_column("user_group", "locked");
# Remove "creator_user_id" from user_group table
#drop_column("user_group", "creator_user_id");
# Remove user_group_connector table
#sql("drop table user_group_connector");

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

$coge->resultset('Log')->create( { user_id => 0, page => $0, description => 'database migrated to coge5.5' } );
print STDERR "All done!\n";
exit;


#-------------------------------------------------------------------------------
sub sql {
	my $cmd = shift;
	print $cmd, "\n" if ($DEBUG);
	$dbh->do($cmd);
}

sub drop_column {
	my $table = shift;
	my $column = shift;

	my $sth = $dbh->prepare("SELECT * FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA='$db' AND TABLE_NAME='$table' AND COLUMN_NAME='$column'");
	my $count = $sth->execute();
	if ($count > 0) {
		sql("alter table $table drop $column");
	}
	$sth->finish();
}

sub add_column {
	my $table = shift;
	my $column = shift;

	my $sth = $dbh->prepare("SELECT * FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA='$db' AND TABLE_NAME='$table' AND COLUMN_NAME='$column'");
	my $count = $sth->execute();
	if ($count == 0) {
		sql("alter table $table add $column");
	}
	$sth->finish();
}


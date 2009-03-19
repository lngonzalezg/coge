package CoGeX::UserGroupConnector;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::UserGroupConnector

=head1 SYNOPSIS

  use CoGeX::UserGroupConnector;
This object uses the DBIx::Class to define an interface to the C<user_group_connector> table in the CoGe database.


=head1 DESCRIPTION


Has columns:
C<user_group_connector_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<user_id>
Type: INT, Default: "", Nullable: no, Size: 10

C<user_group_id>
Type: INT, Default: "", Nullable: no, Size: 10


=head1 USAGE

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("user_group_connector");
__PACKAGE__->add_columns(
  "user_group_connector_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
  "user_group_id",
  { data_type => "INT", default_value => "", is_nullable => 0, size => 10 },
);
__PACKAGE__->set_primary_key("user_group_connector_id");

1;



=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut

package CoGeX::Location;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::

=head1 SYNOPSIS

  use CoGeX::
This object uses the DBIx::Class to define an interface to the C<AnnotationsType> table in the CoGe database.


=head1 DESCRIPTION


Has columns:
C<location_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<start>
Type: INT, Default: 0, Nullable: no, Size: 11

C<stop>
Type: INT, Default: 0, Nullable: no, Size: 11

C<chromosome>
Type: VARCHAR, Default: "", Nullable: no, Size: 255

C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<strand>
Type: TINYINT, Default: "", Nullable: no, Size: 4

Belongs to C<CoGeX::Feature> via C<feature_id>


=head1 USAGE

=head1 METHODS

=cut


__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("location");
__PACKAGE__->add_columns(
  "location_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "start",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "stop",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "chromosome",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "strand",
  { data_type => "TINYINT", default_value => "", is_nullable => 0, size => 4 },
);
__PACKAGE__->set_primary_key("location_id");

__PACKAGE__->belongs_to("feature" => "CoGeX::Feature", 'feature_id');
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

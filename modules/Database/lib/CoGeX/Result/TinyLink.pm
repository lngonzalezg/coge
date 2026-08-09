package CoGeX::Result::TinyLink;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::TinyLink

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<tiny_link> table
in the CoGe database -- the internal link shortener that replaced the external
YOURLS service.

=head1 DESCRIPTION

Has columns:
C<keyword> (Primary Key)
Type: CHAR, Default: undef, Nullable: no, Size: 12

The short key itself, and also the hash of C<rel_url>:
C<lc(base36(md5(rel_url)))> truncated to 12 characters (~62 bits). It is derived,
never allocated, so the same URL always maps to the same key -- across processes,
across time, and across a wiped database. Several consumers depend on that:
C<CoGe::Accessory::Web::get_job> dedups prior jobs by link, SynMap names its result
files after the key (C<< <key>.log >>, C<< dotplot_dots_<key>.cfg >> under DIAGSDIR)
and re-finds them on an identical re-run, and SynFind names its workflow
C<< synfind-<key> >> after extracting the key with C</(\w+)$/>.

C<rel_url>
Type: TEXT, Default: undef, Nullable: no, Size: N/A

The target, stored RELATIVE to the C<SERVER> config value (e.g.
C<"SynMap.pl?dsgid1=..;dsgid2=..">). The redirect re-attaches the current C<SERVER>,
so links survive a domain rename or an http->https switch, scheme/host aliases of the
same page converge on one key, and an off-site target cannot be minted at all.

C<created_at>
Type: TIMESTAMP, Default: CURRENT_TIMESTAMP, Nullable: no

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("tiny_link");
__PACKAGE__->add_columns(
  "keyword",  { data_type => "CHAR", default_value => undef, is_nullable => 0, size => 12 },
  "rel_url",  { data_type => "TEXT", default_value => undef, is_nullable => 0 },
  "created_at",  { data_type => "TIMESTAMP", default_value => \"CURRENT_TIMESTAMP", is_nullable => 0 },
);
__PACKAGE__->set_primary_key("keyword");

1;

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

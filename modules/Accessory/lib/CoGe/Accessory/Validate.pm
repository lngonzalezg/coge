package CoGe::Accessory::Validate;

=head1 NAME

CoGe::Accessory::Validate

=head1 SYNOPSIS

Narrow, fail-closed input validators for security-sensitive boundaries.

    use CoGe::Accessory::Validate qw(valid_id valid_ident valid_filename
                                     valid_relpath valid_enum contained_path);

    my $gid = valid_id($self->param('gid'))
        or return $self->render(status => 400, text => 'bad id');

=head1 DESCRIPTION

Added in security pass 1. These replace the single permissive C<check_taint>
regex (see CoGe::Accessory::Web) with per-context checks that return the
validated value on success and C<undef> on failure, so a caller that does not
test the result fails closed rather than open.

None of these mutate global state and none die; the caller decides how to reject.

=head1 AUTHOR

Security remediation, 2026.

=cut

use strict;
use warnings;

use File::Basename qw(basename);
use File::Spec;
use Cwd qw(realpath);

BEGIN {
    use vars qw($VERSION @ISA @EXPORT_OK);
    require Exporter;
    $VERSION   = 0.1;
    @ISA       = qw(Exporter);
    @EXPORT_OK = qw(
        valid_id valid_ident valid_filename valid_relpath valid_enum contained_path
    );
}

# Non-negative integer id. Returns the (untainted) digit string, or undef.
sub valid_id {
    my ($v) = @_;
    return unless defined $v;
    return unless $v =~ /^(\d+)$/;
    return $1;
}

# Bare identifier: letters, digits, underscore. Returns the value, or undef.
sub valid_ident {
    my ($v) = @_;
    return unless defined $v;
    return unless $v =~ /^(\w+)$/;
    return $1;
}

# Single path component (no directory separators). basename() strips any leading
# path, then the result must be a plain safe name and never '.'/'..'. Returns the
# cleaned basename, or undef.
sub valid_filename {
    my ($v) = @_;
    return unless defined $v && length $v;
    my $name = basename($v);
    return if $name eq '' || $name eq '.' || $name eq '..';
    return unless $name =~ /^([\w.\-]+)$/;
    my $clean = $1;
    return if $clean =~ /\.\./;    # defensive: no '..' anywhere
    return $clean;
}

# Relative multi-segment path where every segment is a valid_filename. Rejects
# absolute paths, '.'/'..' segments, empty segments, backslashes. Needed for
# legitimately slash-bearing request params (e.g. CoGeBlast logfile=tmp/CoGeBlast/x.log).
# Returns the normalized 'a/b/c' string, or undef.
sub valid_relpath {
    my ($v) = @_;
    return unless defined $v && length $v;
    return if $v =~ m{\\};                 # no backslashes
    return if $v =~ m{^/};                 # no absolute paths
    my @segs = split m{/+}, $v;
    return unless @segs;
    my @clean;
    for my $s (@segs) {
        next if $s eq '';                  # tolerate doubled/trailing slashes
        my $c = valid_filename($s);
        return unless defined $c;
        push @clean, $c;
    }
    return unless @clean;
    return join('/', @clean);
}

# Membership test against an explicit allowlist. Returns the value if it matches
# one of @allowed exactly, else undef.
sub valid_enum {
    my ($v, @allowed) = @_;
    return unless defined $v;
    for my $a (@allowed) {
        return $v if $v eq $a;
    }
    return;
}

# Join @parts under $base and confirm the resolved absolute path stays inside the
# resolved $base. Returns the contained absolute path, or undef on escape. The path
# need not exist yet: we resolve the deepest existing ancestor and re-append the rest.
sub contained_path {
    my ($base, @parts) = @_;
    return unless defined $base && length $base;

    my $base_real = realpath($base);
    return unless defined $base_real;      # base must exist

    my $joined = File::Spec->catfile($base, @parts);

    # Resolve as far as the filesystem allows, then rebuild the tail so a
    # not-yet-created target still gets a canonical, symlink-free prefix.
    my @accum = File::Spec->splitdir(File::Spec->canonpath($joined));
    my @tail;
    my $resolved;
    while (@accum) {
        my $cand = File::Spec->catdir(@accum);
        $cand = File::Spec->rootdir() if $cand eq '';
        if (defined(my $r = realpath($cand))) {
            $resolved = $r;
            last;
        }
        unshift @tail, pop @accum;
    }
    return unless defined $resolved;
    my $final = @tail ? File::Spec->catfile($resolved, @tail) : $resolved;

    # Containment: final must equal base_real or sit beneath base_real . '/'.
    my $prefix = $base_real;
    $prefix .= '/' unless $prefix =~ m{/$};
    return unless $final eq $base_real || index($final, $prefix) == 0;

    return $final;
}

1;

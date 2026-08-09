package CoGe::Services::API::Link;

use Mojo::Base 'Mojolicious::Controller';

use CoGeX;
use CoGe::Accessory::Web qw(get_defaults get_tiny_link);
use CoGe::Services::Error;

###############################################################################
# Internal link shortener -- the CoGe-side replacement for YOURLS.
#
#   GET /r/:keyword        resolve() -- 302 to the stored page (proxied from
#                          Apache as /<base>/r/<keyword>, see coge-main.conf)
#   GET /api/v1/links?url= create()  -- mint a link from the browser
#
# Minting itself lives in CoGe::Accessory::Web::get_tiny_link, so the key this
# endpoint returns is byte-identical to the one the server-side Perl call sites
# produce for the same page. SynMap depends on that: the browser puts the link in
# the job payload as `tinylink`, and the pipeline names its result files after the
# key it would have computed itself.
###############################################################################

# The stored URL is relative, so a link can only ever point back into this site --
# no target validation is needed on resolve, and create() needs no authentication.
sub resolve {
    my $self    = shift;
    my $keyword = $self->stash('keyword');

    my $conf = get_defaults();
    my $db   = CoGeX->dbconnect($conf);

    my $tiny_link;
    $tiny_link = $db->resultset('TinyLink')->find($keyword)
      if $db && defined $keyword && length $keyword;

    unless ($tiny_link) {
        return $self->render(
            status => 404,
            format => 'html',
            text   => _not_found_page()
        );
    }

    my $server = _server_base($conf);
    $self->redirect_to( $server . $tiny_link->rel_url );
}

sub create {
    my $self = shift;
    my $url  = $self->param('url');

    unless ( defined $url && length $url ) {
        return $self->render( API_STATUS_BAD_REQUEST('Missing url parameter') );
    }

    # Server-side callers are trusted and always shorten their own pages. This
    # endpoint is public, so an absolute URL must prove it belongs to this site.
    # Relative URLs need no check: they are resolved against SERVER on redirect.
    if ( $url =~ m{^\w+://} ) {
        unless ( _is_local_url( $url, get_defaults() ) ) {
            return $self->render(
                API_STATUS_BAD_REQUEST('URL must point at this CoGe server') );
        }
    }

    my $link = get_tiny_link( url => $url );

    $self->render( json => { link => $link } );
}

sub _server_base {
    my $conf = shift;
    my $server = ( $conf && $conf->{SERVER} ) ? $conf->{SERVER} : '/';
    $server =~ s{/*$}{/};
    return $server;
}

sub _is_local_url {
    my ( $url, $conf ) = @_;
    return 0 unless $conf;

    my ($host_and_path) = $url =~ m{^\w+://(.*)$};
    return 0 unless defined $host_and_path;

    # Scheme is deliberately ignored: http and https forms of the same host are
    # the same site, and they normalize to the same relative URL anyway.
    foreach my $key (qw(SERVER INT_SERVER)) {
        my $base = $conf->{$key};
        next unless $base;
        $base =~ s{^\w+://}{};
        $base =~ s{/*$}{/};
        return 1 if index( lc $host_and_path, lc $base ) == 0;
    }

    return 0;
}

sub _not_found_page {
    return <<'HTML';
<!DOCTYPE html>
<html>
<head><title>CoGe - Link not found</title></head>
<body>
<h1>This short link is no longer valid</h1>
<p>CoGe's short links were reset when the old link service was retired, so links
created before that point no longer resolve.</p>
<p>Please re-run the analysis to get a fresh link.</p>
</body>
</html>
HTML
}

1;

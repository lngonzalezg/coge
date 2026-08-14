#!/usr/bin/perl -w
#
# OIDC redirect URI -- Keycloak sends the browser back here with ?code=&state=
# (plus session_state/iss, which are ignored). This URL is covered by the
# client's registered https://genomevolution.org/coge* wildcard.
#
# A dedicated callback, unlike CAS which returned tickets to every page: the
# auth params never land on content pages, so they can never leak into tiny
# links or logs, and no per-page handling is needed.
#
# All real work happens in CoGe::Accessory::Web::login_oidc -- state/HMAC
# check, code exchange (client secret, server-side only), ID token
# verification, then the SAME add_user/session/cookie machinery CAS uses.

use strict;
use CGI;
use CoGe::Accessory::Web;
use CoGeX;

my $cgi  = CGI->new;
my $conf = CoGe::Accessory::Web::get_defaults();

my $server = $conf->{SERVER} || '/coge/';
$server =~ s{/*$}{/};

# User cancelled at the Keycloak screen (error=access_denied) or the IdP
# reported a failure: go home logged-out, no drama.
if ( my $err = $cgi->param('error') ) {
    print STDERR 'oidc_callback: IdP returned error: ', $err, "\n";
    print $cgi->redirect( -uri => $server );
    exit;
}

my $db = CoGeX->dbconnect($conf);
unless ($db) {
    print STDERR "oidc_callback: cannot connect to database\n";
    print $cgi->redirect( -uri => $server );
    exit;
}

my $result = CoGe::Accessory::Web::login_oidc(
    cgi   => $cgi,
    coge  => $db,
    code  => scalar $cgi->param('code'),
    state => scalar $cgi->param('state'),
);

if ( $result->{user} ) {
    print $cgi->redirect(
        -uri    => $result->{redirect} || $server,
        -cookie => $result->{cookies},
    );
    exit;
}

# Failed login: log the reason server-side, show the user a terse page with a
# retry link. Deliberately no detail client-side -- the reasons (state
# mismatch, audience mismatch, ...) are only useful to an attacker or a log.
print STDERR 'oidc_callback: login failed: ', ( $result->{error} // 'unknown' ), "\n";
print $cgi->header( -status => 401 );
print <<"HTML";
<!DOCTYPE html>
<html><head><title>CoGe: sign-in failed</title></head>
<body style="font-family:sans-serif;margin:4em auto;max-width:36em;">
<h2>Sign-in didn't complete</h2>
<p>Something went wrong while signing you in with CyVerse. This is usually
transient &mdash; please try again.</p>
<p><a href="${server}oidc_login.pl">Try again</a> &middot;
<a href="$server">Continue without signing in</a></p>
</body></html>
HTML

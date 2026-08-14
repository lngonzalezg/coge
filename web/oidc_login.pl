#!/usr/bin/perl -w
#
# Start an OIDC login against CyVerse Keycloak.
#
# The header's "Log in" button points here (tmpl/footer.tmpl login_oidc()).
# oidc_begin() mints state/nonce/PKCE, binds them plus the relative return URL
# into an HMAC-signed cookie, and hands back the Keycloak authorize URL.
#
# FALLBACK: if the OIDC_* keys are absent from coge.conf, this redirects to the
# CAS login instead -- so the feature flag lives here, server-side, in exactly
# one place, rather than being threaded through every page template. Deploying
# this code with the keys unset changes nothing user-visible.

use strict;
use CGI;
use CoGe::Accessory::Web;

my $cgi = CGI->new;
my $return_to = $cgi->param('return_to') // '';

my ($cookie, $auth_url) = CoGe::Accessory::Web::oidc_begin( return_to => $return_to );

if ($cookie && $auth_url) {
    print $cgi->redirect( -uri => $auth_url, -cookie => [$cookie] );
    exit;
}

# OIDC not configured (or its secret unreadable): fall back to CAS, preserving
# the exact pre-OIDC behavior of the login button.
my $conf = CoGe::Accessory::Web::get_defaults();
my $server = $conf->{SERVER} || '';
my $service = $return_to || $server;
print $cgi->redirect( -uri => $conf->{CAS_URL} . '/login?service=' . $service );

#! /usr/bin/perl -w
#
# JBrowse2 genome viewer (Phase B pilot, 2026-08-18). Runs IN PARALLEL with
# the JBrowse1 GenomeView.pl -- nothing links here yet; reach it directly:
#
#   GenomeView2.pl?gid=<genome_id>[&embed=1]
#
# The page is a thin shell: all data access happens in the browser against the
# authenticated file gateway and listing endpoint (CoGe::Services::API::
# DataFiles), which enforce genome/dataset visibility per request. The page
# itself therefore renders for anyone; what appears in it depends on the
# session cookie the browser sends with each fetch.

use strict;

use CGI;
use HTML::Template;

use CoGeX;
use CoGe::Accessory::Web;

use vars qw( $CONF $PAGE_TITLE $USER $DB %FUNCTION $FORM $EMBED $LINK );

$PAGE_TITLE = 'GenomeView2';

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi        => $FORM,
    page_title => $PAGE_TITLE
);

%FUNCTION = ();
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed');
    if ($EMBED) {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            USER          => $USER->display_name || '',
            PAGE_TITLE    => 'Genome Viewer (JBrowse2)',
            PAGE_LINK     => $LINK,
            SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
            TITLE         => 'Genome Viewer (JBrowse2)',
            HOME          => './',
            HELP          => 'GenomeView',
            WIKI_URL      => $CONF->{WIKI_URL} || '',
            ADMIN_ONLY    => $USER->is_admin,
            CAS_URL       => $CONF->{CAS_URL} || '',
            NO_DOCTYPE    => 1,
            COOKIE_NAME   => $CONF->{COOKIE_NAME} || ''
        );
        $template->param( LOGON => 1 ) unless ( $USER->user_name eq 'public' );
    }

    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    return $template->output;
}

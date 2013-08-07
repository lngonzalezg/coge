#! /usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON::XS;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use File::Path;
use Sort::Versions;
no warnings 'redefine';

use vars qw( $P $PAGE_TITLE $USER $coge $FORM %FUNCTION );

$PAGE_TITLE = 'Experiments';

$FORM = new CGI;

( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    ticket     => $FORM->param('ticket'),
    url        => $FORM->url,
    page_title => $PAGE_TITLE
);

%FUNCTION = (
    delete_experiment        => \&delete_experiment,
    get_experiments_for_user => \&get_experiments_for_user,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub delete_experiment {
    my %opts = @_;
    my $eid  = $opts{eid};

    #print STDERR "delete_experiment: eid=$eid\n";
    return "Must have valid experiment id\n" unless ($eid);

    # Check permissions
    return unless ( $USER->is_admin or $USER->is_owner( experiment => $eid ) );

    # Delete the experiment and associated connectors & annotations
    #FIXME add some error checking/logging here
    my $experiment = $coge->resultset('Experiment')->find($eid);
    return 0 unless $experiment;
    $experiment->deleted(1);
    $experiment->update;

    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "$PAGE_TITLE",
        description => 'delete experiment id' . $experiment->id
    );

    return 1;
}

sub get_experiments_for_user {

    #my %opts = @_;

    my %experiments;
    if ( $USER->is_admin ) {
        map { $experiments{ $_->id } = $_ }
          $coge->resultset('Experiment')->all();
    }
    elsif ( $USER->id == 0 ) {    #FIXME: call $USER->is_public
        map { $experiments{ $_->id } = $_ }
          $coge->resultset('Experiment')->search( { restricted => 0 } );
    }
    else {
        map { $experiments{ $_->id } = $_ } $USER->experiments,
          $coge->resultset('Experiment')->search( { restricted => 0 } );
    }

    my @rows;
    foreach my $e ( sort experimentcmp values %experiments ) {
        my $user_can_edit = $USER->is_admin
          || $USER->is_owner_editor( experiment => $e->id );
        push @rows,
          {
            NAME =>
qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=}
              . $e->id
              . qq{")'>}
              . $e->info
              . "</span>",
            VERSION     => $e->version,
            DATE        => $e->date,
            EDIT_BUTTON => $user_can_edit
            ? "<span class='link ui-icon ui-icon-gear' onclick=\"window.open('ExperimentView.pl?eid="
              . $e->id
              . "')\"></span>"
            : "<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be edited.')\"></span>",
            DELETE_BUTTON => $user_can_edit
            ? "<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_experiment({eid: '"
              . $e->id
              . "'});\"></span>"
            : "<span class='link ui-icon ui-icon-locked' onclick=\"alert('This list is locked and cannot be deleted.')\"></span>"
          };
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( DO_EXPERIMENT_TABLE => 1 );
    $template->param( EXPERIMENT_LOOP     => \@rows );

    return $template->output;
}

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => $PAGE_TITLE );
    $template->param( HELP => '/wiki/index.php?title=' . $PAGE_TITLE . '.pl' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= ' ' . $USER->last_name
      if ( $USER->first_name && $USER->last_name );
    $template->param( USER     => $name );
    $template->param( LOGO_PNG => $PAGE_TITLE . "-logo.png" );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
    $link = CoGe::Accessory::Web::get_tiny_link( url => $link );

    $template->param( BODY       => gen_body() );
    $template->param( ADJUST_BOX => 1 );

    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( MAIN             => 1 );
    $template->param( PAGE_NAME        => $PAGE_TITLE );
    $template->param( EXPERIMENT_TABLE => get_experiments_for_user() );

    return $template->output;
}

# FIXME this comparison routine is duplicated elsewhere
sub experimentcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

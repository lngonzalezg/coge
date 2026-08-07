package CoGe::Services::API::JBrowse::Genome;

use Mojo::Base 'Mojolicious::Controller';
use CoGeX;
use CoGe::Services::Auth qw(init);
use Data::Dumper;

sub _add_features {
    my ($name, $chr, $type_ids, $dsid, $hits, $dbh) = @_;
    # Security pass 1 (S2): every value is bound via a placeholder. $type_ids is an
    # arrayref of integer feature_type_ids (the caller validates them); $chr and $name
    # were previously concatenated into single-quoted literals (unauthenticated SQLi).
    my @bind = ($dsid);
    my $query = 'SELECT name,chromosome,start,stop FROM feature JOIN feature_name on feature.feature_id=feature_name.feature_id WHERE dataset_id=?';
    if ($chr) {
        $query .= " AND chromosome=?";
        push @bind, $chr;
    }
    if ($type_ids && @$type_ids) {
        my $ph = join(',', ('?') x scalar @$type_ids);
        $query .= " AND feature_type_id IN($ph)";
        push @bind, @$type_ids;
    }
    my $op = (index($name, '%') != -1) ? 'LIKE' : '=';
    $query .= " AND lower(name) $op ?";
    push @bind, $name;
    my $sth = $dbh->prepare($query);
    $sth->execute(@bind);
    while (my $row = $sth->fetch) {
        push @$hits, { name => $row->[0], location => { ref => $row->[1], start => $row->[2], end => $row->[3] } };
    }
}

sub features {
    my $self = shift;
    my $name = scalar $self->param('name');
    my $chr = $self->param('chr');

    my ($db, $user) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

	my $types = $self->param('features');
	my $type_ids;   # arrayref of integer feature_type_ids
	if (defined $types && $types ne 'all') {
		# Security pass 1 (S2): the 'features' param (a comma list of type names) was
		# injected straight into IN(...). Parse the names, strip any quoting, and bind
		# each via a placeholder.
		my @names = grep { length } map { my $t = $_; $t =~ s/^\s*['"]?//; $t =~ s/['"]?\s*$//; $t } split(/,/, $types);
		if (@names) {
			my $ph = join(',', ('?') x scalar @names);
			$type_ids = $dbh->selectcol_arrayref('SELECT feature_type_id FROM feature_type WHERE name IN(' . $ph . ')', undef, @names);
		}
		$type_ids ||= [];
	}

    my $hits = [];
    # Security pass 1 (S6): bind gid (route is \d+, but never interpolate).
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=?', undef, $self->stash('gid'));
    foreach my $dsid (@$ids) {
        _add_features('%' . lc($name) . '%', $chr, $type_ids, $dsid, $hits, $dbh);
    }
    my @sorted = sort { $a->{name} cmp $b->{name} } @{$hits};
    $self->render(json => \@sorted);
}

sub genes {
    my $self = shift;
    my $name = scalar $self->param('equals');
    if (!$name) {
        $name = scalar $self->param('startswith');
        if ($name) {
            $name .= '%';
        }
    }
    
    if (index($name, ':') != -1 && index($name, '..') != -1) {
        $self->render(json => []);
        return;
    }

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

    my $hits = [];
    # Security pass 1 (S6): bind gid.
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=?', undef, $self->stash('gid'));
    foreach my $dsid (@$ids) {
        _add_features(lc($name), undef, [1], $dsid, $hits, $dbh);
    }
    $self->render(json => $hits);
}

1;

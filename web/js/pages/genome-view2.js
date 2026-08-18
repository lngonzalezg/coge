/*
 * JBrowse2 pilot page (Phase B, 2026-08-18).
 *
 * Everything here consumes the authenticated file gateway + listing endpoint
 * (CoGe::Services::API::DataFiles). All URLs are RELATIVE, hence same-origin,
 * hence the browser attaches the cogec session cookie to every fetch JBrowse2
 * makes -- auth and per-object visibility are enforced server-side on each
 * request, and this page contains no auth logic at all. When the gateway's
 * FILE_STORE driver moves to S3 (302 -> presigned URL), fetch follows the
 * redirects and this page ships unchanged.
 *
 * Track configs use Gff3TabixAdapter with a CSI index: TBI cannot index
 * sequences >2^29-1 bp and CoGe hosts plant chromosomes bigger than that.
 */
(function () {
    'use strict';

    function fail(msg) {
        var e = document.getElementById('jbrowse2_error');
        e.textContent = msg;
        e.style.display = 'block';
        console.error('GenomeView2: ' + msg);
    }

    var params = new URLSearchParams(window.location.search);
    var gid = params.get('gid');
    if (!gid || !/^\d+$/.test(gid)) return fail('Missing or invalid gid parameter.');

    var API = 'api/v1/';

    function fetchJSON(url) {
        return fetch(url).then(function (r) {
            if (!r.ok) throw new Error(url + ' -> HTTP ' + r.status);
            return r.json();
        });
    }

    Promise.all([
        fetchJSON(API + 'genomes/' + gid + '/datasets'),
        // The .fai is tiny text; its first line gives a sensible initial
        // location without any new server code.
        fetch(API + 'genomes/' + gid + '/files/fai').then(function (r) {
            if (!r.ok) throw new Error('fai -> HTTP ' + r.status +
                (r.status === 401 ? ' (restricted genome -- are you logged in?)' : ''));
            return r.text();
        })
    ]).then(function (results) {
        var listing = results[0];
        var fai = results[1];

        if (!listing.genome.files.fasta || !listing.genome.files.fai)
            return fail('This genome has no indexed FASTA in storage.');

        var first = fai.split('\n')[0].split('\t');   // name, length, ...
        var refName = first[0];
        var refLen = parseInt(first[1], 10) || 50000;
        var location = refName + ':1..' + Math.min(refLen, 50000);

        var assemblyName = 'coge_gid_' + gid;
        var assembly = {
            name: assemblyName,
            displayName: listing.genome.name,
            sequence: {
                type: 'ReferenceSequenceTrack',
                trackId: assemblyName + '-seq',
                adapter: {
                    type: 'IndexedFastaAdapter',
                    fastaLocation: { uri: API + 'genomes/' + gid + '/files/fasta', locationType: 'UriLocation' },
                    faiLocation:   { uri: API + 'genomes/' + gid + '/files/fai',   locationType: 'UriLocation' }
                }
            },
            // CoGe FASTAs carry NCBI-style prefixed names (lcl|LL0249_Chr01)
            // while loaded GFFs use bare names -- disjoint refName sets make
            // annotation tracks silently render nothing. The gateway derives
            // this aliases file from the .fai on the fly.
            refNameAliases: {
                adapter: {
                    type: 'RefNameAliasAdapter',
                    location: { uri: API + 'genomes/' + gid + '/files/aliases', locationType: 'UriLocation' }
                }
            }
        };

        // One annotation track per dataset whose preserved tabix pair exists.
        // Datasets without files (loaded before preservation, 2026-08-18) are
        // listed but not renderable yet -- the lazy backfill (Phase D) slots
        // in here later.
        var tracks = [];
        var skipped = [];
        listing.datasets.forEach(function (ds) {
            if (ds.files['gff-tabix'] && ds.files['gff-csi']) {
                tracks.push({
                    type: 'FeatureTrack',
                    trackId: 'dataset-' + ds.id,
                    name: ds.name + (ds.restricted ? ' (restricted)' : ''),
                    assemblyNames: [assemblyName],
                    adapter: {
                        type: 'Gff3TabixAdapter',
                        gffGzLocation: { uri: API + 'datasets/' + ds.id + '/files/gff-tabix', locationType: 'UriLocation' },
                        index: {
                            location: { uri: API + 'datasets/' + ds.id + '/files/gff-csi', locationType: 'UriLocation' },
                            indexType: 'CSI'
                        }
                    }
                });
            }
            else skipped.push(ds.name);
        });
        if (skipped.length)
            console.log('GenomeView2: no preserved files yet (pre-2026 loads): ' + skipped.join(', '));

        var JB = window.JBrowseReactLinearGenomeView;
        var state = JB.createViewState({
            assembly: assembly,
            tracks: tracks,
            location: location,
            // Match CoGe's design system (coge-modern.css green-8/green-9)
            // instead of JBrowse2's default MUI blue.
            configuration: {
                theme: {
                    palette: {
                        primary:   { main: '#2f9e44' },
                        secondary: { main: '#2b8a3e' }
                    }
                }
            }
        });
        // Show the reference sequence and every renderable annotation track by
        // default -- an empty view with a hidden track selector is a bad first
        // impression. Individually guarded: one bad track config should not
        // take down the rest of the view.
        [assemblyName + '-seq'].concat(tracks.map(function (t) { return t.trackId; }))
            .forEach(function (id) {
                try { state.session.view.showTrack(id); }
                catch (e) { console.error('GenomeView2: showTrack(' + id + '): ' + e.message); }
            });

        var root = ReactDOM.createRoot(document.getElementById('jbrowse2_view'));
        root.render(React.createElement(JB.JBrowseLinearGenomeView, { viewState: state }));
    }).catch(function (err) {
        fail(err.message);
    });
})();

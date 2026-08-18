# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

CoGe (Comparative Genomics) — a Perl web platform for comparative genomics (genomevolution.org). It is a long-lived legacy codebase with three cooperating tiers:

1. **CGI web pages** (`web/*.pl`) served by Apache, rendering `HTML::Template` files from `tmpl/` and doing AJAX via `CGI::Ajax`.
2. **REST API** (`web/services/api.pl`) — a `Mojolicious::Lite` app; all routes are declared in that single file and dispatched to controllers in `modules/Services/lib/CoGe/Services/API/`.
3. **JEX** — an external job-execution server reached over ZeroMQ (`CoGe::JEX::Jex`, host/port from `JOBSERVER`/`JOBPORT` in the config). Long-running analyses are compiled into workflows and submitted to it; they are *not* run in the web process.

The current branch is `ubuntu22` (port of the stack to Ubuntu 22 / Perl 5.34); `master` is the upstream main branch.

## Configuration (read this before running anything)

Nothing works without `coge.conf` at the repo root. It is gitignored and not present in a fresh clone.

- `CoGe::Accessory::Web::get_defaults()` (`modules/Accessory/lib/CoGe/Accessory/Web.pm`) is the single config accessor. It reads `$COGE_HOME/coge.conf` — and **falls back to `/opt/apache2/coge/` when `COGE_HOME` is unset**, which is the usual cause of confusing "can't read config" failures.
- The file is a flat `NAME value` (whitespace-separated) list, `#` comments. `BINDIR`, `COGEDIR`, `RESOURCEDIR`, `SCRIPTDIR`, `TMPLDIR` are derived from `COGE_HOME` and must not be set in the file. `setup.sh` shows a representative key set.
- Keys are read all over the code with no schema or validation; a missing key is usually a silent undef (except `PICARD`, which `CoGe::Builder::Buildable::BUILD` throws on).

## Commands

Install/refresh the Perl modules (each `modules/<Name>/` is a separate MakeMaker dist):

```bash
./make_perl.sh                    # installs all modules into the system perl lib (hardcoded to 5.34 x86_64-linux-gnu; runs twice on purpose)
cd modules/Core && perl Makefile.PL && make install    # single module
```

Because modules are *installed* rather than used from the tree, editing a file under `modules/` has no effect on the running app until it is reinstalled.

Run the API in development (morbo, auto-reload):

```bash
./api.sh start|stop|restart       # sets COGE_HOME=$PWD, port from MOJOLICIOUS_PORT in coge.conf (default 3303)
```

API integration tests (Robot Framework + RequestsLibrary; these are the maintained tests):

```bash
cd tests/api && ./run_tests.sh                  # all suites, output to web/test
robot --outputdir ../../web/test genomes.robot   # one suite
robot -t "Genome Fetch" genomes.robot            # one test case
robot -i public genomes.robot                    # only tests tagged public (vs. auth-required)
```

`tests/api/resource.robot` holds the shared variables — note `${API_URL}` points at an external production host by default, so override it to test a sandbox.

Front-end vendor libraries:

```bash
bower install                     # installs into web/js/vendor per .bowerrc
```

The `modules/*/t/*.t` files are mostly `use_ok` load checks plus a few DB tests with hardcoded, long-dead connection strings. Treat them as unmaintained; don't assume a failure there is a regression.

Perl formatting: `.perltidyrc` exists (80 cols, 4-space indent). Note it sets `--backup-and-modify-in-place`, so `perltidy` rewrites files and leaves `.bak` files behind.

## Architecture

### Job submission pipeline

Every analysis, load, and export goes through the same chain, which is the most important thing to understand in this repo:

```
payload {type, parameters, requester?}
  → CoGe::Factory::RequestFactory   maps type → CoGe::Request::* subclass (+ authRequired flag)
  → $request->is_valid / has_access validation and permission check
  → CoGe::Factory::PipelineFactory  maps type → CoGe::Builder::* subclass
  → $builder->pre_build / build / post_build   (extract / transform+load / notify)
  → $builder->submit()              adds tasks to a CoGe::JEX::Workflow, submits over ZeroMQ
```

`modules/Services/lib/CoGe/Services/API/Job.pm` `add()` is the canonical driver of that chain; `PipelineFactory.pm` and `RequestFactory.pm` each hold a `%typeToClass` table that must be kept in sync when adding a new job type. Web pages submit jobs by POSTing the same payload to the API rather than running work inline (`requester.page` is carried through only for logging and tiny-link generation).

`CoGe::Builder::Buildable` is the base class: subclasses override `build()` and push task hashrefs (`cmd`, `args`, `inputs`, `outputs`, `description`) onto `$self->tasks`. JEX resolves inter-task ordering from the declared inputs/outputs — tasks are not sequenced by insertion order. `CoGe::Builder::CommonTasks` supplies reusable task constructors (GFF/BED generation, masking, results registration).

`pre_build()` in `Buildable` creates the workflow with `init => 1` (new ID immediately). **SynMap and SynMap3D override `pre_build` with `init => 0`**, which makes JEX reuse an existing identical workflow and gives them results directories under `DIAGSDIR` instead of the standard staging/results layout — the two paths behave differently and code that assumes one often breaks the other.

### Module layout (`modules/`)

- `Database` — `CoGeX` is the DBIx::Class schema (`CoGeX/Result/*` rows, `CoGeX/ResultSet/*` custom resultsets); `CoGeDBI` holds hand-written SQL for hot paths that were too slow through the ORM.
- `Core` — business logic over the schema (`Genome`, `Experiment`, `Feature`, `Notebook`, `Metadata`, `Search`, `Storage`).
- `Core/Storage.pm` — all on-disk data paths. Genome and experiment data live in a **tiered directory scheme** (`get_tiered_path`: id split into thousands, e.g. `0/0/16/16911`) under `SEQDIR`; workflow staging/results live under `SECTEMPDIR/{staging,results}/<username>/<workflow_id>`.
- `Pipelines` — the `CoGe::Builder::*`, `CoGe::Request::*`, and factory classes described above.
- `Services` — Mojolicious controllers. `API/JBrowse/*` serves the embedded JBrowse genome browser; the router registers both namespaces, so a bare controller name can resolve into either — routes that must hit the non-JBrowse controller pass `namespace => 'CoGe::Services::API'` explicitly.
- `Services/Auth.pm` — accepts a session cookie, `username`+`token` query params, or a JWT header (`x-iplant-de-jwt`, `x-coge-jwt`), and returns `($db, $user, $conf)`. Every controller calls it itself; there is no global before_dispatch hook.
- `Accessory` — `Web.pm` (config, CAS login, CGI page `init`, tiny links, logging) plus parsers for external tool output (BLAST, blastz, LAGAN, dialign, GenBank).
- `Graphics` — GD-based genome/feature image rendering for the older viewers.
- `JEX` — ZeroMQ client and workflow model.
- `Exception` — Throwable-based exceptions; the API turns them into JSON via the `before_render` hook in `api.pl`.
- `Algos`, `ECNCS` — algorithm helpers (Codeml, Ks, PopGen) and a separate legacy ECNCS database.

### CGI page conventions

A page like `web/SynMap.pl` follows a fixed shape: `CoGe::Accessory::Web->init(cgi => ..., page_title => ...)` returns `($db, $user, $config, $tiny_link)`; a `%FUNCTIONS` table maps AJAX function names to subs and is merged with `CoGe::Accessory::Web::ajax_func()`; requests with a `jquery_ajax` param dispatch directly through that table, otherwise `CGI::Ajax->build_html` renders `gen_html`. Pages nest templates: `tmpl/generic_page.tmpl` (chrome) wraps a page template, which pulls fragments from `tmpl/partials/`. Per-page JavaScript lives in `web/js/pages/<page>.js`, shared widgets in `web/js/coge/`.

### Other directories

- `bin/` — bundled third-party binaries and wrappers invoked by pipeline tasks (dagchainer, lagan, blastz, quota-alignment, Tabix, codeml…). Paths to them come from config keys, not from hardcoded relative paths.
- `scripts/` — production scripts called from Apache/JEX, cron jobs, and manual utilities, mixed together; `scripts/README` says which subdirectories are live and which are obsolete (`scripts/old/` is dead).
- `old/` — retired pages and templates kept for reference; do not extend them.
- `web/js/jbrowse/` — JBrowse install (gitignored except `plugins/CoGe`, vendored directly since 2026-08-18 — it was a submodule of the now-dead `LyonsLab/CoGe_plugin`).
- `web/services/jex.py` — small Python WSGI shim that proxies JBrowse/job-status requests to JEX.

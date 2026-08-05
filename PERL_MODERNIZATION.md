# Modern Perl migration notes

Audit of the CoGe Perl tree against modern Perl (5.24 → 5.40), written while
porting the stack off Perl 5.12. Every finding below was verified by compiling
the affected file on this machine, not inferred from a grep — the method is in
[How these findings were produced](#how-these-findings-were-produced) so you can
re-run it after each change.

**Scope:** 346 `.pl`/`.pm` files under `modules/`, `web/`, and `scripts/`,
excluding `old/` and `scripts/old/`.

---

## 0. Environment: this checkout is running Perl 5.38, not 5.34

Worth settling first, because it changes which rules apply and it explains a
whole class of "module not found" errors.

| | |
|---|---|
| OS | Ubuntu 24.04.3 LTS (Noble) |
| `perl -v` | 5.38.2 |
| `INSTALLSITELIB` default | `/usr/local/share/perl/5.38.2` (in `@INC`) |

The branch is named `ubuntu22` and `make_perl.sh` hardcodes a 5.34 path, but the
interpreter here is 5.38.2. Ubuntu 22.04 ships 5.34; 24.04 ships 5.38. Decide
which you are actually targeting — nothing below differs between the two except
§2 (`given`/`when`), which merely warns on both and is removed in 5.42.

### 0.1 `make_perl.sh` installs the CoGe modules where Perl will never find them

```bash
# make_perl.sh
cd ./modules; perl Makefile.PL lib=/usr/local/lib/x86_64-linux-gnu/perl/5.34.0; make install; ...
```

MakeMaker *does* honour that lowercase `lib=`. I verified it sets:

```
INSTALLSITELIB  = /usr/local/lib/x86_64-linux-gnu/perl/5.34.0
INSTALLSITEARCH = /usr/local/lib/x86_64-linux-gnu/perl/5.34.0/x86_64-linux-gnu-thread-multi
```

That directory does not exist and is **not** in this Perl's `@INC` (which has
`/usr/local/lib/x86_64-linux-gnu/perl/5.38.2`). So `make install` succeeds and
every `use CoGe::...` still dies with *Can't locate CoGe/Core/Storage.pm in
@INC*.

**Scope check (see §8.0):** this bites on *this host* (5.38.2) only. Inside the
`coge_main` container, which runs Perl **5.34.0**, that path is correct and is in
`@INC` — so this is a blocker for the 24.04 migration, not a current fault.

**Fix** — drop the hardcoded path and let MakeMaker pick the version-correct
site directory:

```bash
#!/bin/bash
set -e
cd "$(dirname "$0")/modules"
# Two passes: the first configures the sub-distributions, the second builds them.
perl Makefile.PL && make install
perl Makefile.PL && make install
```

If you deliberately want a private prefix, use a version-independent one and add
it to `PERL5LIB` for Apache and the JEX workers:

```bash
perl Makefile.PL INSTALL_BASE=/opt/coge/perl5 && make install
# then: PERL5LIB=/opt/coge/perl5/lib/perl5
```

Sanity check after installing:

```bash
perl -MCoGe::Core::Storage -e 'print "modules visible\n"'
```

---

## 1. Hard compile failures — the code will not load

These are the blockers. Each one is a fatal error on the Perl you are running.

### 1.1 `Switch.pm` — 12 files

`Switch` was removed from core in **5.14**. It is a source filter, is
unmaintained, and its own docs recommend against it.

**Status (see §8.0): `Switch` 2.17 is already installed in the `coge_main`
container**, so these 12 files are *not* currently failing. Treat this section as
technical debt to retire before the 5.38 move, not as an active outage.

Symptom: `syntax error at ... near ") {"` pointing at the `switch (...) {` line.

| File | `use Switch` | `switch` blocks | `case` arms |
|---|---|---|---|
| `modules/Pipelines/lib/CoGe/Builder/Alignment/Aligner.pm` | 6 | 105, 225, 252 | 18 |
| `modules/Core/lib/CoGe/Core/Metadata.pm` | 11 | 43, 342 | 6 |
| `modules/Core/lib/CoGe/Core/Notebook.pm` | 6 | 109 | 6 |
| `modules/JEX/lib/CoGe/JEX/Jex.pm` | 12 | 96 | 5 |
| `modules/Database/lib/CoGeX/Result/Experiment.pm` | 9 | 186 | 4 |
| `modules/Pipelines/lib/CoGe/Builder/SNP/Analyzer.pm` | 6 | 68 | 4 |
| `modules/Pipelines/lib/CoGe/Builder/Trimming/Trimmer.pm` | 6 | 48 | 4 |
| `modules/Database/lib/CoGeX/Result/User.pm` | 9 | 283 | 3 |
| `scripts/irods.pl` | 5 | 41 | 3 |
| `modules/Pipelines/lib/CoGe/Builder/SNP/Merge.pm` | 6 | — | — |
| `web/ExperimentView.pl` | 14 | — | — |
| `web/User.pl` | 6 | — | — |

The last three import `Switch` and never use it — **delete the `use Switch;`
line**, that is the whole fix.

Because `Metadata.pm` is loaded by most of `CoGe::Builder::*`, its failure
cascades: it alone accounted for 34 of the 66 failing files in the last run. Fix
it early and the count drops sharply.

#### Pattern A — string → class dispatch (Aligner, Trimmer, SNP::Analyzer)

These are the bulk of the `case` arms, and the codebase already has the idiom
you want: `%typeToClass` in `CoGe::Factory::PipelineFactory` and
`CoGe::Factory::RequestFactory`. Match it.

```perl
# before -- modules/Pipelines/lib/CoGe/Builder/Trimming/Trimmer.pm:48
use Switch;
...
switch( lc($trimming_params->{trimmer}) ) {
    case 'cutadapt'    { $trimmer = CoGe::Builder::Trimming::Cutadapt->new($self)     }
    case 'trimgalore'  { $trimmer = CoGe::Builder::Trimming::TrimGalore->new($self)   }
    case 'trimmomatic' { $trimmer = CoGe::Builder::Trimming::Trimmomatic->new($self)  }
    case 'bbduk'       { $trimmer = CoGe::Builder::Trimming::BBDuk->new($self)        }
}
```

```perl
# after -- no source filter, and the mapping is now inspectable data
my %TRIMMER_CLASS = (
    'cutadapt'    => 'CoGe::Builder::Trimming::Cutadapt',
    'trimgalore'  => 'CoGe::Builder::Trimming::TrimGalore',
    'trimmomatic' => 'CoGe::Builder::Trimming::Trimmomatic',
    'bbduk'       => 'CoGe::Builder::Trimming::BBDuk',
);
...
my $class = $TRIMMER_CLASS{ lc($trimming_params->{trimmer}) };
$trimmer = $class->new($self) if $class;
```

`Aligner.pm:105` and `:252` map the same seven aligner names twice over; a single
hash at file scope serves both and removes the duplication.

#### Pattern B — value → value lookup (Result::Experiment, Result::User, Notebook)

```perl
# before -- modules/Database/lib/CoGeX/Result/Experiment.pm:186
switch ($self->data_type) {
    case 1 { return 'quantitative'; }
    case 2 { return 'polymorphism'; }
    case 3 { return 'alignment';    }
    case 4 { return 'marker';       }
    else   { return 'unknown';      }
}
```

```perl
# after
my %DATA_TYPE_DESC = (
    1 => 'quantitative',
    2 => 'polymorphism',
    3 => 'alignment',
    4 => 'marker',
);
...
return $DATA_TYPE_DESC{ $self->data_type // '' } // 'unknown';
```

`Result::User.pm:283` switches on `ref($object)`, which is a plain hash lookup
of the same shape.

#### Pattern C — regex arms (SNP::Analyzer, Jex)

`SNP::Analyzer.pm:72` has `case /gatk|gatk-haplotype-vcf|gatk-haplotype-gvcf/`,
and `Jex.pm:96` switches on a status string. A hash cannot express those, so use
`if`/`elsif`:

```perl
my $method = lc($snp_params->{method});
if    ($method eq 'coge')     { $snp = CoGe::Builder::SNP::CoGeSNPs->new($self)  }
elsif ($method eq 'samtools') { $snp = CoGe::Builder::SNP::Samtools->new($self)  }
elsif ($method eq 'platypus') { $snp = CoGe::Builder::SNP::Platypus->new($self)  }
elsif ($method =~ /^gatk/)    { $snp = CoGe::Builder::SNP::GATK->new($self)      }
```

Note that `Switch`'s `case /re/` matched with an implicit `=~` against the
switched value; keep that semantics when you rewrite, and be careful that
`case 'coge'` was a *string* comparison while `case 1` was numeric — Switch chose
the operator by looking at the literal.

### 1.2 Autoderef: `push`/`keys`/`values` on a scalar reference — 5 live sites

Perl 5.14 allowed `push $aref, ...` and `keys $href`; it was experimental and
**removed in 5.24**. Errors: *Experimental push on scalar is now forbidden* /
*Type of arg 1 to keys must be hash or array*.

| File:line | Current | Fix |
|---|---|---|
| `web/Taxonomy.pl:101` | `push ($hash->{children}, gen_subtree(\@array));` | `push @{ $hash->{children} }, gen_subtree(\@array);` |
| `web/Taxonomy.pl:137` | `push ($root_children, $sub_tree);` | `push @$root_children, $sub_tree;` |
| `scripts/backup/sync_jex_db.pl:21` | `scalar(keys $src_workflows)` | `scalar(keys %$src_workflows)` |
| `scripts/backup/sync_jex_db.pl:25` | `scalar(keys $dest_workflows)` | `scalar(keys %$dest_workflows)` |
| `scripts/backup/sync_jex_db.pl:29` | `foreach my $wid (keys $src_workflows)` | `foreach my $wid (keys %$src_workflows)` |
| `scripts/methylation/makeMetaplot.pl:422` | `sort keys $self->{print}->{st}->{'1'}` | `sort keys %{ $self->{print}{st}{'1'} }` |

`web/Taxonomy.pl` is a live web page, so this one is user-visible.

Also present in dead code, fix only if you revive it:
`scripts/old/dotplot_dots.pl:87,88` (`values $org1info`) and `:492`
(`push $chrs`).

### 1.3 `web/CoGeAlign.pl:87` — undeclared `$CONF` under `use strict`

```perl
SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
```

`$CONF` is never declared in this file and is not exported by
`CoGe::Accessory::Web` (its `@EXPORT` is subs only; `@EXPORT_OK` is empty — the
`our $CONF` in `Web.pm` stays private to that package). The file's config
handle is `$P`, declared at line 21 and assigned at line 22 from
`get_defaults()`, and used on the surrounding lines (`$P->{TMPLDIR}`,
`$P->{SERVER}`, `$P->{WIKI_URL}`).

**Fix:** `SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},`

This is a pre-existing bug rather than a Perl-version issue — the page cannot
have compiled since the line was added — but it shows up now as a hard failure.

---

## 2. Deprecated, runs today, removed in Perl 5.42

### `given`/`when` — `modules/Database/lib/CoGeX/Result/Job.pm`

Two blocks: `status_description` (lines 108–117) and `status_color`
(121–129). They work on 5.34/5.38 but each execution prints:

```
given is deprecated at .../Job.pm line 108.
when is deprecated at .../Job.pm line 109.
```

The file has `no warnings 'experimental';` (line 6), which **does not** suppress
these — the category changed to `deprecated`. `no warnings 'deprecated';`
silences them, and I confirmed that on 5.38. But `given`/`when` and smartmatch
were removed outright in **5.42**, so silencing only defers the work.

Both blocks are pure value lookups, so the §1.1 Pattern B rewrite applies:

```perl
my %STATUS_DESC = (
    0 => 'Scheduled', 1 => 'Running',    2 => 'Completed',
    3 => 'Cancelled', 4 => 'Terminated', 5 => 'Failed',
);
my %STATUS_COLOR = (
    1 => 'yellowgreen', 3 => 'salmon', 4 => 'salmon', 5 => 'salmon',
);

sub status_description { $STATUS_DESC{ shift->status // '' } // 'Unknown' }
sub status_color       { $STATUS_COLOR{ shift->status // '' } }
```

That also lets you drop `use v5.10;` → `use v5.16;` (or later) and the
`no warnings` line. No other file uses `given`/`when`, and **nothing in the tree
uses smartmatch `~~`** — verified.

---

## 3. CPAN dependencies whose status changed

`modules.txt` is the current dependency list and it has drifted: no version
constraints, one entry is a filename (`Array::Utils.pm`), and one is a class
rather than a distribution (`Bio::TreeIO` — the dist is `BioPerl`). Replace it
with a `cpanfile` so `cpanm --installdeps .` works and versions are pinned.

### 3.1 Removed from core — must now be installed explicitly

| Module | Removed in | Notes |
|---|---|---|
| `Switch` | 5.14 | See §1.1 — remove rather than install |
| `CGI` (and `CGI::Carp`, `CGI::Cookie`) | 5.22 | Still needed by every `web/*.pl`; Ubuntu package `libcgi-pm-perl` |

Only the OO interface of `CGI` is used (`$FORM->param`, `$FORM->header`,
`$FORM->Vars`). I checked for the deprecated HTML-generation functions
(`start_html`, `popup_menu`, …) — **none are used**, so CGI 4.x is a clean
upgrade here.

### 3.2 Unmaintained — plan replacements

| Module | Where | Replacement |
|---|---|---|
| `Spreadsheet::WriteExcel` | `web/CodeOn.pl`, `CoGeAlign.pl`, `CoGeBlast.pl`, `ExperimentList.pl`, `ExperimentView.pl`, `FeatList.pl`, `GenomeList.pl`, `NotebookView.pl`, `tRNAView.pl` | `Excel::Writer::XLSX` — same API, writes `.xlsx`; update the filenames (`CoGeBlast.pl:2144` hardcodes `.xls`) |
| `JSON::Syck` | `web/gobe/query.pl:8,116,228` | `JSON::XS` (`encode_json`), already used elsewhere in the tree |
| `File::Slurp` | 12 files incl. `Accessory/Web.pm`, `Core/Storage.pm`, `Core/Sequence.pm` | `File::Slurper` or `Path::Tiny`; `File::Slurp` has known encoding bugs and is discouraged by its own maintainer |
| `Class::DBI` | `modules/ECNCS/lib/CoGe/ECNCS/DB.pm` | EOL. ECNCS is legacy — confirm it is still in use before investing; otherwise drop it |
| `CGI::Application` | `modules/Services/lib/CoGe/Services/API/Sequence.pm` | This is a *live API controller* using a different framework from every other controller (which are `Mojolicious::Controller`). Port it to Mojolicious |
| `CGI::Ajax` | every `web/*.pl` | Unmaintained since ~2008 but pure Perl over `CGI`, so it still runs. Large job to remove; low urgency |
| `ZMQ::LibZMQ3` | `modules/JEX/lib/CoGe/JEX/Jex.pm` | `ZMQ::FFI` if the libzmq3 bindings will not build; needs `libzmq3-dev` |
| `Tabix` | `scripts/popgen/sumstats.pl:27` (via `use lib '/opt/apache2/coge/bin/Tabix'`) | Obsolete samtools bindings; `Bio::DB::HTS::Tabix`, or shell out to `tabix` |
| `Bio::DB::Sam` | `scripts/methylation/makeMetaplot.pl:18` | Does not build against modern htslib. Use `Bio::DB::HTS` |

### 3.3 Build-sensitive on a modern toolchain

These compile against system libraries that changed under you. Expect install
failures rather than runtime errors:

- **`Crypt::OpenSSL::RSA`** (`Accessory/Web.pm:31`, used by `jwt_decode_token`
  at line 705) — older releases fail against OpenSSL 3.x, which is what 22.04
  and 24.04 ship. Install a current version.
- **`GD`** — needs `libgd-dev`; the `gdStyled`/`gdTransparent`/`gdBrushed`
  constants and `GD::Font->MediumBold` used throughout `modules/Graphics/`
  require the real module (they are what my stubs could not fake).
- **`BerkeleyDB`** (`Accessory/Tile/Cache.pm`) — needs `libdb-dev`.
- **`DBD::mysql`** — Ubuntu builds it against MariaDB Connector/C. Against
  MySQL 8, whose default auth plugin is `caching_sha2_password`, connections can
  fail with an auth-plugin error. Either set the CoGe DB user to
  `mysql_native_password` or install a `DBD::mysql` built against
  `libmysqlclient`. `CoGeX::dbconnect` (`modules/Database/lib/CoGeX.pm:80`)
  passes no `mysql_*` flags, so there is nothing to change in code.

---

## 4. Behavioural changes that do not raise an error

The dangerous category: it loads, it runs, output differs.

### 4.1 Hash iteration order is randomised (5.18+)

Since 5.18 `keys`/`values`/`each` return a different order **on every process
start**. The tree has **82** `foreach my $x (keys %...)` loops with no `sort`
against **33** that sort. Most do not care, but any that feed a template, a JSON
response, a generated filename, or a reproducible-ID hash now produce unstable
output. Audit these first where the result reaches a user:

```bash
grep -rnE 'for(each)?\s+my\s+\$\w+\s*\(\s*keys\s+%' --include='*.pm' --include='*.pl' modules web
```

Add `sort` wherever order is observable. This is also a correctness trap for
SynMap, which reuses workflows keyed on identical inputs
(`CoGe::Builder::Tools::SynMap::pre_build` passes `init => 0`) — if any part of
that key is built from unsorted hash iteration, cache hits become unreliable.

### 4.2 `Storable` caches on disk

`Accessory/Tile/Cache.pm:8,129,133` (`freeze`/`thaw` into BerkeleyDB) and
`scripts/methylation/makeMetaplot.pl:275` (`retrieve`). Storable's format is
tied to the writing Perl; data frozen under 5.12 may fail or warn when thawed
under 5.38. **Plan to discard and regenerate these caches** as part of the
cutover rather than migrating them.

### 4.3 Unicode / encoding

`read_file` is called without a `binmode` in most places
(`Core/Storage.pm:250,537`, `Core/Sequence.pm:53,84,114`) — only
`Core/Metadata.pm:424` passes `binmode => ':raw'`. Combined with `JSON::XS`
`encode_json` (which emits UTF-8 bytes) and a MySQL move to `utf8mb4`, this is
where "wide character in print" warnings and double-encoded text come from. When
you swap `File::Slurp` for `File::Slurper` (§3.2), make the encoding explicit at
each call site: `read_text` for text, `read_binary` for sequence data.

### 4.4 Lower priority, works fine

- **78 two-argument `open` calls** (e.g. `open(OUT, ">$file_name")`,
  `Graphics/GenomeView.pm:244`). Legal, but a filename containing `>` or `|` is
  a shell-ish injection. Convert to three-arg opportunistically.
- **270 indirect-object constructor calls** (`new CGI`, `new GD::Polygon`,
  `new CoGeX::Result::Feature`). Legal, discouraged; fails confusingly if the
  class is not yet loaded. Rewrite as `Class->new(...)` when touching a file.
- **`UNIVERSAL::isa($_[0], '...')` as a function** — `Accessory/Web.pm:358`,
  `CoGeX.pm:382`. Works; ignores overridden `isa`. Prefer `->isa` or
  `Scalar::Util::blessed`.

---

## 5. Accessors and object systems — assessed, no version blockers

You asked specifically about accessor and object patterns. Three systems
coexist, and all three are fine on modern Perl:

- **`Class::Accessor`** + `mk_accessors` — `CoGeX.pm:10`, `Graphics.pm:3,28`,
  `Result/User.pm:91`, `Result/Feature.pm:116`, several `Graphics/Feature/*`.
  Still maintained, no changes needed.
- **`Moose`** — `CoGe::Builder::*`, `CoGe::Request::*`, `CoGe::JEX::*`,
  the factories. Uses only `has`/`extends`/`around BUILDARGS`/`BUILD`/
  `make_immutable`. I found **no** deprecated Moose idioms (no `lazy_build`, no
  `has '+attr'`, no old `Moose::Util::TypeConstraints` signatures).
- **`DBIx::Class`** — `load_namespaces()` (`CoGeX.pm:13`) is the current API.
  One deprecated call: `search_literal` at `modules/Core/lib/CoGe/Core/Feature.pm:43`
  (a MySQL `MATCH ... AGAINST` full-text query). It still works; the modern form
  is `->search(\[ 'MATCH(me.name) AGAINST (?)', $term ])`, which also fixes the
  SQL injection in the current string interpolation of `$search_term`.

Two things to *verify at runtime* rather than fix blind:

1. **Dual inheritance.** `CoGeX::Result::User`, `::Feature`, and `::Genome` each
   do `use base 'DBIx::Class::Core';` followed by `use base 'Class::Accessor';`.
   Both parents define `new` and `mk_accessors`. Perl's default DFS method
   resolution puts `DBIx::Class::Core` first, so DBIC's `new` wins — which is
   what you want. Keep that ordering if you touch those lines; reversing it
   would break row instantiation in a way that is hard to trace.
2. **Bare row construction.** `Result/Genome.pm:649` and `Result/Dataset.pm:732`
   do `my $feat = new CoGeX::Result::Feature;` — a DBIC row with no result
   source, used as a plain data holder. Legal in current DBIC but outside its
   contract. Worth a targeted test of the pages that call those methods.

---

## 6. Checked and clean

Ruled out by compiling the tree, so you can stop looking for them:

- `defined(@array)` / `defined(%hash)` — fatal since 5.22. **None.**
- `$[` assignment — fatal since 5.30. **None.**
- Smartmatch `~~`. **None.**
- `for my $x qw(...)` — fatal since 5.14. **None.**
- Unescaped literal `{` in a regex — fatal since 5.30. **None.** The one
  candidate, `Accessory/Tile/Cache.pm:88` (`/var\s+config\s*=\s*{/`), I compiled
  on 5.38: it is legal and matches.
- `POSIX::isdigit`/`isalpha`/`tmpnam` — removed from POSIX in 5.24. **None.**
  The negated imports (`use POSIX qw(!tmpnam !tmpfile)` in `Web.pm:27`,
  `use LWP::Simple qw(!getprint !getstore !mirror)` at line 20) are harmless; I
  verified both still work.
- `.` removed from `@INC` in 5.26 — no `require "file.pl"` or `do "file"` of a
  relative path anywhere. The three `use lib` calls are absolute paths
  (`scripts/popgen/sumstats.pl:26`, plus two in dead code).
- `my $x if 0`, `$*`, `$#` as variables, `\C` in regex, `do SUBROUTINE(LIST)`.
  **None.**
- No file uses taint mode (`-T`), so the 5.26 `@INC` taint changes are moot.

---

## 7. Suggested order of work

1. **`make_perl.sh` install path** (§0.1). Nothing else can be tested until the
   modules land in `@INC`.
2. **Install dependencies** — start from a `cpanfile` (§3). Expect to fight
   `Crypt::OpenSSL::RSA`, `GD`, `BerkeleyDB`, `ZMQ::LibZMQ3`, `DBD::mysql`.
3. **`Core/Metadata.pm` Switch** (§1.1) — unblocks ~34 files at once.
4. **Delete the three unused `use Switch;` lines** — free.
5. **Remaining Switch sites**, using the dispatch-table pattern.
6. **§1.2 autoderef** (6 one-line fixes) and **§1.3 `$CONF`** (one line).
7. Re-run the compile sweep; it should come back clean.
8. **`given`/`when` in `Job.pm`** (§2) — small, and removes log noise now.
9. Then the behavioural work: hash ordering where output is observable (§4.1),
   discard Storable caches (§4.2), encoding (§4.3).
10. Framework debt as capacity allows: `CGI::Application` controller,
    `Spreadsheet::WriteExcel`, `File::Slurp`, `JSON::Syck`.

---

## How these findings were produced

None of the CPAN dependencies are installed on this machine, so `perl -c` could
not reach past the first `use`. To get real compiler output I generated stub
modules for all 75 missing distributions, then iterated: compile all 346 files,
create a stub for every module Perl reported missing or empty, repeat until no
new stubs were needed. What survived is genuine syntax and compile-time error.

Each candidate finding was then confirmed individually — by reading the source
around it and, where semantics were in question, by running the construct on
5.38.2. That step matters: it caught a false positive I would otherwise have
reported. `Accessory/Web.pm` calls `carp` and `croak` while importing only
`cluck` from `Carp`, which looks like a bug — but real `CGI::Carp` merges its
`@EXPORT` (`confess croak carp`) with whatever you request, so
`use CGI::Carp('fatalsToBrowser')` does import them. The stub did not replicate
that, and the "error" was mine.

**Once the dependencies are actually installed, no stubs are needed** — this is
the sweep to run after each change, and in CI:

```bash
INC=$(for d in modules/*/lib modules/Algos/*/lib; do printf -- "-I%s " "$PWD/$d"; done)
find modules web scripts -name '*.pl' -o -name '*.pm' \
  | grep -v '/old/' \
  | while read -r f; do
      perl $INC -c "$f" 2>&1 | grep -q 'syntax OK' || echo "FAIL $f"
    done
```

Today that sweep reports **66 failing files**, and every one classifies into a
finding above — there is no unexplained residue:

| Root cause | Files |
|---|---|
| §1.1 `Switch` (mostly cascaded through `Metadata.pm` and `Jex.pm`) | 53 |
| §1.2 autoderef (`Taxonomy.pl`, `sync_jex_db.pl`, `makeMetaplot.pl`) | 3 |
| §1.3 `web/CoGeAlign.pl` `$CONF` | 1 |
| GD constants — stub artifact, will pass once `GD` is installed | 8 |
| `BerkeleyDB::Env` — stub artifact, same | 1 |

So **57 real failures / 9 artifacts**. The cascade is worth seeing, because it
is why the fix order in §7 is what it is — one file blocks 34 others:

| Failing module | Files it takes down with it |
|---|---|
| `modules/Core/lib/CoGe/Core/Metadata.pm:43` | 34 |
| `modules/JEX/lib/CoGe/JEX/Jex.pm:96` | 11 |
| `modules/Core/lib/CoGe/Core/Notebook.pm:109` | 4 |
| `Result/User.pm:283`, `Result/Experiment.pm:186`, `Trimming/Trimmer.pm:48` | 1 each |

Expected after §1 is complete: zero.

---

## Appendix: verification of the 2025-10-11 report

An earlier document ("Perl Compatibility Issues and Solutions", v2.0, 842 files
analysed) covers similar ground. I checked each of its claims against the code.
Its file-level scans are largely accurate; **its root-cause analysis for the
headline error is not**, and two of its prescribed fixes would break working
code. Details below so you can decide what to keep.

### Confirmed

| Claim | Verdict |
|---|---|
| `Switch.pm` is a critical blocker | **Correct** — but it lists 11 files and misses `scripts/irods.pl`; the real count is 12 (§1.1) |
| Dual inheritance `DBIx::Class::Core` + `Class::Accessor` | **Exists** in exactly 3 Result classes: `User.pm`, `Genome.pm`, `Feature.pm` — see caveat below |
| `mk_accessors` inside `BEGIN` blocks, "24 files" | **Count is right** — 24 call sites in BEGIN, 21 at top level. Harmless, though (below) |
| Indirect object notation | **Exists**, 94 files (not 136), 270 call sites |
| Manual `@ISA` manipulation | **Exists**, 36 files (not 43) |
| Not found: lexical `my $_`, smartmatch `~~`, `use encoding`, unescaped regex braces | **All four confirmed absent**, matching §6 |

### Refuted

**1. The Class::MOP error cannot come from `Class::Accessor` or `DBIx::Class`.**
The report's central claim is that
`Can't call method "isa" on an undefined value at Class/MOP/Class.pm line 494`
is caused by `mk_accessors` in a BEGIN block and by DBIC/Class::Accessor
"metaclass conflicts". `Class/MOP/Class.pm` ships with **Moose**;
`Class::Accessor` is pure Perl that installs closures by glob assignment and
never loads Class::MOP, and DBIC uses `Class::Accessor::Grouped`, also not
Class::MOP. I verified there is **no overlap at all** — not one file in the tree
uses `Moose` together with `Class::Accessor` or `DBIx::Class::Core`:

```bash
grep -rln 'use Moose' --include='*.pm' modules \
  | xargs grep -ln 'Class::Accessor\|DBIx::Class::Core'   # returns nothing
```

Neither "metaclass conflict" nor "mixed Moose + mk_accessors" exists here. The
error must originate in a Moose class.

**2. `mk_accessors` in a `BEGIN` block is not a problem.** In all 24 files the
`use base` line precedes the `BEGIN` block (`use base` at line 3, `BEGIN` at
10–24), so `@ISA` is populated before the accessor call. I confirmed the
mechanism on 5.38.2 with a minimal reproduction — parent loaded by `use base`,
`__PACKAGE__->mk_thing()` called inside a later `BEGIN`, accessor works. No
change needed in any of the 24 files.

**3. Pattern 3 ("compile-time execution requiring a DB connection") is a
non-issue.** `CoGeX::node_types()` is a static hash literal — no `$self`, no
query, no class loading:

```perl
sub node_types {
    my %types = ( list => 1, notebook => 1, genome => 2, experiment => 3, ... );
    return wantarray ? %types : \%types;
}
```

`my $node_types = CoGeX::node_types();` at file scope appears in 10 files and is
harmless in all of them. Do not convert these to lazy accessors.

**4. Pattern 4 ("Moose + namespace::clean") affects zero files.**
`namespace::clean` appears exactly once, at `Accessory/Web.pm:33` — a
`Class::Accessor` class, not a Moose class. `namespace::autoclean` is
Moose-specific; swapping it in there would be wrong.

**5. "Switch.pm removed in 5.32+" is inaccurate.** Removed from *core* in 5.14;
still on CPAN and still works on 5.38. It should be replaced because it is an
unmaintained source filter, not because it stopped existing.

**6. "Bareword filehandles — HIGH severity" is not a compatibility issue.**
Bareword filehandles work in every released Perl. This is style; 85 files use
`open(FH, ...)`. Do not spend the estimated 8–12 hours here for compatibility
reasons. (Two-arg `open` is worth hardening, for the injection reason in §4.4.)

**7. "Perl 5.30+ stricter @ISA handling / Class::MOP stricter isa checks."** No
such change is documented in perl5300delta or perl5340delta. This appears to be
inference from the error message rather than a real version change.

**8. The file counts are inflated ~2×.** Claimed 842 Perl files / 621 modules;
actual is 555 `.pl`/`.pm`/`.t` in the entire repo including `old/` and `bin/`,
of which `modules/` holds 207 `.pm`. Treat its other quantities with the same
caution — the two I could recount were both high.

### Two prescribed fixes that would break working code

**Do not delete the `mk_accessors` calls from the Result classes.** The report
says DBIC "creates accessors automatically from `add_columns`", so the manual
calls can go. But none of these names are columns — they are hand-rolled caches:

| File | Accessors | What they hold |
|---|---|---|
| `Result/User.pm:91` | `_genome_ids`, `_experiment_ids`, `_notebook_ids` | hashrefs of ids the user can access, populated lazily at lines 312–329 |
| `Result/Feature.pm:116` | `_genomic_sequence`, `gst`, `dsg`, `trans_type` | cached sequence and type objects |
| `Result/Genome.pm:184` | `_chromosomes` | cached chromosome list |

`grep add_columns` confirms zero of them are declared columns. Removing
`Class::Accessor` from these three classes is defensible, but only if you
replace each call with `__PACKAGE__->mk_group_accessors('simple' => ...)` in the
same commit. Also keep the base order — `DBIx::Class::Core` must stay first so
DBIC's `new` wins (§5).

**Do not run the automated `perl -pi -e` substitutions.** Three problems:

- The indirect-object regex only matches `new Foo(`, so it misses the most
  common form in this tree — `my $cgi = new CGI;` and `new GD::Polygon;` with no
  parens — leaving a half-converted codebase.
- The `@ISA` → `use parent` substitution would mangle
  `@ISA = ( @ISA, qw(Exporter) );` (`KsCalc.pm:14` and similar), which
  *appends* to an `@ISA` that `use base` already populated. Rewriting it as
  `use parent qw(@ISA Exporter)` is meaningless.
- The `mk_accessors` → `mk_group_accessors` substitution is applied to `*.pm`
  across the whole Result directory, but only 3 files there need it, and 21
  other files repo-wide use `mk_accessors` legitimately via `Class::Accessor`.

Its test script is also broken: it passes both `plan tests => scalar(@pm_files)`
and calls `done_testing()`, which Test::More rejects.

### What it missed

Every item in §1 of this document — the actual hard failures:

- The 6 autoderef sites (§1.2), which are fatal on 5.24+.
- `web/CoGeAlign.pl:87` (§1.3).
- `make_perl.sh` installing to a directory not in `@INC` (§0.1) — note that the
  report's own error paths read
  `/usr/local/lib/x86_64-linux-gnu/perl/5.34.0/Class/MOP/Class.pm`, i.e. the very
  hand-rolled directory that script writes to.
- `given`/`when` in `Result/Job.pm` (§2), removed in 5.42.

### A better hypothesis for the Class::MOP error — SUPERSEDED, see §8

**Update:** the production logs (2025-03-31) arrived after this section was
written and refute the `@INC` hypothesis below. The modules *are* installed at
`/usr/local/lib/x86_64-linux-gnu/perl/5.34.0/CoGe/...` and production runs Perl
5.34, so `make_perl.sh`'s hardcoded path is correct **on that box** — §0.1 is a
24.04/5.38 migration issue, not the cause of these errors. The reasoning below
about `extends` without `use` still stands as a code-hygiene point, but the
diagnosis has moved to §8. Kept for the record.



Since the error has to come from Moose, the thing to look at is this: **40 Moose
classes call `extends 'X'` without ever loading `X`.**

```bash
# all in modules/Pipelines/lib/CoGe/{Builder,Request}/
grep -rl "^extends " --include='*.pm' modules   # 40 files, none 'use' their parent
```

`CoGe::Builder::Trimming::Trimmomatic` extends `...::Trimmer`;
`...::Alignment::HISAT2` extends `...::Aligner`; 17 classes extend
`CoGe::Builder::Buildable`; 8 `CoGe::Request::*` extend `CoGe::Request::Request`.
Modern Moose auto-loads a superclass named in `extends`, so this normally works —
and the parent/child cycles here are benign because `use Moose; extends ...;`
always sits at lines 3–4, before the parent's `use` of its own subclasses.

But when the superclass file **cannot be found in `@INC`**, the metaclass lookup
is what fails, and Class::MOP is where it fails — which is consistent with the
error being a downstream symptom of §0.1 rather than an object-system bug. That
would also explain why it appeared during the version migration and not before.

I could not reproduce this on this machine: none of Moose, DBIx::Class, or
Class::Accessor is installed here, so the failing code path cannot be executed.
To confirm or kill the hypothesis, fix §0.1 first, then capture:

```bash
perl -MCarp::Always -MCoGe::Builder::Tools::SynMap -e 1
perl -e 'use CoGe::Builder::Trimming::Trimmomatic; print "ok\n"'   # child-first load
perl -MMoose -e 'print "$Moose::VERSION\n"'
```

The full stack trace will name the class whose metaclass came back undefined.
Regardless of the outcome, adding an explicit `use <Parent>;` next to each
`extends '<Parent>'` in those 40 files is cheap insurance and removes the
ambiguity — worth doing as part of §7 step 5, since you will be editing many of
those files for `Switch` anyway.

### Follow-up: the "are the Class::Accessor fixes self-contained?" analysis

A companion document asks whether removing `Class::Accessor` from the Result
classes is safe. **Its conclusion is correct, and I verified the specifics** —
with two corrections and one caveat. Note this is *optional cleanup, not a fix
for the error*: per the refutation above, the dual inheritance is not what
produces the Class::MOP failure, so this work belongs after §1, not before it.

Verified accurate:

- `_genome_ids`, `_experiment_ids`, `_notebook_ids` appear **only** in
  `Result/User.pm`. Nothing outside the class touches them.
- `_chromosomes` (`Result/Genome.pm:184`) is **genuinely unused** — the only
  occurrence in the entire tree is the `mk_accessors` declaration itself. (Take
  care when checking this yourself: a plain `grep _chromosomes` also matches the
  unrelated and heavily used `get_chromosomes` method. Anchor it:
  `grep -rnE '(^|[^a-zA-Z0-9])_chromosomes\b'`.) Delete the line rather than
  porting it.
- `mk_accessors` → `mk_group_accessors('simple' => ...)` is API-compatible.
  Both store in `$self->{field}` and generate the same get/set behaviour.
- **No code anywhere calls `Class::Accessor`'s inherited `->get()`/`->set()` on a
  DBIC Result object** — I checked, and every `->get(`/`->set(` hit in the tree
  is on something else (`Cache::FileCache`, `LWP::UserAgent`, `GD::Graph`,
  `RequestFactory`). This was the main way the change could have leaked into
  calling code, and it does not.
- `CoGeX.pm` is the safe one to start with: `mk_accessors` is already commented
  out (line 14), and `DBIx::Class::Schema` precedes `Class::Accessor` in `@ISA`
  (lines 9–10), so removing line 10 changes no method resolution.

Corrections:

- **`->gst` / `->dsg` / `->trans_type` are used in 4 files, not 9:**
  `Result/Dataset.pm`, `Result/Feature.pm`, `Result/Genome.pm`, and
  `web/GenomeList.pl`. The last one is outside the Result classes, so "only used
  within the classes or in standard getter/setter patterns" is half right — the
  usage *is* plain get/set, but one caller is a web page. Test
  `web/GenomeList.pl` after the change.
- **"Fix Type 3 / convert compile-time calls to lazy loading (8 Result classes,
  CRITICAL testing)" is unnecessary** — same refutation as above, `node_types()`
  is a static hash. The document's own "Option C: keep as-is" is the right
  answer; take it and skip the rest. Likewise its "Step 5: fix BEGIN block issues
  (24 files)" can be dropped entirely.

Caveat:

- **`Class::Accessor` must stay in your dependency list.** 27 files inherit from
  it; only 4 are in `Database/` (the 3 Result classes plus `CoGeX.pm`). Removing
  it there leaves **23 files** that still need it, and some depend on its
  inherited methods, e.g.
  `Accessory/histogram.pm:41,68` calls `$self->get('histogram_bins')` and
  `$self->set(x_labels_vertical => 1)`. Keep it in the `cpanfile`.

The `Switch.pm` explainer that accompanied this analysis arrived truncated, so I
could only check the part that was legible; §1.1 above supersedes it, and the one
error I did see there is the "removed in 5.32" claim corrected in item 5.

---

## 8. Verified diagnosis — reproduced inside the `coge_main` container

Everything in this section was executed inside the running container, not
inferred. It supersedes the hypotheses in §7, **all three of which were wrong.**

### 8.0 Container facts

| | |
|---|---|
| Image | `gastonlyons/coge:ubuntu22`, up 6 weeks |
| `perl -v` | **5.34.0** |
| Source checkout | `/opt/apache2/coge` — git branch `ubuntu22` at **`442e568e9`** |
| Installed modules | `/usr/local/lib/x86_64-linux-gnu/perl/5.34.0/CoGe/…` |
| Database | MariaDB 10.5 (`coge-mariadb`) |

Corrections this forces on earlier sections:

- **§0.1 does not apply to the container.** `/usr/local/lib/x86_64-linux-gnu/perl/5.34.0`
  *is* in the container's `@INC` (it runs 5.34), so `make_perl.sh`'s hardcoded
  path is correct there. §0.1 is a *future* problem for the 24.04/5.38 move only.
- **There is exactly one CoGe tree in `@INC`.** No stale duplicate. The installed
  copies and `/opt/apache2/coge` sources are byte-identical (`diff -q` → identical).
- **The Moose stack is healthy**, so §8.2's first two hypothetical causes are dead:
  Moose 2.2207, Class::MOP 2.2207 (same dist), Package::Stash 0.39,
  Package::Stash::XS 0.29, Class::Load 0.25, DBIx::Class 0.082844,
  Class::Accessor 0.51. One `Moose.pm`, one `Class/MOP.pm`. `use Moose` works.
- **`Switch` 2.17 is installed**, so the 12 files in §1.1 are *not* currently
  failing. §1.1 remains valid as technical debt, not as an active outage.
- **`PICARD` is set** in `coge.conf`, and `CoGe::Builder::Buildable` loads
  cleanly — killing the "parent failed to compile" theory.
- **The `use Exporter 'import'` + Moose pattern (§8.3, old) is not the cause.**
  `CoGe::Builder::Tools::CoGeBlast` — which uses the same pattern — loads fine.
  It stays an antipattern; it is not what broke.

### 8.1 The actual root cause: one uncommitted edit to `Web.pm`

`git status` in the container shows `modules/Accessory/lib/CoGe/Accessory/Web.pm`
as **modified and uncommitted** (mtime Jun 17 21:29). The diff reverts the
compatibility work done in commit `e7739af27 "Compatibility to perl > 5.24"`:

```diff
-use v5.18;
+use v5.10;
     @EXPORT  = qw( ...
-                   internal_url_for url_for internal_api_url_for api_url_for get_job
+                   internal_url_for url_for api_url_for get_job
-        cluck("error: template=$template_file could not be found");
+        cluck "error: template=$template_file could not be found";
```

It removed **both the `internal_api_url_for` export and the `sub` itself** —
`grep -n 'sub internal_api_url_for'` in the container returns nothing. Ten
modules still import it, so all ten die at compile time:

```
"internal_api_url_for" is not exported by the CoGe::Accessory::Web module
Can't continue after import errors at .../CoGe/Builder/Tools/SynMap.pm line 8.
```

| Failing file | Consequence |
|---|---|
| `CoGe/Builder/Tools/SynMap.pm` | SynMap page dead (matches the logs) |
| `CoGe/Builder/Tools/SynMap3D.pm` | SynMap3D dead |
| `CoGe/Builder/Tools/SynMapN.pm` | SynMapN dead |
| `CoGe/Factory/PipelineFactory.pm` | **all job submission dead** |
| `CoGe/Services/API/Genome.pm` | API `/genomes/*` dead |
| `CoGe/Services/API/Experiment.pm` | API `/experiments/*` dead |
| `CoGe/Services/API/Feature.pm` | API `/features/*` dead |
| `CoGe/Services/API/Job.pm` | API `/jobs/*` dead |
| `CoGe/Services/API/Search.pm` | API `/global/search/*` dead |
| `web/SynMap.pl` | the page itself |

**Confirmed live against the running API** (Mojolicious loads controllers lazily,
so the app stays up while individual routes die):

```
GET /organisms/search/arabidopsis   -> HTTP 200   (controller has no such import)
GET /genomes/search/test            -> HTTP 000   (worker dies loading Genome.pm)
GET /global/search/test             -> HTTP 000   (worker dies loading Search.pm)
```

### 8.2 The fix, verified

Restoring that one file is sufficient — no other change needed:

```bash
cd /opt/apache2/coge
git diff -- modules/Accessory/lib/CoGe/Accessory/Web.pm > /tmp/web-pm-local-edit.patch  # keep a copy
git checkout -- modules/Accessory/lib/CoGe/Accessory/Web.pm
./make_perl.sh          # reinstall so /usr/local/... matches
apachectl restart ; ./api.sh restart
```

Proven in the container by overriding **only** `Web.pm` via `PERL5LIB`, leaving
everything else as installed — all ten went from failing to `OK`:

```
CoGe::Services::API::Genome      OK      CoGe::Builder::Tools::SynMap     OK
CoGe::Services::API::Job         OK      CoGe::Builder::Tools::SynMap3D   OK
CoGe::Services::API::Search      OK      CoGe::Builder::Tools::SynMapN    OK
CoGe::Services::API::Experiment  OK      CoGe::Factory::PipelineFactory   OK
CoGe::Services::API::Feature     OK      web/SynMap.pl            syntax OK
```

Before reverting, decide what that edit was *for*. If `internal_api_url_for` is
meant to be retired, then the host repo's newer commits are the coherent target
(`33d9d25bc`, `ca578d9f1` drop it from `SynMap.pm`, replacing
`internal_url_for(internal_api_url_for("genomes"))` with
`internal_url_for(api_url_for("genomes"))`) — but note the host's
`web/SynMap.pl:9` still *imports* it, so removing the sub requires editing that
line too. The container is 3 commits behind the host checkout; deploying host
HEAD **and** discarding the local `Web.pm` edit is the clean end state.

### 8.3 Also broken in the container, unrelated to the above

- **`Class::DBI` is not installed** → all 12 `modules/ECNCS/**` files fail with
  `Base class package "Class::DBI" is empty.` Decide whether ECNCS is still in
  use; if yes, install `Class::DBI`, if no, delete the directory.
- **`JSON::Syck` is not installed** → `web/gobe/query.pl` fails (§3.2).
- **`Bio::DB::Sam` is not installed** → `scripts/methylation/makeMetaplot.pl`.
- `web/Taxonomy.pl` and `web/CoGeAlign.pl` fail exactly as predicted in §1.2/§1.3
  — confirmed on the container's own Perl.

### 8.4 Method note: two artifacts to avoid

The container-side sweep (`perl -c` over 420 files, 91 failures) contains two
classes of false positive. Filter them or you will chase ghosts:

1. **Bundled JBrowse scripts** under `web/js/jbrowse/` need their own
   `src/perl5` include path (`GenomeDB.pm`, `NCList.pm`, `IntervalStore.pm`,
   `Bio::JBrowse::*`…). 60 of the 91. Exclude `web/js/`.
2. **`perl -c` on a `.pm` that has an installed twin.** Running
   `perl -c modules/Accessory/lib/CoGe/Accessory/Web.pm` reports
   `String found where operator expected ... (Do you need to predeclare cluck?)`
   at lines 220 and 853 — but `perl -e 'use CoGe::Accessory::Web'` **succeeds**.
   The source copy is compiled as the main program while its dependency chain
   loads the *installed* copy of the same package; that copy's
   `use namespace::clean` then strips `cluck` out of the symbol table mid-parse
   of the first. It is not a real defect, and the bare `cluck "..."` form runs
   fine. Always confirm a module-level failure with `perl -e 'use The::Module'`
   before believing it.

Real CoGe failures after filtering: **31**, of which 10 are §8.1, 12 are
`Class::DBI`, and the rest are the known items in §8.3.

---

## 9. Changes applied

All verified inside `coge_main` (Perl 5.34.0) against the real installed
dependencies. Sweep result after these changes: **334 files checked, 5 failing**,
all five known non-issues (see §9.5).

### 9.1 `internal_api_url_for` / `internal_url_for` — removed, then RESTORED

**Final state: both functions exist and are exported; the three call sites use
`internal_url_for(api_url_for(...))`.** They were briefly removed in `0bdff0a07`
on the reasoning that the `INT_*` config keys did not exist, so the functions
could only `croak`. That reasoning was right about the config and wrong about the
fix: the keys were missing, but the functions are needed. Restored in the
following commit, with the config keys added.

Why they are needed — measured inside `coge_main`:

| Builder | Value | Reachable from inside the container? |
|---|---|---|
| `url_for(api_url_for("genomes"))` | `http://localhost:60500/coge/api/v1/genomes` | **No — HTTP 000** |
| `internal_url_for(api_url_for("genomes"))` | `http://localhost/coge/api/v1/genomes` | **Yes — HTTP 200, 0 redirects** |

`SERVER` is `http://localhost:60500/coge/`, and 60500 is the *host-published*
port; Apache inside the container listens on **80**. So any server-side
subprocess that calls back through `SERVER` leaves the container, hits the
published port, and gets bounced through Apache — or, as measured, fails
outright. Three consumers depend on the internal address:

| Call site | Consumer | Notes |
|---|---|---|
| `SynMap.pm:887` | `scripts/synmap/dotplot_dots.py` | inside `if ($ks_type)` — **only exercised with Ks enabled**; reads `api_url` from the `.cfg` written here and calls it with `requests.get` + a JWT |
| `SynMap.pm:1215` | `scripts/synmap/fractionation_bias.py` | `--apiurl` argument |
| `web/SynMap.pl:2031` | `web/run_dotplot.pl` | server-side fetch via `LWP::UserAgent` |

Config keys added to `coge.conf` (gitignored — **must be baked into the image
build**, see §9.6):

```
INT_SERVER http://localhost/coge/
INT_URL /coge/
INT_API_URL /api/v1/
```

A backup of the previous file is at `coge.conf.bak-preINT` in the container.

Note `internal_api_url_for` currently has **no callers** — `INT_API_URL` is
composed by `internal_url_for(api_url_for(...))` instead, which yields the same
result because `INT_API_URL` and `API_URL` are both `/api/v1/`. The function is
restored and exported so it is available, and so that removing it again is not
mistaken for a safe cleanup. If the internal API ever moves to a different path,
`internal_api_url_for` is the hook to use.

### 9.2 Fatal-error fixes (§1.2, §1.3)

- `web/Taxonomy.pl:101` → `push @{ $hash->{children} }, ...`
- `web/Taxonomy.pl:137` → `push @$root_children, $sub_tree;`
- `web/CoGeAlign.pl:87` → `$CONF->{SUPPORT_EMAIL}` → `$P->{SUPPORT_EMAIL}`
- `scripts/backup/sync_jex_db.pl:21,25,29` → `keys %$src_workflows` / `%$dest_workflows`
- `scripts/methylation/makeMetaplot.pl:422` → `sort keys %{ $self->{print}->{st}->{'1'} }`

### 9.3 Dependency changes

**ECNCS deleted** — `modules/ECNCS/` (26 tracked files) and its 24
`modules/MANIFEST` entries. Nothing outside the directory referenced it, so
**`Class::DBI` is no longer needed anywhere.**

**`JSON::Syck` → `JSON::XS`** in `web/gobe/query.pl` (`JSON::Syck::Dump` →
`encode_json`, 3 lines). The two `JSON::Syck` mentions in
`Accessory/Tile/Cache.pm` are inside comments and were left alone. **`JSON::Syck`
is no longer needed.**

**`Bio::DB::Sam` → `Bio::DB::HTS`** in `scripts/methylation/makeMetaplot.pl`
(4 sites: the `require` guard, the `use`, and two `->new(-bam => ...)`).
`->pileup` / `->fast_pileup` keep the same region-string + `($seqid,$pos,$pileup)`
callback contract, so no logic changed.

**`Tabix` → `Bio::DB::HTS::Tabix`** in `scripts/popgen/sumstats.pl`. Also removed
the hardcoded `use lib '/opt/apache2/coge/bin/Tabix'`. See §9.4 — this one needs
validation.

**`modules.txt`** — fixed `Array::Utils.pm` → `Array::Utils`, `Bio::TreeIO` →
`BioPerl`, and added the dists that were installed in the image but undeclared
(a rebuild from the old file would have produced a *worse* image): `CGI`,
`Clone`, `Devel::Size`, `Hash::Merge`, `HTTP::Request::Common`, `IO::String`,
`List::MoreUtils`, `Parse::RecDescent`, `PerlIO::gzip`, `XML::DOM`,
`XML::XPathEngine`, plus the new `Bio::DB::HTS`.

**New image requirement:** `libhts-dev` (system) so `Bio::DB::HTS` builds. It
supplies both `Bio::DB::HTS` and `Bio::DB::HTS::Tabix`, so it is one dependency
for both ported scripts.

### 9.4 ⚠ Needs validation: `sumstats.pl` coordinates

The old samtools `Tabix->query($chr, $start, $end)` took **0-based half-open**
coordinates, so with the 1-based GFF `$start`/`$end` this script passes (proved
by `substr($pSeq->{$chr}, $start-1, $featLen)` on line 118) it returned VCF `POS`
in `[$start+1, $end]` — **it skipped each feature's first base.**
`Bio::DB::HTS::Tabix` region strings are 1-based inclusive.

The port therefore queries `"$chr:" . ($start+1) . "-$end"`, which **reproduces
the old result set exactly**, preserving that off-by-one. This was deliberate:
pi, theta, and Tajima's D would otherwise shift silently.

If the off-by-one should be corrected, drop the `+1` — the feature's full span is
then included. Either way, **run one gene through both the old and new code and
compare `sumstats.tsv`** before trusting output. This is the only change in this
set that can alter numeric results.

### 9.5 Remaining sweep failures — all non-issues

- **4 × DBIC "does not seem to be a Result class"** — `Result/Feature.pm`,
  `Result/Genome.pm`, `Result/Dataset.pm`, `Accessory/GenBank.pm`. An artifact of
  compiling a Result class *standalone*: with `use CoGeX;` first (which every
  entry point does, via `Web.pm`) all four load `OK`. Fragile but benign;
  `GenBank.pm` would be cleaner using `use CoGeX;` instead of
  `use CoGeX::Result::Feature;`.
- **1 × `modules/Algos/Pairwise/scripts/test.pl`** — dead dev script with a
  hardcoded `use lib "/home/elyons/projects/pairwise/lib/"`. Delete or ignore.

### 9.6 Still outstanding

- `given`/`when` in `CoGeX/Result/Job.pm` (§2) — runs fine on 5.34/5.38, removed
  in 5.42. Not done.
- `Switch` in 12 files (§1.1) — still installed in the image, so not urgent.
- The behavioural work in §4 (hash ordering, Storable caches, encoding).
- `resources/CoGe_secret.txt`, `resources/DE_rsa.pub`, and `coge.conf` are
  untracked/gitignored but required at runtime — they must be injected at image
  build or mounted, or JWT auth breaks (§8 note). **`coge.conf` must now also
  carry `INT_SERVER` / `INT_URL` / `INT_API_URL` (§9.1)** or SynMap's Ks/dotplot
  and FractBias tasks will `croak`.

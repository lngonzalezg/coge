#!/usr/bin/perl
#
# mod_perl parent-interpreter preload.
#
# WHO RUNS THIS: the Apache *parent* process, as root, at config-parse time --
# before it forks any children and before privileges drop to www-data. It is pulled
# in by the PerlRequire directive in the vhost (coge-main.conf), and again on every
# `apachectl graceful` / restart.
#
# WHY IT EXISTS: web/*.pl are served by ModPerl::Registry. Without a preload, every
# Perl interpreter compiles the Moose + CoGe module graph lazily, on its first
# request, in whatever order that request happens to trigger. That graph contains a
# circular use:
#
#     modules/Pipelines/lib/CoGe/Builder/Buildable.pm      :19  use CoGe::Builder::Data::Extractor;
#     modules/Pipelines/lib/CoGe/Builder/Data/Extractor.pm :4   extends 'CoGe::Builder::Buildable';
#
# so load order is not deterministic across interpreters. A bad order can leave
# Class::MOP holding a metaclass with an undef package name -- Class::MOP::Package
# ->initialize() does not validate its argument the way Class::MOP::Class->initialize()
# does, so it stores the nameless metaclass silently and only dies much later, with
# the badly misleading:
#
#     Package::Stash->new must be passed the name of the package to access
#         at .../Class/MOP/Package.pm line 218
#     BEGIN failed--compilation aborted at .../CoGe/Builder/Tools/SynMap.pm line 3
#
# Compiling the graph here, once, in a fixed order, removes that nondeterminism.
# Under prefork it additionally means children inherit the compiled modules via
# copy-on-write instead of each compiling its own copy.
#
# NOTE: this makes the known-good load order deterministic. It does NOT fix the
# Buildable <-> Extractor cycle. That cycle currently survives only because
# Extractor.pm never calls make_immutable (a mutable Moose class resolves inherited
# attributes lazily); adding make_immutable there will break it regardless of this
# file. Fixing the cycle properly is a separate change.
#
# RULES FOR THIS FILE:
#   * Modules only. Never open a database handle, a ZeroMQ/JEX socket, or a file
#     handle here. The parent's descriptors are inherited by every forked child,
#     which then interleave on them. `use CoGeX;` is fine (it only loads the schema
#     class); CoGeX->connect(...) is not.
#   * Keep it idempotent. Apache parses its config twice at startup, so this can be
#     executed more than once per boot.
#   * It runs as root. Do not create files or directories here.
#
# TO REMOVE / ROLL BACK:
#   1. Delete or comment out the `PerlRequire /opt/apache2/coge/startup.pl` line in
#      /home/lgonzalez/coge_docker/coge-main.conf (bind-mounted to
#      /etc/apache2/sites-available/coge-main.conf on the host).
#   2. docker exec coge_main apachectl graceful
#   That is the whole rollback -- Apache goes straight back to lazy per-interpreter
#   loading. Nothing else references this file, so leaving it on disk unreferenced is
#   harmless. To remove it completely, also drop the startup.pl bind mount from
#   /home/lgonzalez/coge_docker/docker-compose.yaml and recreate the container.
#   Pre-change backups of all three edited files are in ~ as *.bak-20260806.

use strict;
use warnings;

# CoGe::Accessory::Web::get_defaults() falls back to /opt/apache2/coge/ when
# COGE_HOME is unset. That fallback happens to be correct here, but relying on it is
# how the confusing "can't read config" failures start. Be explicit.
BEGIN { $ENV{COGE_HOME} ||= '/opt/apache2/coge'; }

#-------------------------------------------------------------------------------
# The Moose graph. THIS SECTION IS THE FIX -- order matters, do not alphabetize.
#
# Moose first, so its own Exporter setup happens once here, at the top level of the
# parent, rather than inside a request-time BEGIN block. Then Buildable *before*
# Extractor: Extractor extends Buildable, and Buildable's own BEGIN reaches back for
# Extractor at line 19. Entering the cycle at Buildable is the order that is known to
# work today (`perl -c web/SynMap.pl` passes and the page serves), so pin that one.
#-------------------------------------------------------------------------------
use Moose ();
use CoGe::Builder::Buildable;
use CoGe::Builder::Data::Extractor;
use CoGe::Builder::Tools::SynMap;
use CoGe::Builder::Tools::CoGeBlast;
use CoGe::JEX::Jex;
use CoGe::JEX::Workflow;
use CoGe::Exception::Generic;

#-------------------------------------------------------------------------------
# Shared non-Moose web stack. Optional -- this section is a startup-latency and
# memory optimisation only, not part of the fix. If anything here misbehaves at boot,
# comment out this block first; the section above is the part that matters.
#-------------------------------------------------------------------------------
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Builder::CommonTasks;
use CGI;
use CGI::Ajax;
use HTML::Template;
use JSON::XS;

1;

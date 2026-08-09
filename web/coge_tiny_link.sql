-- CoGe internal link shortener -- replaces the external YOURLS service.
--
-- Apply to the MAIN CoGe database (the one named by DBNAME in coge.conf), e.g.
--
--   mysql -h <DBHOST> -u <DBUSER> -p coge < web/coge_tiny_link.sql
--
-- Idempotent: safe to re-run. No data migration is performed -- the replacement is a
-- clean start, so pre-existing YOURLS keys are deliberately NOT imported and old
-- /r/<key> links stop resolving (see todos/tiny_link_replacement.md).

--
-- Table structure for table `tiny_link`
--

CREATE TABLE IF NOT EXISTS `tiny_link` (
  -- The keyword IS the hash, not a surrogate id: lc(base36(md5(rel_url)))[0..11].
  -- Consequences that the application depends on:
  --   * deterministic -- the same URL yields the same key across processes, across
  --     time, and across a wiped/restored database (SynMap result files under
  --     DIAGSDIR are named after the key and outlive this table);
  --   * matches \w+ -- SynFind.pl extracts it with /(\w+)$/ and SynMap.pm uses it as
  --     a filename component, so no '-' may appear (hence base36, not base64url);
  --   * lowercase only -- keys become filenames, so a case-insensitive filesystem
  --     anywhere in the toolchain must not be able to collide two of them. ascii_bin
  --     keeps MySQL from folding case as well.
  `keyword` char(12) CHARACTER SET ascii COLLATE ascii_bin NOT NULL,
  -- RELATIVE to the SERVER config value, e.g. "SynMap.pl?dsgid1=..;dsgid2=..".
  -- Storing only the relative form makes off-site targets unmintable by construction
  -- (no open redirect), converges scheme/host aliases onto one key, and lets the
  -- redirect re-attach whatever SERVER is current after a domain or http->https move.
  `rel_url` text NOT NULL,
  `created_at` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`keyword`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

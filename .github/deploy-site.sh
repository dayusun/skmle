#!/usr/bin/env bash
# Build the pkgdown site and publish it to gh-pages.
#
# Not pkgdown::deploy_to_branch(): that builds into its own worktree and pushes
# in one step, leaving nowhere to drop CLAUDE.md.  pkgdown renders every *.md in
# the package root apart from a hardcoded list (README, LICENCE, NEWS,
# cran-comments and the issue templates), so the developer instructions would
# otherwise be published at /CLAUDE.html and indexed in the site search.
set -euo pipefail
cd "$(dirname "$0")/.."

Rscript -e '
  pkgdown::build_site(preview = FALSE)

  # Drop the developer instructions from the built site, the search index and
  # the sitemap.  Nothing in it is secret; it is internal guidance, and a
  # student following the reference index should not land on it.
  unlink("docs/CLAUDE.html")
  if (file.exists("docs/CLAUDE.md")) unlink("docs/CLAUDE.md")

  idx <- jsonlite::fromJSON("docs/search.json", simplifyVector = FALSE)
  paths <- vapply(idx, function(e) {
    if (is.null(e$path)) "" else as.character(e$path)[1]
  }, character(1))
  jsonlite::write_json(idx[!grepl("CLAUDE", paths, fixed = TRUE)],
                       "docs/search.json", auto_unbox = TRUE)

  sm <- readLines("docs/sitemap.xml", warn = FALSE)
  one <- paste(sm, collapse = "\n")
  one <- gsub("<url>(?:(?!</url>).)*CLAUDE(?:(?!</url>).)*</url>\\s*", "",
              one, perl = TRUE)
  writeLines(one, "docs/sitemap.xml")

  stopifnot(!file.exists("docs/CLAUDE.html"),
            !any(grepl("CLAUDE", readLines("docs/search.json", warn = FALSE))),
            !any(grepl("CLAUDE", readLines("docs/sitemap.xml", warn = FALSE))))
  cat("site built and cleaned\n")
'

WT=$(mktemp -d)
trap 'git worktree remove --force "$WT" >/dev/null 2>&1 || true' EXIT
git fetch origin gh-pages --quiet
git worktree add --quiet --detach "$WT" origin/gh-pages
rsync -a --delete --exclude='.git' docs/ "$WT"/
touch "$WT/.nojekyll"
git -C "$WT" add -A
if git -C "$WT" diff --cached --quiet; then
  echo "site unchanged"
else
  git -C "$WT" commit --quiet -m "deploy: $(git rev-parse HEAD)"
  git -C "$WT" push --quiet origin HEAD:gh-pages
  echo "deployed $(git -C "$WT" rev-parse --short HEAD)"
fi

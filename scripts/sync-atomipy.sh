#!/usr/bin/env bash
#
# sync-atomipy.sh — mirror the canonical atomipy package into each app's
# embedded copy so all three stay byte-identical.
#
# The two apps (atomipy-web-module, atomipy-topology-generator) deploy from
# their own nested git repos, so they must each ship a real, current copy of
# atomipy (not a symlink). Run this after editing canonical atomipy, then
# review and commit the embedded copy in each app.
#
# Usage:
#   scripts/sync-atomipy.sh          # sync canonical -> both apps
#   scripts/sync-atomipy.sh --check  # verify only; non-zero exit if drifted
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/atomipy"
APPS=(
  "$ROOT/atomipy-web-module/atomipy"
  "$ROOT/atomipy-topology-generator/atomipy"
)
# Caches / OS junk are never part of the package.
EXCLUDES=(--exclude=__pycache__ --exclude='*.pyc' --exclude=.DS_Store)

CHECK_ONLY=0
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

if [ ! -d "$SRC" ]; then
  echo "ERROR: canonical atomipy not found at $SRC" >&2
  exit 1
fi

status=0
for dest in "${APPS[@]}"; do
  app="$(basename "$(dirname "$dest")")"

  if [ "$CHECK_ONLY" -eq 1 ]; then
    if [ -d "$dest" ] && diff -rq "${EXCLUDES[@]}" "$SRC" "$dest" >/dev/null 2>&1; then
      echo "✓ $app: embedded atomipy is identical to canonical"
    else
      echo "✗ $app: embedded atomipy DIFFERS from canonical (run without --check to fix)" >&2
      status=1
    fi
    continue
  fi

  # Replace a stale dev symlink with a real directory if needed.
  [ -L "$dest" ] && rm "$dest"
  mkdir -p "$dest"
  rsync -a --delete "${EXCLUDES[@]}" "$SRC/" "$dest/"

  if diff -rq "${EXCLUDES[@]}" "$SRC" "$dest" >/dev/null; then
    echo "✓ $app: synced and identical to canonical"
  else
    echo "✗ $app: still differs after sync" >&2
    status=1
  fi
done

if [ "$CHECK_ONLY" -eq 0 ] && [ "$status" -eq 0 ]; then
  echo
  echo "Done. Review and commit the embedded copy in each app, e.g.:"
  echo "  (cd atomipy-web-module        && git add atomipy && git commit -m 'chore(vendor): sync atomipy' && git push)"
  echo "  (cd atomipy-topology-generator && git add atomipy && git commit -m 'chore(vendor): sync atomipy' && git push)"
fi

exit "$status"

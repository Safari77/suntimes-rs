#!/bin/sh
set -eu -o pipefail

if [ "$#" -lt 1 ]; then
    printf 'usage: %s <pattern>\n' "$0" >&2
    exit 2
fi

pattern="$1"
tab=$(printf '\t')

cargo metadata --format-version 1 \
  | jaq -r '[.packages[] | select(.name == "rustix")]
           | unique_by(.manifest_path)
           | .[]
           | "\(.version)\t\(.manifest_path | sub("/Cargo\\.toml$"; ""))"' \
  | while IFS="$tab" read -r ver dir; do
      [ -d "$dir/src" ] || continue
      out=$(cd -- "$dir" && rg -n -A4 --heading -- "doc\\(alias\\s*=.*${pattern}" src/) || true
      if [ -n "$out" ]; then
          printf '=== rustix-%s ===\n%s\n' "$ver" "$out"
      fi
    done

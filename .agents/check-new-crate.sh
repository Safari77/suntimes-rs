#!/usr/bin/env sh
set -eu

if [ "$#" -lt 1 ]; then
  printf 'Usage: %s <crate> [cargo add args...]\n' "$0" >&2
  exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT

cargo init --bin --name crate_vet --vcs none -- "$tmp" >/dev/null
cargo add --manifest-path "$tmp/Cargo.toml" "$@"

if [ -f deny.toml ]; then
  cargo deny --config deny.toml --manifest-path "$tmp/Cargo.toml" check advisories
else
  cargo deny --manifest-path "$tmp/Cargo.toml" check advisories
fi

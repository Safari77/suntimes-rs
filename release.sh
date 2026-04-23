#!/bin/sh
cargo test --release && cargo release patch --no-publish --execute

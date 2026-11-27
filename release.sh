#!/bin/sh
cargo test --release && cargo build --release && cargo deny check && cargo release patch --no-publish --execute

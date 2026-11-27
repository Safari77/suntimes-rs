# Rust Workspace Agent Instructions

**The commit is ready IFF `cargo test` passes without errors.**

## Adding New Crates

Before editing `Cargo.toml` or running `cargo add` for a new dependency, run `./.agents/check-new-crate.sh <crate>` without modifying any local files.

- If the last line of output is `advisories FAILED` (or the command exits non-zero), do NOT add the crate. Report the error and ask the user whether the selected crate is still suitable, whether to try a suggested alternative (which must also be checked), or if the user wants to pick a crate themselves.
- Only modify `Cargo.toml` after the check passes or the user explicitly approves the crate.

## Unix Syscalls (`rustix`)

Prefer safe `rustix` wrappers over `unsafe { libc::* }`, subprocesses, or `/proc` workarounds.

### Local API Lookup
- **Never** read `.agents/rustix-api.txt` in full into context.
- If you want '<pattern>', search the local symbol index with `rg`:
```sh
rg -i -- ::'<pattern>' .agents/rustix-api.txt
```
If a C/POSIX syscall name is not found (e.g. `kill` is `kill_process`), search `#[doc(alias` in the pinned crate source:
```sh
./.agent/rustix-alias.sh '<pattern>'
```

**Cargo Configuration**
- Scope Unix-only dependencies under `[target.'cfg(unix)'.dependencies]`.
- Enable only the minimal Cargo features required; each top-level module name maps 1-to-1 to its feature flag (e.g. `process::geteuid` needs `features = ["process"]`).

## Secure File Access Rules

* **temporary files:** Use `tempfile` crate for creating temporary files instead of rolling your own.

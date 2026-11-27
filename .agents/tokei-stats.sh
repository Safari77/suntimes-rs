#!/usr/bin/env bash
set -euo pipefail

DEFAULT_LANGUAGES=(
    Rust
    CSS
    FreeMarker
    HTML
    JavaScript
    Markdown
    Shell
    SQL
    SVG
    TOML
)

measure() {
    local lang=$1
    shift
    tokei -t "$lang" "$@" | awk '
        /^[[:space:]]*Total[[:space:]]/ { print $2, $3, $4, $5, $6; found = 1 }
        END { if (!found) print 0, 0, 0, 0, 0 }
    '
}

collect_stats() {
    local -n out_labels=$1
    local -n out_stats=$2
    shift 2

    out_labels=()
    out_stats=()

    local -A seen=()
    local lang stats
    for lang in "$@"; do
        if [[ -n "${seen[$lang]+x}" ]]; then
            continue
        fi
        seen[$lang]=1

        if [[ "$lang" == 'Rust' ]]; then
            stats='0 0 0 0 0'
            if [[ -e src ]]; then
                stats=$(measure Rust --exclude src/tests -- src)
            fi
            out_labels+=('Rust code')
            out_stats+=("$stats")

            stats='0 0 0 0 0'
            if [[ -e src/tests ]]; then
                stats=$(measure Rust -- src/tests)
            fi
            out_labels+=('Rust tests')
            out_stats+=("$stats")
        else
            stats=$(measure "$lang" -- .)
            out_labels+=("$lang")
            out_stats+=("$stats")
        fi
    done
}

print_diff() {
    local -n labels_ref=$1
    local -n before_ref=$2
    local -n after_ref=$3

    local rule='━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━'
    local t_f=0 t_l=0 t_c=0 t_m=0 t_b=0
    local i label b_f b_l b_c b_m b_b a_f a_l a_c a_m a_b
    local d_f d_l d_c d_m d_b
    local row
    local -a rows=()

    for i in "${!labels_ref[@]}"; do
        label=${labels_ref[i]}
        read -r b_f b_l b_c b_m b_b <<< "${before_ref[$label]}"
        read -r a_f a_l a_c a_m a_b <<< "${after_ref[i]}"

        d_f=$((a_f - b_f))
        d_l=$((a_l - b_l))
        d_c=$((a_c - b_c))
        d_m=$((a_m - b_m))
        d_b=$((a_b - b_b))

        if (( d_f == 0 && d_l == 0 && d_c == 0 && d_m == 0 && d_b == 0 )); then
            continue
        fi

        t_f=$((t_f + d_f))
        t_l=$((t_l + d_l))
        t_c=$((t_c + d_c))
        t_m=$((t_m + d_m))
        t_b=$((t_b + d_b))

        printf -v row ' %-19s %7d %12d %12d %12d %12d' \
            "$label" "$d_f" "$d_l" "$d_c" "$d_m" "$d_b"
        rows+=("$row")
    done

    if (( ${#rows[@]} == 0 )); then
        printf '\nAdded or removed code stats: no changes\n'
        return
    fi

    printf '\nAdded or removed code stats\n%s\n' "$rule"
    printf ' %-19s %7s %12s %12s %12s %12s\n' \
        'Language' 'Files' 'Lines' 'Code' 'Comments' 'Blanks'
    printf '%s\n' "$rule"
    printf '%s\n' "${rows[@]}"
    printf '%s\n' "$rule"
    printf ' %-19s %7d %12d %12d %12d %12d\n' \
        'Total' "$t_f" "$t_l" "$t_c" "$t_m" "$t_b"
    printf '%s\n' "$rule"
}

git_dir=$(git rev-parse --git-dir)
git_dir_real=$(realpath -- "$git_dir")
state_name=$(printf '%s' "$git_dir_real" | sha256sum)
state_name="${state_name%% *}"

state_dir="${TMP:-/tmp}"
state_file="${state_dir}/${state_name}"

cmd="${1:-}"
if (( $# > 0 )); then
    shift
fi
if [[ "${1:-}" == '--' ]]; then
    shift
fi

case "$cmd" in
    start)
        if (( $# > 0 )); then
            languages=("$@")
        else
            languages=("${DEFAULT_LANGUAGES[@]}")
        fi

        collect_stats before_labels before_stats "${languages[@]}"

        # Create 0600 temp file in target directory, flush, and atomically rename
        mkdir -p -m 0700 -- "$state_dir"
        tmp=$(mktemp -p "$state_dir" -- ".${state_name}.XXXXXX")
        trap 'rm -f -- "$tmp"' EXIT
        {
            for i in "${!before_labels[@]}"; do
                printf '%s\t%s\n' "${before_labels[i]}" "${before_stats[i]}"
            done
        } > "$tmp"
        sync -- "$tmp"
        mv -f -- "$tmp" "$state_file"
        trap - EXIT
        ;;
    diff)
        if [[ ! -f "$state_file" || -L "$state_file" || ! -r "$state_file" ]]; then
            printf 'Missing baseline file: %s\n' "$state_file" >&2
            exit 1
        fi

        declare -A before_map=()
        recorded_langs=()
        while IFS=$'\t' read -r label stats; do
            if [[ -z "$label" || -z "$stats" ]]; then
                printf 'Invalid baseline file: %s\n' "$state_file" >&2
                exit 1
            fi
            read -r f l c m b extra <<< "$stats"
            if [[ -z "$b" || -n "$extra" ]]; then
                printf 'Invalid baseline file: %s\n' "$state_file" >&2
                exit 1
            fi
            before_map[$label]=$stats
            if [[ "$label" == 'Rust code' || "$label" == 'Rust tests' ]]; then
                recorded_langs+=('Rust')
            else
                recorded_langs+=("$label")
            fi
        done < "$state_file"

        if (( ${#before_map[@]} == 0 )); then
            printf 'Invalid baseline file: %s\n' "$state_file" >&2
            exit 1
        fi

        if (( $# > 0 )); then
            languages=("$@")
        else
            languages=("${recorded_langs[@]}")
        fi

        collect_stats after_labels after_stats "${languages[@]}"

        for label in "${after_labels[@]}"; do
            if [[ -z "${before_map[$label]+x}" ]]; then
                printf 'Missing baseline entry for %s in %s\n' "$label" "$state_file" >&2
                exit 1
            fi
        done

        print_diff after_labels before_map after_stats
        ;;
    *)
        printf 'Usage: %s {start|diff} [language ...]\n' "$0" >&2
        exit 2
        ;;
esac

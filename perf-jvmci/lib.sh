#!/usr/bin/env bash
# Shared configuration and helpers for the JVMCI-removal compile-cost measurement.
# Sourced by the numbered scripts; not executable on its own.

set -euo pipefail

# ---------------------------------------------------------------------------------------------
# Configuration (override any of these in the environment)
# ---------------------------------------------------------------------------------------------

# Where the TornadoVM clone lives. Worktrees are created as siblings of it.
: "${TORNADO_REPO:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/..}"
TORNADO_REPO="$(cd "$TORNADO_REPO" && pwd)"

# Where built SDKs, probe classes and raw results are kept. Nothing here belongs in git.
: "${PERF_WORK:=$HOME/perf-jvmci-work}"

# Backend under measurement. One backend per campaign run: mixing them makes the table unreadable.
: "${BACKEND:=cuda}"

# Roots searched for JDKs, most specific first. /opt/jenkins/jdks is the CI layout (cyclone),
# ~/.sdkman/candidates/java the usual developer one.
: "${JDK_ROOTS:=/opt/jenkins/jdks $HOME/.sdkman/candidates/java}"

# How many fresh-JVM repetitions per configuration. Medians of >= 5; the campaign guidance calls
# for medians of 3+, and cold numbers are the noisiest thing here.
: "${REPS:=7}"

# Warm-same repetitions inside each JVM run.
: "${WARM_REPS:=20}"

# ---------------------------------------------------------------------------------------------
# The configuration matrix.
#
#   label | git ref | build JDK | make target
#
# Only two of these pairs are true A/Bs, and they answer different questions:
#
#   jvmci-jdk21 vs reflect-jdk21   MECHANISM. upstream/develop is a clean ancestor of the removal
#                                  branch (0 behind / 167 ahead), so the diff IS the JVMCI removal.
#                                  This is the number to defend.
#   jvmci-jdk25 vs reflect-jdk25   PRODUCT. upstream/jdk25 is the maintained JDK 25 line but has
#                                  diverged (730 behind develop, 437 of its own commits), so the
#                                  delta includes unrelated change. Quote it as "what a JDK 25 user
#                                  sees", never as the cost of the removal itself.
#
# reflect-jdk25 is not a separate build: the jdk22plus SDK is built once and RUN on 22..27, which
# is the property the removal branch exists to provide.
# ---------------------------------------------------------------------------------------------
CONFIGS=(
    "jvmci-jdk21|upstream/develop|21|jdk21"
    "reflect-jdk21|jdk27-jvmci-removal|21|jdk21"
    "jvmci-jdk25|upstream/jdk25|25|jdk25"
    "reflect-jdk22plus|jdk27-jvmci-removal|25|jdk22plus"
)

# Which JDKs each built SDK is RUN on. The baselines are pinned to their own JDK; the jdk22plus
# SDK fans out, which is simultaneously the portability check and the cross-JDK variation sweep.
run_jdks_for() {
    case "$1" in
        jvmci-jdk21|reflect-jdk21) echo "21" ;;
        jvmci-jdk25)               echo "25" ;;
        reflect-jdk22plus)         echo "${REFLECT_RUN_JDKS:-22 23 24 25 26 27}" ;;
        *) die "unknown configuration: $1" ;;
    esac
}

# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# Locate a JDK by feature version and VERIFY it, rather than trusting a path pattern: directory
# names lie (zulu27...jdk21 is a JDK 21), and picking the wrong one would silently invalidate a
# whole column of results.
resolve_jdk() {
    local feature="$1" cand home ver
    for root in $JDK_ROOTS; do
        [ -d "$root" ] || continue
        for cand in "$root"/*"$feature"*/Contents/Home "$root"/*/*"$feature"*/Contents/Home \
                    "$root"/*"$feature"* "$root"/*/*"$feature"*; do
            home="${cand%/}"
            [ -x "$home/bin/java" ] || continue
            ver=$("$home/bin/java" -version 2>&1 | head -1 | sed -n 's/.*version "\([0-9][0-9]*\).*/\1/p')
            if [ "$ver" = "$feature" ]; then echo "$home"; return 0; fi
        done
    done
    return 1
}

require_jdk() {
    local home
    home=$(resolve_jdk "$1") || die "no JDK $1 found under: $JDK_ROOTS"
    echo "$home"
}

sdk_dir_for()      { echo "$PERF_WORK/sdks/$1"; }
worktree_dir_for() { echo "$PERF_WORK/worktrees/$1"; }

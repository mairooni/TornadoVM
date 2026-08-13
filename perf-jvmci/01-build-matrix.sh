#!/usr/bin/env bash
#
# Build every (branch, JDK) SDK in the measurement matrix.
#
# Each configuration is built in its OWN git worktree, for two reasons that both bit us before:
#   * graalJars/ is per-checkout state. The removal branch DELETES the raw compiler jar after
#     relocating it, while develop/jdk25 still need that jar (they use --upgrade-module-path).
#     Sharing one checkout means whichever built second fails or, worse, links the wrong Graal.
#   * develop and the removal branch's jdk21 SDK carry the SAME Maven coordinate
#     (5.2.1-jdk21-dev) with different content. Separate trees plus copying each SDK out
#     immediately keeps them from being confused for one another.
#
# Builds run sequentially on purpose: they contend for ~/.m2 and for the machine, and a build
# racing a build is not a build you can attribute a later measurement to.
#
# Usage:  ./01-build-matrix.sh [label ...]      (default: all configurations)

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

mkdir -p "$PERF_WORK"/{sdks,worktrees,logs}

# The JDK 21 image the removal branch needs to source the frozen jdk.vm.ci.* classes from, on any
# build JDK. Harmless for the baseline branches, which ignore it.
export JVMCI_SOURCE_JDK="${JVMCI_SOURCE_JDK:-$(require_jdk 21)}"
log "JVMCI_SOURCE_JDK=$JVMCI_SOURCE_JDK"

wanted=("$@")
[ ${#wanted[@]} -eq 0 ] && wanted=("jvmci-jdk21" "reflect-jdk21" "jvmci-jdk25" "reflect-jdk22plus")

cd "$TORNADO_REPO"
log "fetching refs"
git fetch upstream develop jdk25 >/dev/null 2>&1 || log "WARNING: could not fetch upstream (offline?)"
git fetch origin jdk27-jvmci-removal >/dev/null 2>&1 || log "WARNING: could not fetch origin"

for cfg in "${CONFIGS[@]}"; do
    IFS='|' read -r label ref jdk target <<< "$cfg"
    skip=1
    for w in "${wanted[@]}"; do [ "$w" = "$label" ] && skip=0; done
    [ $skip -eq 1 ] && continue

    wt=$(worktree_dir_for "$label")
    sdk_out=$(sdk_dir_for "$label")
    java_home=$(require_jdk "$jdk")
    log "=== $label : ref=$ref buildJDK=$jdk target=$target backend=$BACKEND"
    log "    JAVA_HOME=$java_home"

    # Fresh worktree every time. A reused one keeps the previous run's graalJars and dist, which
    # is exactly the stale-state class of bug this campaign must not measure.
    if [ -d "$wt" ]; then
        git worktree remove --force "$wt" 2>/dev/null || rm -rf "$wt"
    fi
    git worktree add --detach "$wt" "$ref" >/dev/null 2>&1 || die "cannot create worktree for $ref"

    (
        cd "$wt"
        export JAVA_HOME="$java_home"
        export PATH="$JAVA_HOME/bin:$PATH"
        # Inherited from the environment if the caller set them (CMAKE_ROOT, CUDA_PATH, ...).
        log "    building (log: $PERF_WORK/logs/build-$label.log)"
        if ! make "$target" BACKEND="$BACKEND" > "$PERF_WORK/logs/build-$label.log" 2>&1; then
            tail -25 "$PERF_WORK/logs/build-$label.log" >&2
            die "build failed for $label -- see $PERF_WORK/logs/build-$label.log"
        fi
    )

    built=$(ls -d "$wt"/dist/*/*/ 2>/dev/null | head -1)
    [ -n "$built" ] || die "no SDK produced for $label"
    rm -rf "$sdk_out"
    mkdir -p "$(dirname "$sdk_out")"
    cp -r "${built%/}" "$sdk_out"
    log "    SDK -> $sdk_out"

    # Record what this SDK actually is, so a results table can never be traced to the wrong build.
    {
        echo "label=$label"
        echo "ref=$ref"
        echo "commit=$(git -C "$wt" rev-parse HEAD)"
        echo "build_jdk=$jdk"
        echo "build_java_home=$java_home"
        echo "make_target=$target"
        echo "backend=$BACKEND"
        echo "built_at=$(date -Is)"
    } > "$sdk_out/PERF-PROVENANCE"
done

log "done. SDKs:"
ls -1 "$PERF_WORK/sdks" >&2

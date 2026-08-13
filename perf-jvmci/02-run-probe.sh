#!/usr/bin/env bash
#
# Run the compile-cost probe against every built SDK, on every JDK that SDK supports.
#
# Two measurements per (SDK, JDK):
#   * profiler OFF, REPS fresh JVMs -- the headline wall-clock numbers (cold / warm-new / warm-same)
#   * profiler ON, 1 run            -- phase attribution only (Graal compile vs code gen vs driver)
#
# They are kept apart deliberately. TornadoVM's own profiler perturbs short task graphs, so a
# profiler-on run is a DIFFERENT workload: useful for saying *where* time goes, not *how much*.
#
# Usage:  ./02-run-probe.sh [label ...]      (default: every SDK present)

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

PROBE_SRC="$(dirname "${BASH_SOURCE[0]}")/probe/src/perfprobe"
RESULTS="$PERF_WORK/results"
mkdir -p "$RESULTS" "$PERF_WORK/logs"

# Derive the compile flags from the SDK itself rather than from a hardcoded table: the API jar's
# class-file version says which release to target, and minor 65535 marks preview (JDK 21 FFM).
# Guessing here is how the sample apps ended up compiling at release 21 with --enable-preview on
# JDK 22 and failing outright.
probe_javac_flags() {
    # Scan MANY classes, not one. javac stamps the preview flag only on the classes that actually
    # use a preview feature, so TaskGraph.class is minor=0 while FloatArray.class (FFM-backed) is
    # minor=65535 in the very same jar. Sniffing a single class picked the wrong answer and made
    # the probe fail to compile against the JDK 21 baseline.
    local api_jar="$1" tmp cls minor major max_major=0 preview=0
    tmp=$(mktemp -d)
    ( cd "$tmp" && unzip -qo "$api_jar" 'uk/ac/manchester/tornado/api/*.class' \
                                        'uk/ac/manchester/tornado/api/types/arrays/*.class' 2>/dev/null )
    while IFS= read -r cls; do
        minor=$(od -An -tu2 -j4 -N2 --endian=big "$cls" | tr -d ' ')
        major=$(od -An -tu2 -j6 -N2 --endian=big "$cls" | tr -d ' ')
        [ -n "$major" ] && [ "$major" -gt "$max_major" ] && max_major=$major
        [ "$minor" = "65535" ] && preview=1
    done < <(find "$tmp" -name '*.class' | head -200)
    rm -rf "$tmp"
    [ "$max_major" -gt 0 ] || die "cannot read any class-file version from $api_jar"
    local release=$((max_major - 44))
    if [ "$preview" = "1" ]; then echo "--release $release --enable-preview"; else echo "--release $release"; fi
}

labels=("$@")
if [ ${#labels[@]} -eq 0 ]; then
    mapfile -t labels < <(ls -1 "$PERF_WORK/sdks" 2>/dev/null)
fi
[ ${#labels[@]} -gt 0 ] || die "no SDKs in $PERF_WORK/sdks -- run 01-build-matrix.sh first"

for label in "${labels[@]}"; do
    sdk=$(sdk_dir_for "$label")
    [ -d "$sdk" ] || die "missing SDK: $sdk"
    api_jar=$(ls "$sdk"/share/java/tornado/tornado-api-*.jar | head -1)

    # Compile the probe against THIS SDK's own jars. Never against ~/.m2: develop and the removal
    # branch publish the same coordinate (5.2.1-jdk21-dev) with different content, so a repository
    # lookup would silently pick whichever was installed last.
    flags=$(probe_javac_flags "$api_jar")
    classes="$PERF_WORK/probe-classes/$label"
    rm -rf "$classes"; mkdir -p "$classes"

    for jdk in $(run_jdks_for "$label"); do
        java_home=$(resolve_jdk "$jdk") || { log "SKIP $label on JDK $jdk (not installed)"; continue; }

        # Compile once per SDK, using a JDK that can accept those flags.
        if [ ! -f "$classes/perfprobe/CompilePerf.class" ]; then
            log "compiling probe for $label ($flags)"
            # -g is mandatory, not cosmetic. The JVMCI code path reads the real LocalVariableTable
            # to name kernel parameters and NPEs without it (CUDANodeLIRBuilder.emitPrologue); the
            # reflection path synthesises one and does not care. Maven compiles with debug info by
            # default, so -g is also what a real application ships -- without it the baseline
            # cannot run at all and the A/B is impossible.
            "$java_home/bin/javac" $flags -g -nowarn -cp "$api_jar" -d "$classes" "$PROBE_SRC"/*.java \
                > "$PERF_WORK/logs/probe-compile-$label.log" 2>&1 \
                || { tail -15 "$PERF_WORK/logs/probe-compile-$label.log" >&2; die "probe compile failed for $label"; }
        fi

        out="$RESULTS/${label}__on-jdk${jdk}.tsv"
        : > "$out"

        # One discarded warm-up run. "Cold JVM" is not the same as "cold machine": the first run
        # after a build also faults the SDK jars and the GPU driver in from disk, which measured
        # ~4x higher warm-new times than every subsequent run. That is page cache, not compilation,
        # and it belongs in neither column.
        log "=== $label on JDK $jdk : warm-up (discarded) + $REPS fresh-JVM runs"
        JAVA_HOME="$java_home" TORNADOVM_HOME="$sdk" \
            "$sdk/bin/tornado" --jvm="-Dtornado.recover.bailout=False" \
            --classpath "$classes" perfprobe.CompilePerf 3 >/dev/null 2>&1 || true

        for rep in $(seq 1 "$REPS"); do
            # A fresh JVM per repetition is the whole point: "cold" means cold caches AND cold JIT,
            # and neither can be reset inside a running VM.
            if ! JAVA_HOME="$java_home" TORNADOVM_HOME="$sdk" \
                 "$sdk/bin/tornado" --jvm="-Dtornado.recover.bailout=False" \
                 --classpath "$classes" perfprobe.CompilePerf "$WARM_REPS" \
                 > "$PERF_WORK/logs/run-$label-jdk$jdk-$rep.out" 2>"$PERF_WORK/logs/run-$label-jdk$jdk-$rep.err"; then
                log "  rep $rep FAILED (see logs/run-$label-jdk$jdk-$rep.err)"
                tail -5 "$PERF_WORK/logs/run-$label-jdk$jdk-$rep.err" >&2
                continue
            fi
            awk -v rep="$rep" 'NR>1 && NF==4 {print $1"\t"$2"\t"$3"\t"$4"\t"rep}' \
                "$PERF_WORK/logs/run-$label-jdk$jdk-$rep.out" >> "$out"
            printf '.' >&2
        done
        echo >&2
        log "  -> $out ($(wc -l < "$out") samples)"

        # Compile-phase measurement, profiler on, REPS times.
        #
        # This is the PRIMARY signal for per-kernel compile cost, not the wall clock above. A
        # first-execution wall clock bundles compile with GPU execution, and the execution variance
        # (interquartile ranges of 3-5 ms) completely swamps a sub-millisecond compile delta -- every
        # tier came out statistically indistinguishable that way. The profiler's phase timers
        # exclude execution, which is exactly the noise that needed removing.
        #
        # One JSON per repetition so the aggregator can take medians and, crucially, measure the
        # spread of the CONTROL. TOTAL_DRIVER_COMPILE_TIME must not move between the two SDKs; how
        # much it moves anyway is this machine's noise floor, and no smaller delta is believable.
        for prep in $(seq 1 "$REPS"); do
            prof="$RESULTS/${label}__on-jdk${jdk}.profiler.$prep.json"
            rm -f "$prof"
            JAVA_HOME="$java_home" TORNADOVM_HOME="$sdk" \
                "$sdk/bin/tornado" --dumpProfiler "$prof" --jvm="-Dtornado.recover.bailout=False" \
                --classpath "$classes" perfprobe.CompilePerf 3 \
                > "$PERF_WORK/logs/prof-$label-jdk$jdk-$prep.out" 2>&1 || log "  (profiler rep $prep failed, non-fatal)"
        done
    done
done

log "results in $RESULTS"

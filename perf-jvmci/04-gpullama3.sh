#!/usr/bin/env bash
#
# GPULlama3 as the large real-world workload for the JVMCI-removal comparison.
#
# The synthetic probe prices ONE kernel at a time. An LLM is the opposite end of the scale: dozens
# of distinct kernels compiled at start-up, then a long steady-state decode loop. That splits the
# question in two, and the two halves are expected to move in OPPOSITE directions:
#
#   startup / time-to-first-token   pays every kernel compile, so this is where a per-kernel
#                                   compile regression accumulates -- and also where the removal
#                                   branch's smaller module path wins time back
#   tokens/second (steady state)    pure GPU decode, no compilation at all. This is the CONTROL:
#                                   it must NOT move, because nothing about metadata access can
#                                   reach a kernel that is already compiled and resident.
#
# Both SDKs run the SAME application jar, built once against the removal branch's API. That is
# deliberate: the removal branch made TornadoFunctions.TaskN Serializable, so a jar compiled
# against it carries writeReplace() and runs on both runtimes, whereas a jar compiled against
# develop's API lacks it and cannot run on the removal branch at all. One jar, two runtimes, is
# the only way the comparison isolates the runtime.
#
# Usage:  ./04-gpullama3.sh [reps]        (default 5)

set -euo pipefail

: "${PERF_WORK:=$HOME/perf-jvmci-work}"
: "${LLAMA_ROOT:=$PERF_WORK/gpullama3}"
: "${JDK21_HOME:?set JDK21_HOME to the JDK 21 both SDKs were built with}"

REPS="${1:-5}"
OUT="$PERF_WORK/gpullama3-results"
mkdir -p "$OUT"

# Fixed prompt, fixed seed, fixed temperature, fixed token count. Decode speed depends on how many
# tokens are produced and sampling is stochastic by default (the seed defaults to a timestamp), so
# without pinning all four, two runs are not comparable measurements of the same work.
PROMPT="Explain in detail how a modern GPU executes thousands of threads in parallel, and why that suits matrix multiplication."
SEED=42
TEMPERATURE=0.1
MAX_TOKENS=256

MODELS=(
    "llama3.2-1B-F16|$HOME/GPULlama3.java/Llama-3.2-1B-Instruct-F16.gguf"
    "llama3.2-1B-Q8|$HOME/GPULlama3.java/Llama-3.2-1B-Instruct-Q8_0.gguf"
    "qwen3-0.6B-F16|/opt/models/Qwen3-0.6B-f16.gguf"
    "qwen3-0.6B-Q8|/opt/models/Qwen3-0.6B-Q8_0.gguf"
)
SDKS=(jvmci-jdk21 reflect-jdk21)

{
    echo "captured_at=$(date -Is)"
    echo "host=$(hostname)"
    echo "gpu=$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader | head -1)"
    echo "jdk=$("$JDK21_HOME/bin/java" -version 2>&1 | head -1)"
    echo "llama_commit=$(git -C "$LLAMA_ROOT" rev-parse --short HEAD 2>/dev/null)"
    echo "prompt=$PROMPT"
    echo "seed=$SEED temperature=$TEMPERATURE max_tokens=$MAX_TOKENS reps=$REPS"
    for s in "${SDKS[@]}"; do
        echo "--- sdk $s"; sed 's/^/    /' "$PERF_WORK/sdks/$s/PERF-PROVENANCE" 2>/dev/null
    done
} > "$OUT/config.txt"
echo "config -> $OUT/config.txt"

export JAVA_HOME="$JDK21_HOME"
export LLAMA_ROOT
cd "$LLAMA_ROOT"

for m in "${MODELS[@]}"; do
    mname="${m%%|*}"; mpath="${m##*|}"
    [ -f "$mpath" ] || { echo "SKIP $mname (missing $mpath)"; continue; }
    for sdk in "${SDKS[@]}"; do
        tsv="$OUT/${mname}__${sdk}.tsv"
        : > "$tsv"
        echo "=== $mname on $sdk"
        # One discarded warm-up: the first touch of a multi-GB model file reads it from disk, and
        # that page-cache cost is not what we are measuring.
        TORNADOVM_HOME="$PERF_WORK/sdks/$sdk" timeout 900 ./llama-tornado --gpu \
            --model "$mpath" --prompt "$PROMPT" --seed "$SEED" --temperature "$TEMPERATURE" \
            -n "$MAX_TOKENS" > /dev/null 2>&1 || true
        for r in $(seq 1 "$REPS"); do
            log="$OUT/${mname}__${sdk}__run${r}.log"
            start=$(date +%s.%N)
            TORNADOVM_HOME="$PERF_WORK/sdks/$sdk" timeout 900 ./llama-tornado --gpu \
                --model "$mpath" --prompt "$PROMPT" --seed "$SEED" --temperature "$TEMPERATURE" \
                -n "$MAX_TOKENS" > "$log" 2>&1 || { echo "  rep $r FAILED"; continue; }
            wall=$(echo "$(date +%s.%N) - $start" | bc)
            # "achieved tok/s: 110.71. Tokens: 256, seconds: 2.31"
            toks=$(grep -oE "achieved tok/s: [0-9.]+" "$log" | tail -1 | grep -oE "[0-9.]+$")
            ntok=$(grep -oE "Tokens: [0-9]+" "$log" | tail -1 | grep -oE "[0-9]+$")
            secs=$(grep -oE "seconds: [0-9.]+" "$log" | tail -1 | grep -oE "[0-9.]+$")
            [ -n "$toks" ] || { echo "  rep $r: no tok/s in output"; continue; }
            # wall - decode = everything before steady state: JVM start, model load, and every
            # kernel compile. That is the half a compile regression would show up in.
            startup=$(echo "$wall - $secs" | bc)
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$r" "$toks" "$ntok" "$secs" "$wall" "$startup" >> "$tsv"
            printf '  rep %s: %s tok/s, decode %ss, wall %.2fs, startup %.2fs\n' "$r" "$toks" "$secs" "$wall" "$startup"
        done
    done
done
echo "results -> $OUT"

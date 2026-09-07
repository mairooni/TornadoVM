#!/bin/bash
# Build mma_async_stage2.ptx by splicing the new head (cp.async pipelined K-loop)
# with the unchanged BLOCK_9 output epilogue from the validated mma_async_stage1.ptx.
#
# Usage: ./splice_stage2.sh
# Run from this directory (alongside mma_async_stage2_head.ptx and your stage 1 file).
set -euo pipefail

HEAD="${HEAD:-mma_async_stage2_head.ptx}"
STAGE1="${STAGE1:-/home/mary/Projects/TornadoVM/mma_async_stage1.ptx}"
OUT="${OUT:-mma_async_stage2.ptx}"

if [ ! -f "$HEAD" ]; then
  echo "ERROR: head file not found: $HEAD" >&2; exit 1
fi
if [ ! -f "$STAGE1" ]; then
  echo "ERROR: stage 1 file not found: $STAGE1" >&2
  echo "Set STAGE1 env var to the correct path." >&2
  exit 1
fi

# Find line number where BLOCK_9 starts in stage 1
B9_LINE=$(grep -n "^BLOCK_9:" "$STAGE1" | head -1 | cut -d: -f1)
if [ -z "$B9_LINE" ]; then
  echo "ERROR: could not find BLOCK_9: in $STAGE1" >&2; exit 1
fi
echo "BLOCK_9 starts at line $B9_LINE in stage 1"

# Concatenate: head + BLOCK_9 onward from stage 1
cat "$HEAD" > "$OUT"
tail -n +"$B9_LINE" "$STAGE1" >> "$OUT"

echo "Wrote $OUT ($(wc -l < "$OUT") lines)"
echo "Tail of output (should end with 'ret;' and '}'):"
tail -3 "$OUT"

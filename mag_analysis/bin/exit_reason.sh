#!/bin/bash
# Describe a tool's exit status for a warning: exit_reason.sh <code>
# Codes above 128 are a signal (128 + N): the tool was stopped from outside,
# which says nothing about the input.
case "$1" in
    137) echo "killed (SIGKILL: out of memory, or the job ended)" ;;
    143) echo "terminated (SIGTERM: a time limit, a cancelled job or a stopped run)" ;;
    *)   if [ "$1" -gt 128 ] 2>/dev/null; then echo "stopped by signal $(($1 - 128))"
         else echo "failed; see its output above"; fi ;;
esac

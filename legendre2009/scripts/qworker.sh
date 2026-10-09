#!/bin/bash
# qworker.sh SLOTS -- shared job queue for the campaign (queue/*.job on the
# common shared filesystem).  Starts jobs while fewer than SLOTS of this
# user's mpiruns run on this host; a job is claimed by an atomic mkdir.
# When the queue is empty and this host is idle the worker exits; the last
# worker to exit writes ALLDONE into logs/campaign_*.log.
P=$(cd "$(dirname "$0")/.." && pwd); Q=$P/queue; LOG=$P/logs/queue.log
SLOTS=$1; H=$(hostname -s)
while :; do
  n=$(pgrep -u "$USER" -x mpirun | wc -l)
  if [ "$n" -lt "$SLOTS" ]; then
    job=""
    for f in $(cd "$Q" && ls *.job 2>/dev/null | sort); do
      if mkdir "$Q/${f%.job}.claim" 2>/dev/null; then job=$f; break; fi
    done
    if [ -n "$job" ]; then
      echo "$(date +%H:%M:%S) $H start ${job%.job}" >> "$LOG"
      ( bash -c "$(cat "$Q/$job")"; echo "$(date +%H:%M:%S) $H end   ${job%.job}" >> "$LOG" ) &
      sleep 10; continue
    fi
    [ "$n" -eq 0 ] && break
  fi
  sleep 30
done
echo "$(date +%H:%M:%S) $H worker exits" >> "$LOG"
if [ "$(grep -c 'worker exits' "$LOG")" -ge 2 ]; then
  for l in "$P"/logs/campaign_*.log; do echo ALLDONE >> "$l"; done
fi

#!/bin/bash
# ============================================================
# Start ParaView server on the cluster for remote rendering.
# Use with local ParaView client via SSH tunnel.
#
# Setup (do once per session):
#   1. Reconnect SSH with port forwarding:
#        ssh -L 11111:localhost:11111 cluster
#
#   2. On cluster, in your work directory:
#        util/paraview/pvserver_start.sh
#
#   3. On local machine, open ParaView:
#        File → Connect → Add Server
#          Name: cluster
#          Server Type: Client / Server
#          Host: localhost
#          Port: 11111
#        → Connect
#
#   4. In connected ParaView: File → Open → browse cluster files
# ============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- ParaView environment ---
if [ -f "$HOME/usr/paraview511.sh" ]; then
  source "$HOME/usr/paraview511.sh"
else
  echo "ERROR: ParaView setup script not found at $HOME/usr/paraview511.sh"
  exit 1
fi

# --- Anti-flicker for Mesa ---
export vblank_mode=0

PORT="${1:-11111}"

echo "========================================="
echo " ParaView Server (pvserver)"
echo " Port: $PORT"
echo "========================================="
echo ""
echo "On your LOCAL machine, open ParaView and:"
echo "  File → Connect → Add Server"
echo "    Name: cluster"
echo "    Host: localhost"
echo "    Port: $PORT"
echo "  → Connect"
echo ""
echo "IMPORTANT: SSH must have port forwarding:"
echo "  ssh -L ${PORT}:localhost:${PORT} cluster"
echo ""

env LD_LIBRARY_PATH="$PARAVIEW_HOME/lib:/opt/ohpc/pub/compiler/gcc/8.3.0/lib64:$LD_LIBRARY_PATH" \
    pvserver --mesa --server-port="$PORT"

#!/bin/bash
# ============================================================
# Open a JOREK VTK file with ParaView 5.11.1
#
# Usage:
#   util/open_paraview.sh [vtk_file]
#
#   Default: opens jorek_tmp.vtk in current directory
#
# Examples:
#   util/open_paraview.sh                          # open ./jorek_tmp.vtk
#   util/open_paraview.sh jorek_tmp_step10.vtk     # open specific file
#   util/open_paraview.sh ../case2/jorek_tmp.vtk   # open from another dir
# ============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOREK_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# --- ParaView environment ---
if [ -f "$HOME/usr/paraview511.sh" ]; then
  source "$HOME/usr/paraview511.sh"
else
  echo "ERROR: ParaView setup script not found at $HOME/usr/paraview511.sh"
  exit 1
fi

# --- Check DISPLAY ---
if [ -z "$DISPLAY" ]; then
  echo "ERROR: DISPLAY not set. Reconnect with: ssh -X cluster"
  exit 1
fi

# --- Anti-flicker: disable vsync for Mesa software rendering ---
export vblank_mode=0

# --- VTK file ---
VTK_FILE="${1:-jorek_tmp.vtk}"

if [ ! -f "$VTK_FILE" ]; then
  echo "ERROR: VTK file not found: $VTK_FILE"
  echo ""
  echo "To generate a VTK file from a JOREK restart file, run:"
  echo "  util/convert2vtk.sh jorek_restart_XXXX.h5"
  exit 1
fi

# --- One-time defaults setup (color map, etc.) ---
PARAVIEW_SETUP="$JOREK_ROOT/util/paraview/setup_defaults.py"
if [ -f "$PARAVIEW_SETUP" ]; then
  python3 "$PARAVIEW_SETUP" 2>/dev/null || true
fi

echo "Opening: $VTK_FILE"
echo "ParaView: $(paraview --version 2>/dev/null || echo '5.11.1')"
echo ""

# Set ParaView's own libs locally (first in LD_LIBRARY_PATH) to avoid
# polluting the global environment for Intel MPI-dependent programs.
export PARAVIEW_DATA_FILE="$VTK_FILE"
PARAVIEW_STARTUP="$JOREK_ROOT/util/paraview/startup.py"

env LD_LIBRARY_PATH="$PARAVIEW_HOME/lib:/opt/ohpc/pub/compiler/gcc/8.3.0/lib64:$LD_LIBRARY_PATH" \
    paraview --mesa --backend=swr --script="$PARAVIEW_STARTUP" &

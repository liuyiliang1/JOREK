"""
ParaView startup script — load data and configure per-session view properties.
Called by util/open_paraview.sh via paraview --script=.

Persistent defaults (color map etc.) are handled by setup_defaults.py.
"""
import os
from paraview.simple import *


def configure_views():
    """Apply per-session view settings to all render views."""
    views = GetRenderViews()
    for view in views:
        # ---- Disable LightKit ----
        try:
            view.UseLight = 0
        except Exception:
            pass

        # ---- Hide orientation axes ----
        try:
            view.OrientationAxesVisibility = 0
        except Exception:
            pass


def main():
    vtk_file = os.environ.get("PARAVIEW_DATA_FILE", "")

    if vtk_file and os.path.isfile(vtk_file):
        reader = OpenDataFile(vtk_file)
        Show(reader)

    configure_views()
    Render()


if __name__ == "__main__":
    main()

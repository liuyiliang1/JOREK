#!/usr/bin/env python3
"""
One-time ParaView settings setup — sets default color map to Rainbow Desaturated.
Run once; settings persist across ParaView sessions.
"""
import json
import os
import sys


SETTINGS_FILE = os.path.expanduser("~/.config/ParaView/ParaView-UserSettings.json")
TARGET_DEFAULT = "Rainbow Desaturated"


def setup_default_color_map():
    if not os.path.exists(SETTINGS_FILE):
        print(f"Settings file not found: {SETTINGS_FILE}")
        print("Run ParaView once first to generate it, then re-run this script.")
        return False

    with open(SETTINGS_FILE, "r") as f:
        settings = json.load(f)

    # Parse the TransferFunctionPresets groups
    groups_key = "TransferFunctionPresets"
    if groups_key not in settings:
        print(f"No '{groups_key}' in settings, cannot configure.")
        return False

    groups = json.loads(settings[groups_key]["Groups"])

    modified = False
    for group in groups:
        if group.get("groupName") == "Default":
            presets = group["presets"]
            if TARGET_DEFAULT in presets and presets[0] != TARGET_DEFAULT:
                presets.remove(TARGET_DEFAULT)
                presets.insert(0, TARGET_DEFAULT)
                modified = True
                print(f"Moved '{TARGET_DEFAULT}' to front of Default group.")
            elif presets[0] == TARGET_DEFAULT:
                print(f"'{TARGET_DEFAULT}' is already the default — nothing to do.")
                return True
            else:
                print(f"WARNING: '{TARGET_DEFAULT}' not found in Default group!")
                return False
            break

    if modified:
        settings[groups_key]["Groups"] = json.dumps(groups)
        # Backup original
        backup = SETTINGS_FILE + ".bak"
        if not os.path.exists(backup):
            os.rename(SETTINGS_FILE, backup)
            print(f"Backup saved to {backup}")
        with open(SETTINGS_FILE, "w") as f:
            json.dump(settings, f, indent="\t")
        print("Settings updated. Restart ParaView for changes to take effect.")

    return True


if __name__ == "__main__":
    success = setup_default_color_map()
    sys.exit(0 if success else 1)

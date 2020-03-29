"""Management of pygimli user configuration file."""

import json
import os
import sys


# pyGIMLi default configuration
rc = {
    'lang': 'eng', # 'eng', 'german', 'de', 'ger'
    'unitStyle': 2,
    # 1: german style: 'value in unit'
    # 2: default style: 'value (unit)'
    'view3D': 'auto',
    # auto: Use pyvista if installed or set it to 'fallback' to force fallback mode
    'globalCache': True
}


def getConfigPath():
    r"""Return the user space path for configuration files.

    - Windows: 'C:\Documents and Settings\username\Appdata\pygimli'
    - Linu: '~/.config/pygimli' (if not overwritten by $XDG_CONFIG_HOME)
    - Mac: '~/Library/Preferences/pygimli'

    See Also
    --------
    pygimli.getCachePath
    """
    appname = "pygimli"
    system = sys.platform

    if system == "win32":
        path = os.path.join(os.environ['APPDATA'])
    elif system == "darwin":
        path = os.path.expanduser('~/Library/Preferences/')
    else:
        path = os.getenv('XDG_CONFIG_HOME', os.path.expanduser("~/.config"))

    return os.path.join(path, appname)


configpath = getConfigPath()
configfile = os.path.join(configpath, "config.json")

if os.path.isfile(configfile):
    # print("Loading user configuration file at " + configfile)
    with open(configfile) as cfg:
        userrc = json.load(cfg)

    # Check if user config holds all keys and update if necessary
    if len(userrc.keys()) != len(rc.keys()):
        for key in rc:
            if key not in userrc:
                print("Updating user configuration with", key, "=", rc[key])
                userrc[key] = rc[key]
        with open(configfile, "w") as cfg:
            json.dump(userrc, cfg, indent=4, sort_keys=True)
    rc = userrc
else:
    print("Creating default user configuration file at " + configfile)
    with open(configfile, "w") as cfg:
        json.dump(rc, cfg, indent=4, sort_keys=True)

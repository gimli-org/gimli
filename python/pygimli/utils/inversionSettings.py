# -*- coding: utf-8 -*-
"""TODO documentme!.

SORT ME.
"""
from pygimli.io import opt_import

json = opt_import("json", "read and write inversion settings")


class InversionSettings(dict):
    """Extends the built-in dict with methods for saving and loading a file."""

    def __init__(self, *args, **kwargs):
        """Constructor.

        Initialize the settings object by either using a file or sepcifying the
        settings directly. Works exatly like a dict(), only with extended
        capabilities like saving and loading to disk.
        """
        try:
            self.filename = kwargs.pop('filename')
        except KeyError:
            print('Creating new settings object.')
        else:
            print('Loading settings from: "{}".'.format(self.filename))
        finally:
            super(InversionSettings, self).__init__(*args, **kwargs)

        if hasattr(self, 'filename'):
            self.update(InversionSettings.load(self.filename))

    def save(self, f):
        """Save settings.

        Saves the settings object to disk.

        Parameters
        ----------
        f : str
            filename or file handle to use.
        """
        with open(f, 'w') as thefile:
            json.dump(self, thefile)

    @classmethod
    def load(cls, f):
        """Load settings from.

        Load a saved settings object from disk.
        Parameters
        ----------
        f : str
            filename or file handle to use.
        """
        with open(f, 'r') as thefile:
            return json.load(thefile)

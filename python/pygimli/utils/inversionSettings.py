import cPickle as pickle


class InversionSettings(dict):
    """
    Extends the built-in dict with methods for saving and loading a file.
    """
    def __init__(self, *args, **kwargs):
        """
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
            self.update(InversionSettings.loadtxt(self.filename))

    def savetxt(self, f):
        """
        Saves the settings object to disk.

        ::input
        f : the filename or file handle to use
        """

        with open(f, 'w') as thefile:
            pickle.dump(self, thefile)

    @classmethod
    def loadtxt(cls, f):
        """
        Class method that loads a saved settings object from disk.

        ::input
        f : the filename or file handle to use
        """

        with open(f, 'r') as thefile:
            return pickle.load(thefile)

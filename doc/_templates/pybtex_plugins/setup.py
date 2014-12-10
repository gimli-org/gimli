from setuptools import setup

setup(name='plugins',
      version='0.1.0',
      entry_points={'pybtex.style.formatting': ['mystyle = plugins:MyStyle'],
                    'pybtex.style.labels': ['alpha = plugins:Alpha']},
      py_modules=['plugins'])
      

#!/ussr/bin/env python
# -*- coding: utf-8 -*-
'''
    
TUTORIAL NAME
-------------

*This introductory sentence should state the intent and goal of the tutorial. Keep it brief.*

*This next block should state any assumptions that you the writer are making. Present them in list form.*

Cite something :cite:`Zienkiewicz1977` or [Zienkiewicz1977]_

glossary ref :term:`numpy`

GIMLi api ref :gimliapi:`GIMLI::Cell`

pygimli ref ??????

'''

import pygimli as g

print((g.versionStr()))

"""
Last output

.. lastcout::

"""

"""
Invoking :term:`Matplotlib`.

Please see http://matplotlib.org/1.3.0/faq/usage_faq.html#general-concepts
"""

import matplotlib.pyplot as plt
plt.plot(1, 1, 'x')


plt.show()

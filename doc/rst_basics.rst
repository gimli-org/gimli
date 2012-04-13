Section How to document with reStructuredText
=============================================

.. only:: latex

    This is the latex-only version for how to document :term:`GIMLi` 
    
This is the common version for how to document :term:`GIMLi`. 
Can we uses cites [RueckerGueSpi2006]_, [GuentherRueSpi2006]_, [RueckerGuen2011]_.

reStructuredText is a plaintext markup syntax that is interpreted by sphinx to generate this html or pdf document

* http://docutils.sourceforge.net/rst.html
*`bibstuff.sphinxext.bibref`

this may be `a link to sphinx`_

one asterisk: *text* for emphasis (italics),

two asterisks: **text** for strong emphasis (boldface), and

backquotes: ``text`` for code samples.

.. _a link to sphinx: http://sphinx.pocoo.org/


SubSection 
----------

SubSubSection
.............



.. this is a comment
    comment comment comment

...please read :ref:`INSTALLATION`, :ref:`USAGE` and :ref:`FAQ`, as well

special box::

    this looks like a listing

::

    this looks like a listing

.. note::

   to indicate something that may be of interest

or a

.. warning::

   to indicate something that could cause serious problems.


manually reference:

.. [GuentherRueSpi2006] T. Günther, C. Rücker and K. Spitzer. Three-dimensional modelling and inversion of dc resistivity data incorporating topography -- I. Modelling. Geophys. J. Int. 166, 495--505 

.. [RueckerGueSpi2006] C. Rücker, T. Günther and K. Spitzer. Three-dimensional modelling and inversion of dc resistivity data incorporating topography -- I. Modelling. Geophys. J. Int. 166, 495--505 


.. warning::
    automatic reference by bibtex using `bibstuff.sphinxext.bibref` does not support utf8 so far

.. bibmissing:: libgimli.bib
    :sort:





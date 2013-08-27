Section reStructuredText Cockbock
=================================

.. this is a comment
    comment comment comment

.. only:: latex

    This is the latex-only version for how to document :term:`GIMLi` 
    

SubSection 
----------

This is the common version for how to document :term:`GIMLi`. 
Can we uses cites [RueckerGueSpi2006]_, [GuentherRueSpi2006]_, [RueckerGuen2011]_.

Check for autocite :cite:`RueckerGuen2011`, :cite:`GuentherRueSpi2006`, :cite:`RueckerGueSpi2006`


    * Terms, e.g., :term:`GIMLi` are defined in glossery.rst

    * References, e.g., :term:`Matplotlib` are defined in glossery.rst

    * Modules, e.g., :mod:`Matplotlib` are defined in glossery.rst


SubSubSection
.............

automatic show api documentation

.. autofunction:: pygimli.randN

reStructuredText is a plaintext markup syntax that is interpreted by sphinx to generate this html or pdf document

* http://docutils.sourceforge.net/rst.html
* bibstuff.sphinxext.bibref

this may be `a link to sphinx`_

one asterisk: *text* for emphasis (italics),

two asterisks: **text** for strong emphasis (boldface), and

backquotes: ``text`` for code samples.

.. _a link to sphinx: http://sphinx.pocoo.org/


...please read :ref:`sec:install`, and :ref:`sec:faq`, as well

special box::

    this looks like a listing

::

    this looks like a listing

Admonitions
...........

(Most themes style only “note” and “warning” specially.)

.. note::

   to indicate something that may be of interest

or a

.. warning::

   to indicate something that could cause serious problems.

.. admonition:: Whatever we need for captionname
    
    simple box with personal annotation
    also possible boxes are: attention, caution, danger, error, hint, important, note, tip, warning and the generic admonition. 


.. note::

    ::

        pip install sphinxcontrib-bibtex
    
    install http://sphinxcontrib-bibtex.readthedocs.org/en/0.2.8/quickstart.html#installation






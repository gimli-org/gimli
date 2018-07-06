{{ fullname | escape | underline }}

.. currentmodule:: {{ fullname }}

.. automodule:: {{ fullname }}

Overview
--------

{% block functions %}
{% if functions %}

.. rubric:: Functions

.. autosummary::

{% for function in functions %}
  {{ function }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block classes %}
{% if classes %}

.. rubric:: Classes

.. autosummary::

{% for class in classes %}
  {{ class }}
{%- endfor %}
{% endif %}
{% endblock %}


{% block exceptions %}
{% if exceptions %}

Exceptions
----------

.. autosummary::

{% for exception in exceptions %}
  {{ exception }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block functionapi %}
{% if functions %}

Functions
---------

{% for function in functions %}
.. rubric:: {{ function }}

.. autofunction:: {{ function }}

.. include:: {{fullname}}.{{function}}.examples
    :end-line: 3

.. include:: {{fullname}}.{{function}}.examples
    :start-line: 4

.. raw:: html

     <div style='clear:both'></div>
{%- endfor %}
{% endif %}
{% endblock %}

{% block classesapi %}
{% if classes %}

Classes
-------

{% for class in classes %}
.. rubric:: {{ class }}
.. autoclass:: {{ class }}
  :members:
  :undoc-members:
  :special-members: __init__
  :inherited-members:

{%- endfor %}
{% endif %}
{% endblock %}

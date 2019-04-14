:mod:`{{module}}`.{{objname}}
{{ underline }}===============

.. meta::
   :robots: noindex

.. warning::
   **DEPRECATED**


.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__
   .. automethod:: __call__
   {% endblock %}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div class="clearer"></div>

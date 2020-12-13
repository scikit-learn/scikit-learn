:mod:`{{module}}`.{{objname}}
{{ underline }}===============

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __call__
   {% endblock %}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div class="clearer"></div>

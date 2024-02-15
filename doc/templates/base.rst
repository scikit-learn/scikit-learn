{{ objname | escape | underline(line="=") }}

{% if objtype == "module" %}
.. automodule:: {{ fullname }}
{% else %}
.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}
{% endif %}

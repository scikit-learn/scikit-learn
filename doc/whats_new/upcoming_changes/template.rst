{% if render_title %}
{% if versiondata.name %}
{{ versiondata.name }} {{ versiondata.version }} ({{ versiondata.date }})
{{ top_underline * ((versiondata.name + versiondata.version + versiondata.date)|length + 4)}}
{% else %}
{{ versiondata.version }} ({{ versiondata.date }})
{{ top_underline * ((versiondata.version + versiondata.date)|length + 3)}}
{% endif %}
{% endif %}

{% set underline = underlines[0] %}
{% for section, content_per_category in sections.items() %}
{% if content_per_category and section %}
{{ section }}
{{ underline * section|length }}
{% endif %}
{% for category, content in content_per_category.items() %}
{% for text, issue_links in content.items() %}
{% set tag = definitions[category]['name'] %}
{# TODO a bit hacky replace first space with tag like |Fix|, is there a cleaner way? #}
{% set text_with_tag = text if category == 'other' else text.replace(' ', ' ' + tag + ' ', 1) %}
{# TODO I don't really understand how you can have multiple issue_links but this is a list so need to join #}
{{ text_with_tag }} {{ issue_links|join(', ') }}

{% endfor %}
{% endfor %}
{% endfor %}

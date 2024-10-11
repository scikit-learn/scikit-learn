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
{% for section, content_per_category in sections.items() if content_per_category and section %}
{{ section }}
{{ underline * section|length }}
{# We loop over definitions because, contrary to content_per_category, it follow the category order as defined in pyproject.toml #}
{% for category in definitions if category in content_per_category %}
{% set content = content_per_category[category] %}
{% for text, issue_links in content.items() %}
{% set tag = definitions[category]['name'] %}
{# TODO a bit hacky replace first space with tag like |Fix|, is there a cleaner way? #}
{% set text_with_tag = text if category == 'other' else text.replace(' ', ' ' + tag + ' ', 1) %}
{# issue_links is a list so need to join. For our purposes, issue_links is always of length 1 #}
{{ text_with_tag }} {{ issue_links|join(', ') }}

{% endfor %}
{% endfor %}
{% endfor %}

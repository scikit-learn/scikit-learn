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
{# content_per_category does not respect the category order in definitions so need reordering ... #}
{% set ordered_content_per_category = dict() %}
{% for cat in definitions %}
{% if cat in content_per_category %}
{% set _ = ordered_content_per_category.update({cat: content_per_category[cat]}) %}
{% endif %}
{% endfor %}
{% for category, content in ordered_content_per_category.items() %}
{% for text, issue_links in content.items() %}
{% set tag = definitions[category]['name'] %}
{# TODO a bit hacky replace first space with tag like |Fix|, is there a cleaner way? #}
{% set text_with_tag = text if category == 'other' else text.replace(' ', ' ' + tag + ' ', 1) %}
{# TODO I don't really understand how you can have multiple issue_links but this is a list so need to join #}
{{ text_with_tag }} {{ issue_links|join(', ') }}

{% endfor %}
{% endfor %}
{% endfor %}

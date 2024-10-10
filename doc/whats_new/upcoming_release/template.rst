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
{% if content_per_category %}
{{ section }}
{{ underline * section|length }}
{% endif %}

{% for category, content in content_per_category.items() %}
{% for text, issue_links in content.items() %}
{# TODO for later we can generate tags like |Fix|, |API|, etc ... but need to special case top-level section
|{{ definitions[category]['name'] }}| {{ text }} {{ issue_links|join(', ') }}
#}
{# TODO I don't really understand how you can have multiple issue_links but this is a list so need join
#}
{{ text }} {{ issue_links|join(', ') }}

{% endfor %}
{% endfor %}
{% endfor %}

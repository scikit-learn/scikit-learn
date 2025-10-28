{# Displays a link to the .rst source of the current page. #}
{% if show_source and has_source and sourcename %}
  <div class="tocsection sourcelink">
    <a href="{{ pathto('_sources/' + sourcename, true)|e }}">
      <i class="fa-solid fa-file-lines"></i> {{ _("Show Source") }}
    </a>
  </div>
{% endif %}

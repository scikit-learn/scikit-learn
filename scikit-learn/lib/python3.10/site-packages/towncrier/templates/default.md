{#-
══════════════════════════════════════════════════════════════════════════════
TOWNCRIER MARKDOWN TEMPLATE
══════════════════════════════════════════════════════════════════════════════

─── Macro: heading ─────────────────────────────────────────────────────────
Purpose:
   Generates Markdown headings with the appropriate number of # characters.
   Based on header_prefix (default: "#") and the level argument.

Arguments:
   level             The relative heading level (1=#, 2=##, 3=###, etc.)
-#}
{%- macro heading(level) -%}
    {{- "#" * ( header_prefix | length + level -1 ) }}
{%- endmacro -%}

{#-
─── Variable: newline ──────────────────────────────────────────────────────
Purpose:
    Consistent newline handling. -#}
{%- set newline = "\n" -%}

{#- ════════════════════════ TEMPLATE GENERATION ════════════════════════ -#}
{#- ─── TITLE HEADING ─── #}
{#- render_title is false when title_format is specified in the config #}
{%- if render_title %}
    {%- if versiondata.name %}
        {{- heading(1) ~ " " ~ versiondata.name ~ " " ~ versiondata.version ~ " (" ~ versiondata.date ~ ")" ~ newline }}
    {%- else %}
        {{- heading(1) ~ " " ~ versiondata.version ~ " (" ~ versiondata.date ~ ")" ~ newline }}
    {%- endif %}
{%- endif %}
{#- If title_format is specified, we start with a new line #}
{{- newline }}

{%- for section, _ in sections.items() %}
    {#- ─── SECTION HEADING ─── #}
    {%- if section %}
        {{- newline }}
        {{- heading(2) ~ " " ~ section ~ newline }}
        {{- newline }}
    {%- endif %}

    {%- if sections[section] %}

        {%- for category, val in definitions.items() if category in sections[section] %}
            {%- set issue_pks = [] %}
            {#- ─── CATEGORY HEADING ─── #}
            {#- Increase heading level if section is not present #}
            {{- heading(3 if section else 2) ~" " ~ definitions[category]['name'] ~ newline }}
            {{- newline }}

            {#- ─── RENDER ENTRIES ─── #}
            {%- for text, values in sections[section][category].items() %}
                {#- Prepare the string of issue numbers (e.g., "#1, #9, #142") #}
                {%- set issue_pks = [] %}
                {%- for v_issue in values %}
                    {%- set _= issue_pks.append(v_issue.split(": ", 1)[0]) %}
                {%- endfor %}
                {%- set issues_list = issue_pks | join(", ") %}

                {#- Check if text contains a sublist #}
                {%- set text_has_sublist = (("\n  - " in text) or ("\n  * " in text)) %}

                {#- CASE 1: No text, only issues #}
                {#- Output: -  #1, #9, #142 #}
                {%- if not text and issues_list %}
                    {{- "- " ~ issues_list ~ newline }}

                {#- Cases where both text and issues exist #}
                {%- elif text and issues_list %}
                    {%- if text_has_sublist %}
                        {#- CASE 3: Text with sublist #}
                        {#- Output: - TEXT\n\n  (#1, #9, #142) #}
                        {{- "- " ~ text ~ newline ~ newline ~ "  (" ~ issues_list ~ ")" ~ newline }}
                    {%- else %}
                        {#- CASE 2: Text, no sublist #}
                        {#- Output: - TEXT (#1, #9, #142) #}
                        {{- "- " ~ text ~ " (" ~ issues_list ~ ")" ~ newline }}
                    {%- endif %}

                {%- elif text %}
                    {#- Implicit Case: Text, but no issues #}
                    {#- Output: - TEXT #}
                    {{- "- " ~ text ~ newline }}
                {%- endif %}
            {%- endfor %}

            {#- New line between list and link references #}
            {{- newline }}

            {#- Link references #}
            {%- if issues_by_category[section][category] and "]: " in issues_by_category[section][category][0] %}
                {%- for issue in issues_by_category[section][category] %}
                    {{- issue ~ newline }}
                {%- endfor %}
                {{- newline }}
            {%- endif %}

            {#- No changes in this category #}
            {%- if sections[section][category]|length == 0 %}
                {{- newline }}
                {{- "No significant changes." ~ newline * 2 }}
            {%- endif %}
        {%- endfor %}
    {%- else %}
        {#- No changes in this section #}
        {{- "No significant changes." ~ newline * 2 }}
    {%- endif %}
{%- endfor %}
{#-
Newline at the end of the rendered newsfile content.
In this way the there are 2 newlines between the latest release and the previous release content.
-#}
{{- newline -}}

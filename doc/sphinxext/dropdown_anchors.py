import re

from docutils import nodes
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx_design.dropdown import dropdown_main, dropdown_title


class DropdownAnchorAdder(SphinxPostTransform):
    """Insert anchor links to the sphinx-design dropdowns.

    Some of the dropdowns were originally headers that had automatic anchors, so we
    need to make sure that the old anchors still work. See the original implementation
    (in JS): https://github.com/scikit-learn/scikit-learn/pull/27409

    The structure of each sphinx-design dropdown node is expected to be:

    <dropdown_main ...>
        <dropdown_title ...>
            ...icon      <-- This exists if the "icon" option of the sphinx-design
                             dropdown is set; we do not use it in our documentation

            ...title     <-- This may contain multiple nodes, e.g. literal nodes if
                             there are inline codes; we use the concatenated text of
                             all these nodes to generate the anchor ID

            Here we insert the anchor link!

            <container ...>  <-- The "dropdown closed" marker
            <container ...>  <-- The "dropdown open" marker
        </dropdown_title>
        <container...>
            ...main contents
        </container>
    </dropdown_main>
    """

    default_priority = 9999  # Apply later than everything else
    formats = ["html"]

    def run(self):
        """Run the post transformation."""
        # Counter to store the duplicated summary text to add it as a suffix in the
        # anchor ID
        anchor_id_counters = {}

        for sd_dropdown in self.document.findall(dropdown_main):
            # Grab the dropdown title
            sd_dropdown_title = sd_dropdown.next_node(dropdown_title)

            # Concatenate the text of relevant nodes as the title text
            # Since we do not have the prefix icon, the relevant nodes are the very
            # first child node until the third last node (last two are markers)
            title_text = "".join(
                node.astext() for node in sd_dropdown_title.children[:-2]
            )

            # The ID uses the first line, lowercased, with spaces replaced by dashes;
            # suffix the anchor ID with a counter if it already exists
            anchor_id = re.sub(r"\s+", "-", title_text.strip().split("\n")[0]).lower()
            if anchor_id in anchor_id_counters:
                anchor_id_counters[anchor_id] += 1
                anchor_id = f"{anchor_id}-{anchor_id_counters[anchor_id]}"
            else:
                anchor_id_counters[anchor_id] = 1
            sd_dropdown["ids"].append(anchor_id)

            # Create the anchor element and insert after the title text; we do this
            # directly with raw HTML
            anchor_html = (
                f'<a class="headerlink" href="#{anchor_id}" '
                'title="Link to this dropdown">#</a>'
            )
            anchor_node = nodes.raw("", anchor_html, format="html")
            sd_dropdown_title.insert(-2, anchor_node)  # before the two markers


def setup(app):
    app.add_post_transform(DropdownAnchorAdder)

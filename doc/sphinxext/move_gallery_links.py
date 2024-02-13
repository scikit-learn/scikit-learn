from pathlib import Path

from bs4 import BeautifulSoup
from sphinx.util.display import status_iterator
from sphinx.util.logging import getLogger

logger = getLogger(__name__)


def move_gallery_links(app, exception):
    if exception is not None:
        return

    for gallery_dir in app.config.sphinx_gallery_conf["gallery_dirs"]:
        html_gallery_dir = Path(app.builder.outdir, gallery_dir)

        # Get all gallery example files to be tweaked; tuples (file, docname)
        flat = []
        for file in html_gallery_dir.rglob("*.html"):
            if file.name in ("index.html", "sg_execution_times.html"):
                # These are not gallery example pages, skip
                continue

            # Extract the documentation name from the path
            docname = file.relative_to(app.builder.outdir).with_suffix("").as_posix()
            if docname in app.config.html_context["redirects"]:
                # This is a redirected page, skip
                continue
            if docname not in app.project.docnames:
                # This should not happen, warn
                logger.warning(f"Document {docname} not found but {file} exists")
                continue
            flat.append((file, docname))

        for html_file, _ in status_iterator(
            flat,
            length=len(flat),
            summary="Tweaking gallery links... ",
            verbosity=app.verbosity,
            stringify_func=lambda x: x[1],  # display docname
        ):
            with html_file.open("r", encoding="utf-8") as f:
                html = f.read()
            soup = BeautifulSoup(html, "html.parser")

            # Find the secondary sidebar; it should exist in all gallery example pages
            secondary_sidebar = soup.find("div", class_="sidebar-secondary-items")
            if secondary_sidebar is None:
                logger.warning(f"Secondary sidebar not found in {html_file}")
                continue

            def _create_secondary_sidebar_component(items):
                """Create a new component in the secondary sidebar.

                `items` should be a list of dictionaries with "element" being the bs4
                tag of the component and "title" being the title (None if not needed).
                """
                component = soup.new_tag("div", **{"class": "sidebar-secondary-item"})
                for item in items:
                    item_wrapper = soup.new_tag("div")
                    item_wrapper.append(item["element"])
                    if item["title"]:
                        item_wrapper["title"] = item["title"]
                    component.append(item_wrapper)
                secondary_sidebar.append(component)

            def _create_download_link(link, is_jupyter=False):
                """Create a download link to be appended to a component.

                `link` should be the bs4 tag of the original download link, either for
                the Python source code (is_jupyter=False) of for the Jupyter notebook
                (is_jupyter=True). `link` will not be removed; instead the whole
                footnote would be removed where `link` is located.

                This returns a dictionary with "element" being the bs4 tag of the new
                download link and "title" being the name of the file to download.
                """
                new_link = soup.new_tag("a", href=link["href"], download="")

                # Place a download icon at the beginning of the new link
                download_icon = soup.new_tag("i", **{"class": "fa-solid fa-download"})
                new_link.append(download_icon)

                # Create the text of the new link; it is shortend to fit better into
                # the secondary sidebar
                link_type = "Jupyter notebook" if is_jupyter else "source code"
                new_text = soup.new_string(f"Download {link_type}")
                new_link.append(new_text)

                # Get the file name to download and use it as the title of the new link
                # which will show up when hovering over the link; the file name is
                # expected to be in the last span of `link`
                link_spans = link.find_all("span")
                title = link_spans[-1].text if link_spans else None

                return {"element": new_link, "title": title}

            def _create_badge_link(link):
                """Create a badge link to be appended to a component.

                `link` should be the bs4 tag of the original badge link, either for
                binder or JupyterLite. `link` will not be removed; instead the whole
                footnote would be removed where `link` is located.

                This returns a dictionary with "element" being the bs4 tag of the new
                download link and "title" being `None` (no need).
                """
                new_link = soup.new_tag("a", href=link["href"])

                # The link would essentially be an anchor wrapper outside the image of
                # the badge; we get the src and alt attributes by finding the original
                # image and limit the height to 20px (fixed) so that the secondary
                # sidebar will appear neater
                badge_img = link.find("img")
                new_img = soup.new_tag(
                    "img", src=badge_img["src"], alt=badge_img["alt"], height=20
                )
                new_link.append(new_img)

                return {"element": new_link, "title": None}

            try:
                # `sg_note` is the "go to the end" note at the top of the page
                # `sg_footer` is the footer with the download links and badge links
                # These will be removed at the end if new links are successfully created
                sg_note = soup.find("div", class_="sphx-glr-download-link-note")
                sg_footer = soup.find("div", class_="sphx-glr-footer")

                # Move the download links into the secondary sidebar
                py_link = sg_footer.find("div", class_="sphx-glr-download-python").a
                ipy_link = sg_footer.find("div", class_="sphx-glr-download-jupyter").a
                _create_secondary_sidebar_component(
                    [
                        _create_download_link(py_link, is_jupyter=False),
                        _create_download_link(ipy_link, is_jupyter=True),
                    ]
                )

                # Move the badge links into the secondary sidebar
                lite_link = sg_footer.find("div", class_="lite-badge").a
                binder_link = sg_footer.find("div", class_="binder-badge").a
                _create_secondary_sidebar_component(
                    [_create_badge_link(lite_link), _create_badge_link(binder_link)]
                )

                # Remove the sourcelink component from the secondary sidebar; the reason
                # we do not remove it by configuration is that we need the secondary
                # sidebar to be present for this script to work, while in-page toc alone
                # could have been empty
                sourcelink = secondary_sidebar.find("div", class_="sourcelink")
                if sourcelink is not None:
                    sourcelink.parent.extract()  # because sourcelink has a wrapper div

                # Remove the the top note and the whole footer
                sg_note.extract()
                sg_footer.extract()

            except Exception as e:
                # If any step fails we directly skip the file
                logger.warning(f"Failed to tweak gallery links in {html_file}: {e}")
                continue

            # Write the modified file back
            with html_file.open("w", encoding="utf-8") as f:
                f.write(str(soup))


def setup(app):
    # Default priority is 500 which sphinx-gallery uses for its build-finished events;
    # we need a larger (i.e. lower) priority to run after sphinx-gallery
    app.connect("build-finished", move_gallery_links, priority=900)

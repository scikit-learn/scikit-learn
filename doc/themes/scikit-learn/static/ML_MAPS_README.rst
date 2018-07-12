Machine Learning Cheat Sheet (for scikit-learn)
===============================================

This document is intended to explain how to edit
the machine learning cheat sheet, originally created
by Andreas Mueller:

(http://peekaboo-vision.blogspot.de/2013/01/machine-learning-cheat-sheet-for-scikit.html)

The image is made interactive using an imagemap, and uses the jQuery Map Hilight plugin module
by David Lynch (http://davidlynch.org/projects/maphilight/docs/) to highlight
the different items on the image upon mouseover.

Modifying the map on the docs is currently a little bit tedious,
so I'll try to make it as simple as possible.

1. Editing the layout of the map and its paths.
------------------------------------------------

Use a Graphics editor like Inkscape Vector Graphics Editor
to open the ml_map.svg file, in this folder. From there
you can move objects around, etc. as you need.

Save when done, and make sure to export a .PNG file
to replace the old-outdated ml_map.png, as that file
is used as a background image.

2. Accessing the paths of the SVG file and exporting them.
----------------------------------------------------------

Use an image manipulation package like GIMP Image Editor to open
the ml_map.svg file, in this folder. With GIMP, make sure
to select 'Import paths'.

Once the image has been opened, you can see all imported paths on the paths tab.
You can edit these or create new paths. In GIMP, right-clicking one of the
paths and choosing: Path Tool will allow you to see the paths on
the image. The paths will be exported later and will be used to
make the click able regions on our image map.

3. Export paths as SVG files
----------------------------

After you've edited a path or created a new one, right click it on
the paths menu and choose 'Export Path..'. This way we extract just
that path on its own as 'new_area.svg' for example.

4. Edit the SVG file
---------------------
Using a script made by David Lynch, we will convert the svg files into
html maps. To do this, open the svg file in question in any text editor.
Make sure that the 'width' and 'height' are not in 'in' or 'px', i.e
"100" is OK, but "100px" or "1.25in" are not.

Then wrap the <path> tags in <g> and </g> tags.
Then the file is ready for the script.

5. From SVG to HTML map
-----------------------

Use the provided svg2imagemap.py script on your edited svg file:

$ python svg2imagemap.py new_area.svg

where new_area.svg is our file.

6. Add the new map to the main html file
------------------------------------------

Copy the code from the newly created 'new_area.html'
file. Open the ml_map.html file.

Add the <area href=....... ></area> that you copied
after the last </area> tag in the ml_map.html file.

Add the link address to 'href' and a tooltip to
'title' within your <area ...> tag.

If you wish to add the green and blue hover effect
to the area, add
data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'

to your  area tag, as done in the other <area..> tags above.

Save the file, and you're done.

-----------------------------------------------------

I'll take some time to make some scripts to automate this process
a bit more at some point, as it is not difficult to do,
but tedious.

-Jaques Grobler

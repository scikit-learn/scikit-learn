#!/usr/local/bin/python

"""
This script converts a subset of SVG into an HTML imagemap

Note *subset*.  It only handles <path> elements, for which it only pays
attention to the M and L commands.  Futher, it only notices the "translate"
transform.

It was written to generate the examples in the documentation for maphilight,
and thus is very squarely aimed at handling several SVG maps from wikipedia.
It *assumes* that all the <path>s it will need are inside a <g>.  Any <path>
outside of a <g> will be ignored.

It takes several possible arguments, in the form:
$ svn2imagemap.py FILENAME [x y [group1 group2 ... groupN]]

FILENAME must be the name of an SVG file.  All other arguments are optional.

x and y, if present, are the dimensions of the image you'll be creating from
the SVG.  If not present, it assumes the values of the width and height
attributes in the SVG file.

group1 through groupN are group ids.  If only want particular groups used,
enter their ids here and all others will be ignored.
"""

import os
import re
import sys
import xml.dom.minidom

import parse_path

if len(sys.argv) == 1:
    sys.exit("svn2imagemap.py FILENAME [x y [group1 group2 ... groupN]]")
if not os.path.exists(sys.argv[1]):
    sys.exit("Input file does not exist")
x, y, groups = None, None, None
if len(sys.argv) >= 3:
    x = float(sys.argv[2])
    y = float(sys.argv[3])
    if len(sys.argv) > 3:
        groups = sys.argv[4:]

svg_file = xml.dom.minidom.parse(sys.argv[1])
svg = svg_file.getElementsByTagName('svg')[0]

raw_width = float(svg.getAttribute('width'))
raw_height = float(svg.getAttribute('height'))
width_ratio = x and (x / raw_width) or 1
height_ratio = y and (y / raw_height) or 1

if groups:
    elements = [g for g in svg.getElementsByTagName('g') if (g.hasAttribute('id') and g.getAttribute('id') in groups)]
    elements.extend([p for p in svg.getElementsByTagName('path') if (p.hasAttribute('id') and p.getAttribute('id') in groups)])
else:
    elements = svg.getElementsByTagName('g')

parsed_groups = {}
for e in elements:
    paths = []
    if e.nodeName == 'g':
        for path in e.getElementsByTagName('path'):
            points = parse_path.get_points(path.getAttribute('d'))
            for pointset in points:
                paths.append([path.getAttribute('id'), pointset])
    else:
        points = parse_path.get_points(e.getAttribute('d'))
        for pointset in points:
            paths.append([e.getAttribute('id'), pointset])
    if e.hasAttribute('transform'):
        print e.getAttribute('id'), e.getAttribute('transform')
        for transform in re.findall(r'(\w+)\((-?\d+.?\d*),(-?\d+.?\d*)\)', e.getAttribute('transform')):
            if transform[0] == 'translate':
                x_shift = float(transform[1])
                y_shift = float(transform[2])
                for path in paths:
                    path[1] = [(p[0] + x_shift, p[1] + y_shift) for p in path[1]]

    parsed_groups[e.getAttribute('id')] = paths

out = []
for g in parsed_groups:
    for path in parsed_groups[g]:
        out.append('<area href="#" title="%s" shape="poly" coords="%s"></area>' %
            (path[0], ', '.join([("%d,%d" % (p[0]*width_ratio, p[1]*height_ratio)) for p in path[1]])))

outfile = open(sys.argv[1].replace('.svg', '.html'), 'w')
outfile.write('\n'.join(out))

# tifffile.py

# Copyright (c) 2008-2025, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

r"""Read and write TIFF files.

Tifffile is a Python library to

(1) store NumPy arrays in TIFF (Tagged Image File Format) files, and
(2) read image and metadata from TIFF-like files used in bioimaging.

Image and metadata can be read from TIFF, BigTIFF, OME-TIFF, GeoTIFF,
Adobe DNG, ZIF (Zoomable Image File Format), MetaMorph STK, Zeiss LSM,
ImageJ hyperstack, Micro-Manager MMStack and NDTiff, SGI, NIHImage,
Olympus FluoView and SIS, ScanImage, Molecular Dynamics GEL,
Aperio SVS, Leica SCN, Roche BIF, PerkinElmer QPTIFF (QPI, PKI),
Hamamatsu NDPI, Argos AVS, and Philips DP formatted files.

Image data can be read as NumPy arrays or Zarr 2 arrays/groups from strips,
tiles, pages (IFDs), SubIFDs, higher-order series, and pyramidal levels.

Image data can be written to TIFF, BigTIFF, OME-TIFF, and ImageJ hyperstack
compatible files in multi-page, volumetric, pyramidal, memory-mappable,
tiled, predicted, or compressed form.

Many compression and predictor schemes are supported via the imagecodecs
library, including LZW, PackBits, Deflate, PIXTIFF, LZMA, LERC, Zstd,
JPEG (8 and 12-bit, lossless), JPEG 2000, JPEG XR, JPEG XL, WebP, PNG, EER,
Jetraw, 24-bit floating-point, and horizontal differencing.

Tifffile can also be used to inspect TIFF structures, read image data from
multi-dimensional file sequences, write fsspec ReferenceFileSystem for
TIFF files and image file sequences, patch TIFF tag values, and parse
many proprietary metadata formats.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD-3-Clause
:Version: 2025.5.10
:DOI: `10.5281/zenodo.6795860 <https://doi.org/10.5281/zenodo.6795860>`_

Quickstart
----------

Install the tifffile package and all dependencies from the
`Python Package Index <https://pypi.org/project/tifffile/>`_::

    python -m pip install -U tifffile[all]

Tifffile is also available in other package repositories such as Anaconda,
Debian, and MSYS2.

The tifffile library is type annotated and documented via docstrings::

    python -c "import tifffile; help(tifffile)"

Tifffile can be used as a console script to inspect and preview TIFF files::

    python -m tifffile --help

See `Examples`_ for using the programming interface.

Source code and support are available on
`GitHub <https://github.com/cgohlke/tifffile>`_.

Support is also provided on the
`image.sc <https://forum.image.sc/tag/tifffile>`_ forum.

Requirements
------------

This revision was tested with the following requirements and dependencies
(other versions may work):

- `CPython <https://www.python.org>`_ 3.10.11, 3.11.9, 3.12.10, 3.13.3 64-bit
- `NumPy <https://pypi.org/project/numpy/>`_ 2.2.5
- `Imagecodecs <https://pypi.org/project/imagecodecs/>`_ 2025.3.30
  (required for encoding or decoding LZW, JPEG, etc. compressed segments)
- `Matplotlib <https://pypi.org/project/matplotlib/>`_ 3.10.3
  (required for plotting)
- `Lxml <https://pypi.org/project/lxml/>`_ 5.3.1
  (required only for validating and printing XML)
- `Zarr <https://pypi.org/project/zarr/>`_ 2.18.7
  (required only for opening Zarr stores; Zarr 3 is not compatible)
- `Fsspec <https://pypi.org/project/fsspec/>`_ 2025.3.2
  (required only for opening ReferenceFileSystem files)

Revisions
---------

2025.5.10

- Pass 5111 tests.
- Raise ValueError when using zarr 3 (#296).
- Fall back to compression.zstd on Python >= 3.14 if no imagecodecs.
- Remove doctest command line option.
- Support Python 3.14.

2025.3.30

- Fix for imagecodecs 2025.3.30.

2025.3.13

- Change bytes2str to decode only up to first NULL character (breaking).
- Remove stripnull function calls to reduce overhead (#285).
- Deprecate stripnull function.

2025.2.18

- Fix julian_datetime milliseconds (#283).
- Remove deprecated dtype arguments from imread and FileSequence (breaking).
- Remove deprecated imsave and TiffWriter.save function/method (breaking).
- Remove deprecated option to pass multiple values to compression (breaking).
- Remove deprecated option to pass unit to resolution (breaking).
- Remove deprecated enums from TIFF namespace (breaking).
- Remove deprecated lazyattr and squeeze_axes functions (breaking).

2025.1.10

- Improve type hints.
- Deprecate Python 3.10.

2024.12.12

- Read PlaneProperty from STK UIC1Tag (#280).
- Allow 'None' as alias for COMPRESSION.NONE and PREDICTOR.NONE (#274).
- Zarr 3 is not supported (#272).

2024.9.20

- Fix writing colormap to ImageJ files (breaking).
- Improve typing.
- Drop support for Python 3.9.

2024.8.30

- Support writing OME Dataset and some StructuredAnnotations elements.

2024.8.28

- Fix LSM scan types and dimension orders (#269, breaking).
- Use IO[bytes] instead of BinaryIO for typing (#268).

2024.8.24

- Do not remove trailing length-1 dimension writing non-shaped file (breaking).
- Fix writing OME-TIFF with certain modulo axes orders.
- Make imshow NaN aware.

2024.8.10

- Relax bitspersample check for JPEG, JPEG2K, and JPEGXL compression (#265).

2024.7.24

- Fix reading contiguous multi-page series via Zarr store (#67).

2024.7.21

- Fix integer overflow in product function caused by numpy types.
- Allow tag reader functions to fail.

2024.7.2

- Enable memmap to create empty files with non-native byte order.
- Deprecate Python 3.9, support Python 3.13.

2024.6.18

- Ensure TiffPage.nodata is castable to dtype (breaking, #260).
- Support Argos AVS slides.

2024.5.22

- Derive TiffPages, TiffPageSeries, FileSequence, StoredShape from Sequence.
- Truncate circular IFD chain, do not raise TiffFileError (breaking).
- Deprecate access to TiffPages.pages and FileSequence.files.
- Enable DeprecationWarning for enums in TIFF namespace.
- Remove some deprecated code (breaking).
- Add iccprofile property to TiffPage and parameter to TiffWriter.write.
- Do not detect VSI as SIS format.
- Limit length of logged exception messages.
- Fix docstring examples not correctly rendered on GitHub (#254, #255).

- …

Refer to the CHANGES file for older revisions.

Notes
-----

TIFF, the Tagged Image File Format, was created by the Aldus Corporation and
Adobe Systems Incorporated.

Tifffile supports a subset of the TIFF6 specification, mainly 8, 16, 32, and
64-bit integer, 16, 32, and 64-bit float, grayscale and multi-sample images.
Specifically, CCITT and OJPEG compression, chroma subsampling without JPEG
compression, color space transformations, samples with differing types, or
IPTC, ICC, and XMP metadata are not implemented.

Besides classic TIFF, tifffile supports several TIFF-like formats that do not
strictly adhere to the TIFF6 specification. Some formats allow file and data
sizes to exceed the 4 GB limit of the classic TIFF:

- **BigTIFF** is identified by version number 43 and uses different file
  header, IFD, and tag structures with 64-bit offsets. The format also adds
  64-bit data types. Tifffile can read and write BigTIFF files.
- **ImageJ hyperstacks** store all image data, which may exceed 4 GB,
  contiguously after the first IFD. Files > 4 GB contain one IFD only.
  The size and shape of the up to 6-dimensional image data can be determined
  from the ImageDescription tag of the first IFD, which is Latin-1 encoded.
  Tifffile can read and write ImageJ hyperstacks.
- **OME-TIFF** files store up to 8-dimensional image data in one or multiple
  TIFF or BigTIFF files. The UTF-8 encoded OME-XML metadata found in the
  ImageDescription tag of the first IFD defines the position of TIFF IFDs in
  the high-dimensional image data. Tifffile can read OME-TIFF files (except
  multi-file pyramidal) and write NumPy arrays to single-file OME-TIFF.
- **Micro-Manager NDTiff** stores multi-dimensional image data in one
  or more classic TIFF files. Metadata contained in a separate NDTiff.index
  binary file defines the position of the TIFF IFDs in the image array.
  Each TIFF file also contains metadata in a non-TIFF binary structure at
  offset 8. Downsampled image data of pyramidal datasets are stored in
  separate folders. Tifffile can read NDTiff files. Version 0 and 1 series,
  tiling, stitching, and multi-resolution pyramids are not supported.
- **Micro-Manager MMStack** stores 6-dimensional image data in one or more
  classic TIFF files. Metadata contained in non-TIFF binary structures and
  JSON strings define the image stack dimensions and the position of the image
  frame data in the file and the image stack. The TIFF structures and metadata
  are often corrupted or wrong. Tifffile can read MMStack files.
- **Carl Zeiss LSM** files store all IFDs below 4 GB and wrap around 32-bit
  StripOffsets pointing to image data above 4 GB. The StripOffsets of each
  series and position require separate unwrapping. The StripByteCounts tag
  contains the number of bytes for the uncompressed data. Tifffile can read
  LSM files of any size.
- **MetaMorph Stack, STK** files contain additional image planes stored
  contiguously after the image data of the first page. The total number of
  planes is equal to the count of the UIC2tag. Tifffile can read STK files.
- **ZIF**, the Zoomable Image File format, is a subspecification of BigTIFF
  with SGI's ImageDepth extension and additional compression schemes.
  Only little-endian, tiled, interleaved, 8-bit per sample images with
  JPEG, PNG, JPEG XR, and JPEG 2000 compression are allowed. Tifffile can
  read and write ZIF files.
- **Hamamatsu NDPI** files use some 64-bit offsets in the file header, IFD,
  and tag structures. Single, LONG typed tag values can exceed 32-bit.
  The high bytes of 64-bit tag values and offsets are stored after IFD
  structures. Tifffile can read NDPI files > 4 GB.
  JPEG compressed segments with dimensions >65530 or missing restart markers
  cannot be decoded with common JPEG libraries. Tifffile works around this
  limitation by separately decoding the MCUs between restart markers, which
  performs poorly. BitsPerSample, SamplesPerPixel, and
  PhotometricInterpretation tags may contain wrong values, which can be
  corrected using the value of tag 65441.
- **Philips TIFF** slides store padded ImageWidth and ImageLength tag values
  for tiled pages. The values can be corrected using the DICOM_PIXEL_SPACING
  attributes of the XML formatted description of the first page. Tile offsets
  and byte counts may be 0. Tifffile can read Philips slides.
- **Ventana/Roche BIF** slides store tiles and metadata in a BigTIFF container.
  Tiles may overlap and require stitching based on the TileJointInfo elements
  in the XMP tag. Volumetric scans are stored using the ImageDepth extension.
  Tifffile can read BIF and decode individual tiles but does not perform
  stitching.
- **ScanImage** optionally allows corrupted non-BigTIFF files > 2 GB.
  The values of StripOffsets and StripByteCounts can be recovered using the
  constant differences of the offsets of IFD and tag values throughout the
  file. Tifffile can read such files if the image data are stored contiguously
  in each page.
- **GeoTIFF sparse** files allow strip or tile offsets and byte counts to be 0.
  Such segments are implicitly set to 0 or the NODATA value on reading.
  Tifffile can read GeoTIFF sparse files.
- **Tifffile shaped** files store the array shape and user-provided metadata
  of multi-dimensional image series in JSON format in the ImageDescription tag
  of the first page of the series. The format allows multiple series,
  SubIFDs, sparse segments with zero offset and byte count, and truncated
  series, where only the first page of a series is present, and the image data
  are stored contiguously. No other software besides Tifffile supports the
  truncated format.

Other libraries for reading, writing, inspecting, or manipulating scientific
TIFF files from Python are
`aicsimageio <https://pypi.org/project/aicsimageio>`_,
`apeer-ometiff-library
<https://github.com/apeer-micro/apeer-ometiff-library>`_,
`bigtiff <https://pypi.org/project/bigtiff>`_,
`fabio.TiffIO <https://github.com/silx-kit/fabio>`_,
`GDAL <https://github.com/OSGeo/gdal/>`_,
`imread <https://github.com/luispedro/imread>`_,
`large_image <https://github.com/girder/large_image>`_,
`openslide-python <https://github.com/openslide/openslide-python>`_,
`opentile <https://github.com/imi-bigpicture/opentile>`_,
`pylibtiff <https://github.com/pearu/pylibtiff>`_,
`pylsm <https://launchpad.net/pylsm>`_,
`pymimage <https://github.com/ardoi/pymimage>`_,
`python-bioformats <https://github.com/CellProfiler/python-bioformats>`_,
`pytiff <https://github.com/FZJ-INM1-BDA/pytiff>`_,
`scanimagetiffreader-python
<https://gitlab.com/vidriotech/scanimagetiffreader-python>`_,
`SimpleITK <https://github.com/SimpleITK/SimpleITK>`_,
`slideio <https://gitlab.com/bioslide/slideio>`_,
`tiffslide <https://github.com/bayer-science-for-a-better-life/tiffslide>`_,
`tifftools <https://github.com/DigitalSlideArchive/tifftools>`_,
`tyf <https://github.com/Moustikitos/tyf>`_,
`xtiff <https://github.com/BodenmillerGroup/xtiff>`_, and
`ndtiff <https://github.com/micro-manager/NDTiffStorage>`_.

References
----------

- TIFF 6.0 Specification and Supplements. Adobe Systems Incorporated.
  https://www.adobe.io/open/standards/TIFF.html
  https://download.osgeo.org/libtiff/doc/
- TIFF File Format FAQ. https://www.awaresystems.be/imaging/tiff/faq.html
- The BigTIFF File Format.
  https://www.awaresystems.be/imaging/tiff/bigtiff.html
- MetaMorph Stack (STK) Image File Format.
  http://mdc.custhelp.com/app/answers/detail/a_id/18862
- Image File Format Description LSM 5/7 Release 6.0 (ZEN 2010).
  Carl Zeiss MicroImaging GmbH. BioSciences. May 10, 2011
- The OME-TIFF format.
  https://docs.openmicroscopy.org/ome-model/latest/
- UltraQuant(r) Version 6.0 for Windows Start-Up Guide.
  http://www.ultralum.com/images%20ultralum/pdf/UQStart%20Up%20Guide.pdf
- Micro-Manager File Formats.
  https://micro-manager.org/wiki/Micro-Manager_File_Formats
- ScanImage BigTiff Specification.
  https://docs.scanimage.org/Appendix/ScanImage+BigTiff+Specification.html
- ZIF, the Zoomable Image File format. https://zif.photo/
- GeoTIFF File Format https://gdal.org/drivers/raster/gtiff.html
- Cloud optimized GeoTIFF.
  https://github.com/cogeotiff/cog-spec/blob/master/spec.md
- Tags for TIFF and Related Specifications. Digital Preservation.
  https://www.loc.gov/preservation/digital/formats/content/tiff_tags.shtml
- CIPA DC-008-2016: Exchangeable image file format for digital still cameras:
  Exif Version 2.31.
  http://www.cipa.jp/std/documents/e/DC-008-Translation-2016-E.pdf
- The EER (Electron Event Representation) file format.
  https://github.com/fei-company/EerReaderLib
- Digital Negative (DNG) Specification. Version 1.7.1.0, September 2023.
  https://helpx.adobe.com/content/dam/help/en/photoshop/pdf/DNG_Spec_1_7_1_0.pdf
- Roche Digital Pathology. BIF image file format for digital pathology.
  https://diagnostics.roche.com/content/dam/diagnostics/Blueprint/en/pdf/rmd/Roche-Digital-Pathology-BIF-Whitepaper.pdf
- Astro-TIFF specification. https://astro-tiff.sourceforge.io/
- Aperio Technologies, Inc. Digital Slides and Third-Party Data Interchange.
  Aperio_Digital_Slides_and_Third-party_data_interchange.pdf
- PerkinElmer image format.
  https://downloads.openmicroscopy.org/images/Vectra-QPTIFF/perkinelmer/PKI_Image%20Format.docx
- NDTiffStorage. https://github.com/micro-manager/NDTiffStorage
- Argos AVS File Format.
  https://github.com/user-attachments/files/15580286/ARGOS.AVS.File.Format.pdf

Examples
--------

Write a NumPy array to a single-page RGB TIFF file:

>>> import numpy
>>> data = numpy.random.randint(0, 255, (256, 256, 3), 'uint8')
>>> imwrite('temp.tif', data, photometric='rgb')

Read the image from the TIFF file as NumPy array:

>>> image = imread('temp.tif')
>>> image.shape
(256, 256, 3)

Use the `photometric` and `planarconfig` arguments to write a 3x3x3 NumPy
array to an interleaved RGB, a planar RGB, or a 3-page grayscale TIFF:

>>> data = numpy.random.randint(0, 255, (3, 3, 3), 'uint8')
>>> imwrite('temp.tif', data, photometric='rgb')
>>> imwrite('temp.tif', data, photometric='rgb', planarconfig='separate')
>>> imwrite('temp.tif', data, photometric='minisblack')

Use the `extrasamples` argument to specify how extra components are
interpreted, for example, for an RGBA image with unassociated alpha channel:

>>> data = numpy.random.randint(0, 255, (256, 256, 4), 'uint8')
>>> imwrite('temp.tif', data, photometric='rgb', extrasamples=['unassalpha'])

Write a 3-dimensional NumPy array to a multi-page, 16-bit grayscale TIFF file:

>>> data = numpy.random.randint(0, 2**12, (64, 301, 219), 'uint16')
>>> imwrite('temp.tif', data, photometric='minisblack')

Read the whole image stack from the multi-page TIFF file as NumPy array:

>>> image_stack = imread('temp.tif')
>>> image_stack.shape
(64, 301, 219)
>>> image_stack.dtype
dtype('uint16')

Read the image from the first page in the TIFF file as NumPy array:

>>> image = imread('temp.tif', key=0)
>>> image.shape
(301, 219)

Read images from a selected range of pages:

>>> images = imread('temp.tif', key=range(4, 40, 2))
>>> images.shape
(18, 301, 219)

Iterate over all pages in the TIFF file and successively read images:

>>> with TiffFile('temp.tif') as tif:
...     for page in tif.pages:
...         image = page.asarray()
...

Get information about the image stack in the TIFF file without reading
any image data:

>>> tif = TiffFile('temp.tif')
>>> len(tif.pages)  # number of pages in the file
64
>>> page = tif.pages[0]  # get shape and dtype of image in first page
>>> page.shape
(301, 219)
>>> page.dtype
dtype('uint16')
>>> page.axes
'YX'
>>> series = tif.series[0]  # get shape and dtype of first image series
>>> series.shape
(64, 301, 219)
>>> series.dtype
dtype('uint16')
>>> series.axes
'QYX'
>>> tif.close()

Inspect the "XResolution" tag from the first page in the TIFF file:

>>> with TiffFile('temp.tif') as tif:
...     tag = tif.pages[0].tags['XResolution']
...
>>> tag.value
(1, 1)
>>> tag.name
'XResolution'
>>> tag.code
282
>>> tag.count
1
>>> tag.dtype
<DATATYPE.RATIONAL: 5>

Iterate over all tags in the TIFF file:

>>> with TiffFile('temp.tif') as tif:
...     for page in tif.pages:
...         for tag in page.tags:
...             tag_name, tag_value = tag.name, tag.value
...

Overwrite the value of an existing tag, for example, XResolution:

>>> with TiffFile('temp.tif', mode='r+') as tif:
...     _ = tif.pages[0].tags['XResolution'].overwrite((96000, 1000))
...

Write a 5-dimensional floating-point array using BigTIFF format, separate
color components, tiling, Zlib compression level 8, horizontal differencing
predictor, and additional metadata:

>>> data = numpy.random.rand(2, 5, 3, 301, 219).astype('float32')
>>> imwrite(
...     'temp.tif',
...     data,
...     bigtiff=True,
...     photometric='rgb',
...     planarconfig='separate',
...     tile=(32, 32),
...     compression='zlib',
...     compressionargs={'level': 8},
...     predictor=True,
...     metadata={'axes': 'TZCYX'},
... )

Write a 10 fps time series of volumes with xyz voxel size 2.6755x2.6755x3.9474
micron^3 to an ImageJ hyperstack formatted TIFF file:

>>> volume = numpy.random.randn(6, 57, 256, 256).astype('float32')
>>> image_labels = [f'{i}' for i in range(volume.shape[0] * volume.shape[1])]
>>> imwrite(
...     'temp.tif',
...     volume,
...     imagej=True,
...     resolution=(1.0 / 2.6755, 1.0 / 2.6755),
...     metadata={
...         'spacing': 3.947368,
...         'unit': 'um',
...         'finterval': 1 / 10,
...         'fps': 10.0,
...         'axes': 'TZYX',
...         'Labels': image_labels,
...     },
... )

Read the volume and metadata from the ImageJ hyperstack file:

>>> with TiffFile('temp.tif') as tif:
...     volume = tif.asarray()
...     axes = tif.series[0].axes
...     imagej_metadata = tif.imagej_metadata
...
>>> volume.shape
(6, 57, 256, 256)
>>> axes
'TZYX'
>>> imagej_metadata['slices']
57
>>> imagej_metadata['frames']
6

Memory-map the contiguous image data in the ImageJ hyperstack file:

>>> memmap_volume = memmap('temp.tif')
>>> memmap_volume.shape
(6, 57, 256, 256)
>>> del memmap_volume

Create a TIFF file containing an empty image and write to the memory-mapped
NumPy array (note: this does not work with compression or tiling):

>>> memmap_image = memmap(
...     'temp.tif', shape=(256, 256, 3), dtype='float32', photometric='rgb'
... )
>>> type(memmap_image)
<class 'numpy.memmap'>
>>> memmap_image[255, 255, 1] = 1.0
>>> memmap_image.flush()
>>> del memmap_image

Write two NumPy arrays to a multi-series TIFF file (note: other TIFF readers
will not recognize the two series; use the OME-TIFF format for better
interoperability):

>>> series0 = numpy.random.randint(0, 255, (32, 32, 3), 'uint8')
>>> series1 = numpy.random.randint(0, 255, (4, 256, 256), 'uint16')
>>> with TiffWriter('temp.tif') as tif:
...     tif.write(series0, photometric='rgb')
...     tif.write(series1, photometric='minisblack')
...

Read the second image series from the TIFF file:

>>> series1 = imread('temp.tif', series=1)
>>> series1.shape
(4, 256, 256)

Successively write the frames of one contiguous series to a TIFF file:

>>> data = numpy.random.randint(0, 255, (30, 301, 219), 'uint8')
>>> with TiffWriter('temp.tif') as tif:
...     for frame in data:
...         tif.write(frame, contiguous=True)
...

Append an image series to the existing TIFF file (note: this does not work
with ImageJ hyperstack or OME-TIFF files):

>>> data = numpy.random.randint(0, 255, (301, 219, 3), 'uint8')
>>> imwrite('temp.tif', data, photometric='rgb', append=True)

Create a TIFF file from a generator of tiles:

>>> data = numpy.random.randint(0, 2**12, (31, 33, 3), 'uint16')
>>> def tiles(data, tileshape):
...     for y in range(0, data.shape[0], tileshape[0]):
...         for x in range(0, data.shape[1], tileshape[1]):
...             yield data[y : y + tileshape[0], x : x + tileshape[1]]
...
>>> imwrite(
...     'temp.tif',
...     tiles(data, (16, 16)),
...     tile=(16, 16),
...     shape=data.shape,
...     dtype=data.dtype,
...     photometric='rgb',
... )

Write a multi-dimensional, multi-resolution (pyramidal), multi-series OME-TIFF
file with optional metadata. Sub-resolution images are written to SubIFDs.
Limit parallel encoding to 2 threads. Write a thumbnail image as a separate
image series:

>>> data = numpy.random.randint(0, 255, (8, 2, 512, 512, 3), 'uint8')
>>> subresolutions = 2
>>> pixelsize = 0.29  # micrometer
>>> with TiffWriter('temp.ome.tif', bigtiff=True) as tif:
...     metadata = {
...         'axes': 'TCYXS',
...         'SignificantBits': 8,
...         'TimeIncrement': 0.1,
...         'TimeIncrementUnit': 's',
...         'PhysicalSizeX': pixelsize,
...         'PhysicalSizeXUnit': 'µm',
...         'PhysicalSizeY': pixelsize,
...         'PhysicalSizeYUnit': 'µm',
...         'Channel': {'Name': ['Channel 1', 'Channel 2']},
...         'Plane': {'PositionX': [0.0] * 16, 'PositionXUnit': ['µm'] * 16},
...         'Description': 'A multi-dimensional, multi-resolution image',
...         'MapAnnotation': {  # for OMERO
...             'Namespace': 'openmicroscopy.org/PyramidResolution',
...             '1': '256 256',
...             '2': '128 128',
...         },
...     }
...     options = dict(
...         photometric='rgb',
...         tile=(128, 128),
...         compression='jpeg',
...         resolutionunit='CENTIMETER',
...         maxworkers=2,
...     )
...     tif.write(
...         data,
...         subifds=subresolutions,
...         resolution=(1e4 / pixelsize, 1e4 / pixelsize),
...         metadata=metadata,
...         **options,
...     )
...     # write pyramid levels to the two subifds
...     # in production use resampling to generate sub-resolution images
...     for level in range(subresolutions):
...         mag = 2 ** (level + 1)
...         tif.write(
...             data[..., ::mag, ::mag, :],
...             subfiletype=1,
...             resolution=(1e4 / mag / pixelsize, 1e4 / mag / pixelsize),
...             **options,
...         )
...     # add a thumbnail image as a separate series
...     # it is recognized by QuPath as an associated image
...     thumbnail = (data[0, 0, ::8, ::8] >> 2).astype('uint8')
...     tif.write(thumbnail, metadata={'Name': 'thumbnail'})
...

Access the image levels in the pyramidal OME-TIFF file:

>>> baseimage = imread('temp.ome.tif')
>>> second_level = imread('temp.ome.tif', series=0, level=1)
>>> with TiffFile('temp.ome.tif') as tif:
...     baseimage = tif.series[0].asarray()
...     second_level = tif.series[0].levels[1].asarray()
...     number_levels = len(tif.series[0].levels)  # includes base level
...

Iterate over and decode single JPEG compressed tiles in the TIFF file:

>>> with TiffFile('temp.ome.tif') as tif:
...     fh = tif.filehandle
...     for page in tif.pages:
...         for index, (offset, bytecount) in enumerate(
...             zip(page.dataoffsets, page.databytecounts)
...         ):
...             _ = fh.seek(offset)
...             data = fh.read(bytecount)
...             tile, indices, shape = page.decode(
...                 data, index, jpegtables=page.jpegtables
...             )
...

Use Zarr 2 to read parts of the tiled, pyramidal images in the TIFF file:

>>> import zarr
>>> store = imread('temp.ome.tif', aszarr=True)
>>> z = zarr.open(store, mode='r')
>>> z
<zarr.hierarchy.Group '/' read-only>
>>> z[0]  # base layer
<zarr.core.Array '/0' (8, 2, 512, 512, 3) uint8 read-only>
>>> z[0][2, 0, 128:384, 256:].shape  # read a tile from the base layer
(256, 256, 3)
>>> store.close()

Load the base layer from the Zarr 2 store as a dask array:

>>> import dask.array
>>> store = imread('temp.ome.tif', aszarr=True)
>>> dask.array.from_zarr(store, 0)
dask.array<...shape=(8, 2, 512, 512, 3)...chunksize=(1, 1, 128, 128, 3)...
>>> store.close()

Write the Zarr 2 store to a fsspec ReferenceFileSystem in JSON format:

>>> store = imread('temp.ome.tif', aszarr=True)
>>> store.write_fsspec('temp.ome.tif.json', url='file://')
>>> store.close()

Open the fsspec ReferenceFileSystem as a Zarr group:

>>> import fsspec
>>> import imagecodecs.numcodecs
>>> imagecodecs.numcodecs.register_codecs()
>>> mapper = fsspec.get_mapper(
...     'reference://', fo='temp.ome.tif.json', target_protocol='file'
... )
>>> z = zarr.open(mapper, mode='r')
>>> z
<zarr.hierarchy.Group '/' read-only>

Create an OME-TIFF file containing an empty, tiled image series and write
to it via the Zarr 2 interface (note: this does not work with compression):

>>> imwrite(
...     'temp2.ome.tif',
...     shape=(8, 800, 600),
...     dtype='uint16',
...     photometric='minisblack',
...     tile=(128, 128),
...     metadata={'axes': 'CYX'},
... )
>>> store = imread('temp2.ome.tif', mode='r+', aszarr=True)
>>> z = zarr.open(store, mode='r+')
>>> z
<zarr.core.Array (8, 800, 600) uint16>
>>> z[3, 100:200, 200:300:2] = 1024
>>> store.close()

Read images from a sequence of TIFF files as NumPy array using two I/O worker
threads:

>>> imwrite('temp_C001T001.tif', numpy.random.rand(64, 64))
>>> imwrite('temp_C001T002.tif', numpy.random.rand(64, 64))
>>> image_sequence = imread(
...     ['temp_C001T001.tif', 'temp_C001T002.tif'], ioworkers=2, maxworkers=1
... )
>>> image_sequence.shape
(2, 64, 64)
>>> image_sequence.dtype
dtype('float64')

Read an image stack from a series of TIFF files with a file name pattern
as NumPy or Zarr 2 arrays:

>>> image_sequence = TiffSequence('temp_C0*.tif', pattern=r'_(C)(\d+)(T)(\d+)')
>>> image_sequence.shape
(1, 2)
>>> image_sequence.axes
'CT'
>>> data = image_sequence.asarray()
>>> data.shape
(1, 2, 64, 64)
>>> store = image_sequence.aszarr()
>>> zarr.open(store, mode='r')
<zarr.core.Array (1, 2, 64, 64) float64 read-only>
>>> image_sequence.close()

Write the Zarr 2 store to a fsspec ReferenceFileSystem in JSON format:

>>> store = image_sequence.aszarr()
>>> store.write_fsspec('temp.json', url='file://')

Open the fsspec ReferenceFileSystem as a Zarr 2 array:

>>> import fsspec
>>> import tifffile.numcodecs
>>> tifffile.numcodecs.register_codec()
>>> mapper = fsspec.get_mapper(
...     'reference://', fo='temp.json', target_protocol='file'
... )
>>> zarr.open(mapper, mode='r')
<zarr.core.Array (1, 2, 64, 64) float64 read-only>

Inspect the TIFF file from the command line::

    $ python -m tifffile temp.ome.tif

"""

from __future__ import annotations

__version__ = '2025.5.10'

__all__ = [
    '__version__',
    'TiffFile',
    'TiffFileError',
    'TiffFrame',
    'TiffPage',
    'TiffPages',
    'TiffPageSeries',
    'TiffReader',
    'TiffSequence',
    'TiffTag',
    'TiffTags',
    'TiffTagRegistry',
    'TiffWriter',
    'TiffFormat',
    'ZarrFileSequenceStore',
    'ZarrStore',
    'ZarrTiffStore',
    'imread',
    'imshow',
    'imwrite',
    'lsm2bin',
    'memmap',
    'read_ndtiff_index',
    'read_gdal_structural_metadata',
    'read_micromanager_metadata',
    'read_scanimage_metadata',
    'tiff2fsspec',
    'tiffcomment',
    'TIFF',
    'DATATYPE',
    'CHUNKMODE',
    'COMPRESSION',
    'EXTRASAMPLE',
    'FILETYPE',
    'FILLORDER',
    'OFILETYPE',
    'ORIENTATION',
    'PHOTOMETRIC',
    'PLANARCONFIG',
    'PREDICTOR',
    'RESUNIT',
    'SAMPLEFORMAT',
    'OmeXml',
    'OmeXmlError',
    'FileCache',
    'FileHandle',
    'FileSequence',
    'StoredShape',
    'TiledSequence',
    'NullContext',
    'Timer',
    'askopenfilename',
    'astype',
    'create_output',
    'enumarg',
    'enumstr',
    'format_size',
    'hexdump',
    'imagej_description',
    'imagej_metadata_tag',
    'logger',
    'matlabstr2py',
    'natural_sorted',
    'nullfunc',
    'parse_filenames',
    'parse_kwargs',
    'pformat',
    'product',
    'repeat_nd',
    'reshape_axes',
    'reshape_nd',
    'strptime',
    'transpose_axes',
    'update_kwargs',
    'validate_jhove',
    'xml2dict',
    '_TIFF',  # private
    # deprecated
    'stripnull',
]

import binascii
import collections
import enum
import glob
import io
import json
import logging
import math
import os
import re
import struct
import sys
import threading
import time
import warnings
from collections.abc import (
    Callable,
    Iterable,
    Mapping,
    MutableMapping,
    Sequence,
)
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as DateTime
from datetime import timedelta as TimeDelta
from functools import cached_property

import numpy

try:
    import imagecodecs
except ImportError:
    # load pure Python implementation of some codecs
    try:
        from . import _imagecodecs as imagecodecs  # type: ignore[no-redef]
    except ImportError:
        import _imagecodecs as imagecodecs  # type: ignore[no-redef]

from typing import IO, TYPE_CHECKING, cast, final, overload

if TYPE_CHECKING:
    from collections.abc import (
        Collection,
        Container,
        ItemsView,
        Iterator,
        KeysView,
        ValuesView,
    )
    from typing import Any, Literal, Optional, TextIO, Union

    from numpy.typing import ArrayLike, DTypeLike, NDArray

    ByteOrder = Literal['>', '<']
    OutputType = Union[str, IO[bytes], NDArray[Any], None]
    TagTuple = tuple[
        Union[int, str], Union[int, str], Optional[int], Any, bool
    ]


@overload
def imread(
    files: (
        str
        | os.PathLike[Any]
        | FileHandle
        | IO[bytes]
        | Sequence[str | os.PathLike[Any]]
        | None
    ) = None,
    *,
    selection: Any | None = None,  # TODO: type this
    aszarr: Literal[False] = ...,
    key: int | slice | Iterable[int] | None = None,
    series: int | None = None,
    level: int | None = None,
    squeeze: bool | None = None,
    maxworkers: int | None = None,
    buffersize: int | None = None,
    mode: Literal['r', 'r+'] | None = None,
    name: str | None = None,
    offset: int | None = None,
    size: int | None = None,
    pattern: str | None = None,
    axesorder: Sequence[int] | None = None,
    categories: dict[str, dict[str, int]] | None = None,
    imread: Callable[..., NDArray[Any]] | None = None,
    sort: Callable[..., Any] | bool | None = None,
    container: str | os.PathLike[Any] | None = None,
    chunkshape: tuple[int, ...] | None = None,
    dtype: DTypeLike | None = None,
    axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
    ioworkers: int | None = 1,
    chunkmode: CHUNKMODE | int | str | None = None,
    fillvalue: int | float | None = None,
    zattrs: dict[str, Any] | None = None,
    multiscales: bool | None = None,
    omexml: str | None = None,
    out: OutputType = None,
    out_inplace: bool | None = None,
    _multifile: bool | None = None,
    _useframes: bool | None = None,
    **kwargs: Any,
) -> NDArray[Any]: ...


@overload
def imread(
    files: (
        str
        | os.PathLike[Any]
        | FileHandle
        | IO[bytes]
        | Sequence[str | os.PathLike[Any]]
        | None
    ) = None,
    *,
    selection: Any | None = None,  # TODO: type this
    aszarr: Literal[True],
    key: int | slice | Iterable[int] | None = None,
    series: int | None = None,
    level: int | None = None,
    squeeze: bool | None = None,
    maxworkers: int | None = None,
    buffersize: int | None = None,
    mode: Literal['r', 'r+'] | None = None,
    name: str | None = None,
    offset: int | None = None,
    size: int | None = None,
    pattern: str | None = None,
    axesorder: Sequence[int] | None = None,
    categories: dict[str, dict[str, int]] | None = None,
    imread: Callable[..., NDArray[Any]] | None = None,
    imreadargs: dict[str, Any] | None = None,
    sort: Callable[..., Any] | bool | None = None,
    container: str | os.PathLike[Any] | None = None,
    chunkshape: tuple[int, ...] | None = None,
    chunkdtype: DTypeLike | None = None,
    axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
    ioworkers: int | None = 1,
    chunkmode: CHUNKMODE | int | str | None = None,
    fillvalue: int | float | None = None,
    zattrs: dict[str, Any] | None = None,
    multiscales: bool | None = None,
    omexml: str | None = None,
    out: OutputType = None,
    out_inplace: bool | None = None,
    _multifile: bool | None = None,
    _useframes: bool | None = None,
    **kwargs: Any,
) -> ZarrTiffStore | ZarrFileSequenceStore: ...


@overload
def imread(
    files: (
        str
        | os.PathLike[Any]
        | FileHandle
        | IO[bytes]
        | Sequence[str | os.PathLike[Any]]
        | None
    ) = None,
    *,
    selection: Any | None = None,  # TODO: type this
    aszarr: bool = False,
    key: int | slice | Iterable[int] | None = None,
    series: int | None = None,
    level: int | None = None,
    squeeze: bool | None = None,
    maxworkers: int | None = None,
    buffersize: int | None = None,
    mode: Literal['r', 'r+'] | None = None,
    name: str | None = None,
    offset: int | None = None,
    size: int | None = None,
    pattern: str | None = None,
    axesorder: Sequence[int] | None = None,
    categories: dict[str, dict[str, int]] | None = None,
    imread: Callable[..., NDArray[Any]] | None = None,
    imreadargs: dict[str, Any] | None = None,
    sort: Callable[..., Any] | bool | None = None,
    container: str | os.PathLike[Any] | None = None,
    chunkshape: tuple[int, ...] | None = None,
    chunkdtype: DTypeLike | None = None,
    axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
    ioworkers: int | None = 1,
    chunkmode: CHUNKMODE | int | str | None = None,
    fillvalue: int | float | None = None,
    zattrs: dict[str, Any] | None = None,
    multiscales: bool | None = None,
    omexml: str | None = None,
    out: OutputType = None,
    out_inplace: bool | None = None,
    _multifile: bool | None = None,
    _useframes: bool | None = None,
    **kwargs: Any,
) -> NDArray[Any] | ZarrTiffStore | ZarrFileSequenceStore: ...


def imread(
    files: (
        str
        | os.PathLike[Any]
        | FileHandle
        | IO[bytes]
        | Sequence[str | os.PathLike[Any]]
        | None
    ) = None,
    *,
    selection: Any | None = None,  # TODO: type this
    aszarr: bool = False,
    key: int | slice | Iterable[int] | None = None,
    series: int | None = None,
    level: int | None = None,
    squeeze: bool | None = None,
    maxworkers: int | None = None,
    buffersize: int | None = None,
    mode: Literal['r', 'r+'] | None = None,
    name: str | None = None,
    offset: int | None = None,
    size: int | None = None,
    pattern: str | None = None,
    axesorder: Sequence[int] | None = None,
    categories: dict[str, dict[str, int]] | None = None,
    imread: Callable[..., NDArray[Any]] | None = None,
    imreadargs: dict[str, Any] | None = None,
    sort: Callable[..., Any] | bool | None = None,
    container: str | os.PathLike[Any] | None = None,
    chunkshape: tuple[int, ...] | None = None,
    chunkdtype: DTypeLike | None = None,
    axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
    ioworkers: int | None = 1,
    chunkmode: CHUNKMODE | int | str | None = None,
    fillvalue: int | float | None = None,
    zattrs: dict[str, Any] | None = None,
    multiscales: bool | None = None,
    omexml: str | None = None,
    out: OutputType = None,
    out_inplace: bool | None = None,
    _multifile: bool | None = None,
    _useframes: bool | None = None,
    **kwargs: Any,
) -> NDArray[Any] | ZarrTiffStore | ZarrFileSequenceStore:
    """Return image from TIFF file(s) as NumPy array or Zarr 2 store.

    The first image series in the file(s) is returned by default.

    Parameters:
        files:
            File name, seekable binary stream, glob pattern, or sequence of
            file names. May be *None* if `container` is specified.
        selection:
            Subset of image to be extracted.
            If not None, a Zarr 2 array is created, indexed with the
            `selection` value, and returned as a NumPy array. Only segments
            that are part of the selection will be read from file.
            Refer to the Zarr 2 documentation for valid selections.
            Depending on selection size, image size, and storage properties,
            it may be more efficient to read the whole image from file and
            then index it.
        aszarr:
            Return file sequences, series, or single pages as Zarr 2 store
            instead of NumPy array if `selection` is None.
        mode, name, offset, size, omexml, _multifile, _useframes:
            Passed to :py:class:`TiffFile`.
        key, series, level, squeeze, maxworkers, buffersize:
            Passed to :py:meth:`TiffFile.asarray`
            or :py:meth:`TiffFile.aszarr`.
        imread, container, sort, pattern, axesorder, axestiled, categories,\
        ioworkers:
            Passed to :py:class:`FileSequence`.
        chunkmode, fillvalue, zattrs, multiscales:
            Passed to :py:class:`ZarrTiffStore`
            or :py:class:`ZarrFileSequenceStore`.
        chunkshape, chunkdtype:
            Passed to :py:meth:`FileSequence.asarray` or
            :py:class:`ZarrFileSequenceStore`.
        out_inplace:
            Passed to :py:meth:`FileSequence.asarray`
        out:
            Passed to :py:meth:`TiffFile.asarray`,
            :py:meth:`FileSequence.asarray`, or :py:func:`zarr_selection`.
        imreadargs:
            Additional arguments passed to :py:attr:`FileSequence.imread`.
        **kwargs:
            Additional arguments passed to :py:class:`TiffFile` or
            :py:attr:`FileSequence.imread`.

    Returns:
        Images from specified files, series, or pages.
        Zarr 2 store instances must be closed after use.
        See :py:meth:`TiffPage.asarray` for operations that are applied
        (or not) to the image data stored in the file.

    """
    store: ZarrStore
    aszarr = aszarr or (selection is not None)
    is_flags = parse_kwargs(kwargs, *(k for k in kwargs if k[:3] == 'is_'))

    if imread is None and kwargs:
        raise TypeError(
            'imread() got unexpected keyword arguments '
            + ', '.join(f"'{key}'" for key in kwargs)
        )

    if container is None:
        if isinstance(files, str) and ('*' in files or '?' in files):
            files = glob.glob(files)
        if not files:
            raise ValueError('no files found')

        if (
            isinstance(files, Sequence)
            and not isinstance(files, str)
            and len(files) == 1
        ):
            files = files[0]

        if isinstance(files, str) or not isinstance(files, Sequence):
            with TiffFile(
                files,
                mode=mode,
                name=name,
                offset=offset,
                size=size,
                omexml=omexml,
                _multifile=_multifile,
                _useframes=_useframes,
                **is_flags,
            ) as tif:
                if aszarr:
                    assert key is None or isinstance(key, int)
                    store = tif.aszarr(
                        key=key,
                        series=series,
                        level=level,
                        squeeze=squeeze,
                        maxworkers=maxworkers,
                        buffersize=buffersize,
                        chunkmode=chunkmode,
                        fillvalue=fillvalue,
                        zattrs=zattrs,
                        multiscales=multiscales,
                    )
                    if selection is None:
                        return store
                    return zarr_selection(store, selection, out=out)
                return tif.asarray(
                    key=key,
                    series=series,
                    level=level,
                    squeeze=squeeze,
                    maxworkers=maxworkers,
                    buffersize=buffersize,
                    out=out,
                )

    elif isinstance(files, (FileHandle, IO)):
        raise ValueError('BinaryIO not supported')

    imread_kwargs = kwargs_notnone(
        key=key,
        series=series,
        level=level,
        squeeze=squeeze,
        maxworkers=maxworkers,
        buffersize=buffersize,
        imreadargs=imreadargs,
        _multifile=_multifile,
        _useframes=_useframes,
        **is_flags,
        **kwargs,
    )

    with TiffSequence(
        files,
        pattern=pattern,
        axesorder=axesorder,
        categories=categories,
        container=container,
        sort=sort,
        **kwargs_notnone(imread=imread),
    ) as imseq:
        if aszarr:
            store = imseq.aszarr(
                axestiled=axestiled,
                chunkmode=chunkmode,
                chunkshape=chunkshape,
                chunkdtype=chunkdtype,
                fillvalue=fillvalue,
                zattrs=zattrs,
                **imread_kwargs,
            )
            if selection is None:
                return store
            return zarr_selection(store, selection, out=out)
        return imseq.asarray(
            axestiled=axestiled,
            chunkshape=chunkshape,
            chunkdtype=chunkdtype,
            ioworkers=ioworkers,
            out=out,
            out_inplace=out_inplace,
            **imread_kwargs,
        )


def imwrite(
    file: str | os.PathLike[Any] | FileHandle | IO[bytes],
    /,
    data: (
        ArrayLike | Iterator[NDArray[Any] | None] | Iterator[bytes] | None
    ) = None,
    *,
    mode: Literal['w', 'x', 'r+'] | None = None,
    bigtiff: bool | None = None,
    byteorder: ByteOrder | None = None,
    imagej: bool = False,
    ome: bool | None = None,
    shaped: bool | None = None,
    append: bool = False,
    shape: Sequence[int] | None = None,
    dtype: DTypeLike | None = None,
    photometric: PHOTOMETRIC | int | str | None = None,
    planarconfig: PLANARCONFIG | int | str | None = None,
    extrasamples: Sequence[EXTRASAMPLE | int | str] | None = None,
    volumetric: bool = False,
    tile: Sequence[int] | None = None,
    rowsperstrip: int | None = None,
    bitspersample: int | None = None,
    compression: COMPRESSION | int | str | None = None,
    compressionargs: dict[str, Any] | None = None,
    predictor: PREDICTOR | int | str | bool | None = None,
    subsampling: tuple[int, int] | None = None,
    jpegtables: bytes | None = None,
    iccprofile: bytes | None = None,
    colormap: ArrayLike | None = None,
    description: str | bytes | None = None,
    datetime: str | bool | DateTime | None = None,
    resolution: (
        tuple[float | tuple[int, int], float | tuple[int, int]] | None
    ) = None,
    resolutionunit: RESUNIT | int | str | None = None,
    subfiletype: FILETYPE | int | None = None,
    software: str | bytes | bool | None = None,
    # subifds: int | Sequence[int] | None = None,
    metadata: dict[str, Any] | None = {},
    extratags: Sequence[TagTuple] | None = None,
    contiguous: bool = False,
    truncate: bool = False,
    align: int | None = None,
    maxworkers: int | None = None,
    buffersize: int | None = None,
    returnoffset: bool = False,
) -> tuple[int, int] | None:
    """Write NumPy array to TIFF file.

    A BigTIFF file is written if the data size is larger than 4 GB less
    32 MB for metadata, and `bigtiff` is not *False*, and `imagej`,
    `truncate` and `compression` are not enabled.
    Unless `byteorder` is specified, the TIFF file byte order is determined
    from the dtype of `data` or the `dtype` argument.

    Parameters:
        file:
            Passed to :py:class:`TiffWriter`.
        data, shape, dtype:
            Passed to :py:meth:`TiffWriter.write`.
        mode, append, byteorder, bigtiff, imagej, ome, shaped:
            Passed to :py:class:`TiffWriter`.
        photometric, planarconfig, extrasamples, volumetric, tile,\
        rowsperstrip, bitspersample, compression, compressionargs, predictor,\
        subsampling, jpegtables, iccprofile, colormap, description, datetime,\
        resolution, resolutionunit, subfiletype, software,\
        metadata, extratags, maxworkers, buffersize, \
        contiguous, truncate, align:
            Passed to :py:meth:`TiffWriter.write`.
        returnoffset:
            Return offset and number of bytes of memory-mappable image data
            in file.

    Returns:
        If `returnoffset` is *True* and the image data in the file are
        memory-mappable, the offset and number of bytes of the image
        data in the file.

    """
    if data is None:
        # write empty file
        if shape is None or dtype is None:
            raise ValueError("missing required 'shape' or 'dtype' argument")
        dtype = numpy.dtype(dtype)
        shape = tuple(shape)
        datasize = product(shape) * dtype.itemsize
        if byteorder is None:
            byteorder = dtype.byteorder  # type: ignore[assignment]
    else:
        #
        try:
            datasize = data.nbytes  # type: ignore[union-attr]
            if byteorder is None:
                byteorder = data.dtype.byteorder  # type: ignore[union-attr]
        except Exception:
            datasize = 0

    if bigtiff is None:
        bigtiff = (
            datasize > 2**32 - 2**25
            and not imagej
            and not truncate
            and compression in {None, 0, 1, 'NONE', 'None', 'none'}
        )

    with TiffWriter(
        file,
        mode=mode,
        bigtiff=bigtiff,
        byteorder=byteorder,
        append=append,
        imagej=imagej,
        ome=ome,
        shaped=shaped,
    ) as tif:
        result = tif.write(
            data,
            shape=shape,
            dtype=dtype,
            photometric=photometric,
            planarconfig=planarconfig,
            extrasamples=extrasamples,
            volumetric=volumetric,
            tile=tile,
            rowsperstrip=rowsperstrip,
            bitspersample=bitspersample,
            compression=compression,
            compressionargs=compressionargs,
            predictor=predictor,
            subsampling=subsampling,
            jpegtables=jpegtables,
            iccprofile=iccprofile,
            colormap=colormap,
            description=description,
            datetime=datetime,
            resolution=resolution,
            resolutionunit=resolutionunit,
            subfiletype=subfiletype,
            software=software,
            metadata=metadata,
            extratags=extratags,
            contiguous=contiguous,
            truncate=truncate,
            align=align,
            maxworkers=maxworkers,
            buffersize=buffersize,
            returnoffset=returnoffset,
        )
    return result


def memmap(
    filename: str | os.PathLike[Any],
    /,
    *,
    shape: Sequence[int] | None = None,
    dtype: DTypeLike | None = None,
    page: int | None = None,
    series: int = 0,
    level: int = 0,
    mode: Literal['r+', 'r', 'c'] = 'r+',
    **kwargs: Any,
) -> numpy.memmap[Any, Any]:
    """Return memory-mapped NumPy array of image data stored in TIFF file.

    Memory-mapping requires the image data stored in native byte order,
    without tiling, compression, predictors, etc.
    If `shape` and `dtype` are provided, existing files are overwritten or
    appended to depending on the `append` argument.
    Else, the image data of a specified page or series in an existing
    file are memory-mapped. By default, the image data of the first
    series are memory-mapped.
    Call `flush` to write any changes in the array to the file.

    Parameters:
        filename:
            Name of TIFF file which stores array.
        shape:
            Shape of empty array.
        dtype:
            Datatype of empty array.
        page:
            Index of page which image data to memory-map.
        series:
            Index of page series which image data to memory-map.
        level:
            Index of pyramid level which image data to memory-map.
        mode:
            Memory-map file open mode. The default is 'r+', which opens
            existing file for reading and writing.
        **kwargs:
            Additional arguments passed to :py:func:`imwrite` or
            :py:class:`TiffFile`.

    Returns:
        Image in TIFF file as memory-mapped NumPy array.

    Raises:
        ValueError: Image data in TIFF file are not memory-mappable.

    """
    filename = os.fspath(filename)
    if shape is not None:
        shape = tuple(shape)
    if shape is not None and dtype is not None:
        # create a new, empty array
        dtype = numpy.dtype(dtype)
        if 'byteorder' in kwargs:
            dtype = dtype.newbyteorder(kwargs['byteorder'])
        kwargs.update(
            data=None,
            shape=shape,
            dtype=dtype,
            align=TIFF.ALLOCATIONGRANULARITY,
            returnoffset=True,
        )
        result = imwrite(filename, **kwargs)
        if result is None:
            # TODO: fail before creating file or writing data
            raise ValueError('image data are not memory-mappable')
        offset = result[0]
    else:
        # use existing file
        with TiffFile(filename, **kwargs) as tif:
            if page is None:
                tiffseries = tif.series[series].levels[level]
                if tiffseries.dataoffset is None:
                    raise ValueError('image data are not memory-mappable')
                shape = tiffseries.shape
                dtype = tiffseries.dtype
                offset = tiffseries.dataoffset
            else:
                tiffpage = tif.pages[page]
                if not tiffpage.is_memmappable:
                    raise ValueError('image data are not memory-mappable')
                offset = tiffpage.dataoffsets[0]
                shape = tiffpage.shape
                dtype = tiffpage.dtype
                assert dtype is not None
            dtype = numpy.dtype(tif.byteorder + dtype.char)
    return numpy.memmap(filename, dtype, mode, offset, shape, 'C')


class TiffFileError(Exception):
    """Exception to indicate invalid TIFF structure."""


@final
class TiffWriter:
    """Write NumPy arrays to TIFF file.

    TiffWriter's main purpose is saving multi-dimensional NumPy arrays in
    TIFF containers, not to create any possible TIFF format.
    Specifically, ExifIFD and GPSIFD tags are not supported.

    TiffWriter instances must be closed with :py:meth:`TiffWriter.close`,
    which is automatically called when using the 'with' context manager.

    TiffWriter instances are not thread-safe. All attributes are read-only.

    Parameters:
        file:
            Specifies file to write.
        mode:
            Binary file open mode if `file` is file name.
            The default is 'w', which opens files for writing, truncating
            existing files.
            'x' opens files for exclusive creation, failing on existing files.
            'r+' opens files for updating, enabling `append`.
        bigtiff:
            Write 64-bit BigTIFF formatted file, which can exceed 4 GB.
            By default, a classic 32-bit TIFF file is written, which is
            limited to 4 GB.
            If `append` is *True*, the existing file's format is used.
        byteorder:
            Endianness of TIFF format. One of '<', '>', '=', or '|'.
            The default is the system's native byte order.
        append:
            If `file` is existing standard TIFF file, append image data
            and tags to file.
            Parameters `bigtiff` and `byteorder` set from existing file.
            Appending does not scale well with the number of pages already in
            the file and may corrupt specifically formatted TIFF files such as
            OME-TIFF, LSM, STK, ImageJ, or FluoView.
        imagej:
            Write ImageJ hyperstack compatible file if `ome` is not enabled.
            This format can handle data types uint8, uint16, or float32 and
            data shapes up to 6 dimensions in TZCYXS order.
            RGB images (S=3 or S=4) must be `uint8`.
            ImageJ's default byte order is big-endian, but this
            implementation uses the system's native byte order by default.
            ImageJ hyperstacks do not support BigTIFF or compression.
            The ImageJ file format is undocumented.
            Use FIJI's Bio-Formats import function for compressed files.
        ome:
            Write OME-TIFF compatible file.
            By default, the OME-TIFF format is used if the file name extension
            contains '.ome.', `imagej` is not enabled, and the `description`
            argument in the first call of :py:meth:`TiffWriter.write` is not
            specified.
            The format supports multiple, up to 9 dimensional image series.
            The default axes order is TZC(S)YX(S).
            Refer to the OME model for restrictions of this format.
        shaped:
            Write tifffile "shaped" compatible file.
            The shape of multi-dimensional images is stored in JSON format in
            a ImageDescription tag of the first page of a series.
            This is the default format used by tifffile unless `imagej` or
            `ome` are enabled or ``metadata=None`` is passed to
            :py:meth:`TiffWriter.write`.

    Raises:
        ValueError:
            The TIFF file cannot be appended to. Use ``append='force'`` to
            force appending, which may result in a corrupted file.

    """

    tiff: TiffFormat
    """Format of TIFF file being written."""

    _fh: FileHandle
    _omexml: OmeXml | None
    _ome: bool | None  # writing OME-TIFF format
    _imagej: bool  # writing ImageJ format
    _tifffile: bool  # writing Tifffile shaped format
    _truncate: bool
    _metadata: dict[str, Any] | None
    _colormap: NDArray[numpy.uint16] | None
    _tags: list[tuple[int, bytes, Any, bool]] | None
    _datashape: tuple[int, ...] | None  # shape of data in consecutive pages
    _datadtype: numpy.dtype[Any] | None  # data type
    _dataoffset: int | None  # offset to data
    _databytecounts: list[int] | None  # byte counts per plane
    _dataoffsetstag: int | None  # strip or tile offset tag code
    _descriptiontag: TiffTag | None  # TiffTag for updating comment
    _ifdoffset: int
    _subifds: int  # number of subifds
    _subifdslevel: int  # index of current subifd level
    _subifdsoffsets: list[int]  # offsets to offsets to subifds
    _nextifdoffsets: list[int]  # offsets to offset to next ifd
    _ifdindex: int  # index of current ifd
    _storedshape: StoredShape | None  # normalized shape in consecutive pages

    def __init__(
        self,
        file: str | os.PathLike[Any] | FileHandle | IO[bytes],
        /,
        *,
        mode: Literal['w', 'x', 'r+'] | None = None,
        bigtiff: bool = False,
        byteorder: ByteOrder | None = None,
        append: bool | str = False,
        imagej: bool = False,
        ome: bool | None = None,
        shaped: bool | None = None,
    ) -> None:
        if mode in {'r+', 'r+b'} or (
            isinstance(file, FileHandle) and file._mode == 'r+b'
        ):
            mode = 'r+'
            append = True
        if append:
            # determine if file is an existing TIFF file that can be extended
            try:
                with FileHandle(file, mode='rb', size=0) as fh:
                    pos = fh.tell()
                    try:
                        with TiffFile(fh) as tif:
                            if append != 'force' and not tif.is_appendable:
                                raise ValueError(
                                    'cannot append to file containing metadata'
                                )
                            byteorder = tif.byteorder
                            bigtiff = tif.is_bigtiff
                            self._ifdoffset = cast(
                                int, tif.pages.next_page_offset
                            )
                    finally:
                        fh.seek(pos)
                    append = True
            except (OSError, FileNotFoundError):
                append = False

        if append:
            if mode not in {None, 'r+', 'r+b'}:
                raise ValueError("append mode must be 'r+'")
            mode = 'r+'
        elif mode is None:
            mode = 'w'

        if byteorder is None or byteorder in {'=', '|'}:
            byteorder = '<' if sys.byteorder == 'little' else '>'
        elif byteorder not in {'<', '>'}:
            raise ValueError(f'invalid byteorder {byteorder}')

        if byteorder == '<':
            self.tiff = TIFF.BIG_LE if bigtiff else TIFF.CLASSIC_LE
        else:
            self.tiff = TIFF.BIG_BE if bigtiff else TIFF.CLASSIC_BE

        self._truncate = False
        self._metadata = None
        self._colormap = None
        self._tags = None
        self._datashape = None
        self._datadtype = None
        self._dataoffset = None
        self._databytecounts = None
        self._dataoffsetstag = None
        self._descriptiontag = None
        self._subifds = 0
        self._subifdslevel = -1
        self._subifdsoffsets = []
        self._nextifdoffsets = []
        self._ifdindex = 0
        self._omexml = None
        self._storedshape = None

        self._fh = FileHandle(file, mode=mode, size=0)
        if append:
            self._fh.seek(0, os.SEEK_END)
        else:
            assert byteorder is not None
            self._fh.write(b'II' if byteorder == '<' else b'MM')
            if bigtiff:
                self._fh.write(struct.pack(byteorder + 'HHH', 43, 8, 0))
            else:
                self._fh.write(struct.pack(byteorder + 'H', 42))
            # first IFD
            self._ifdoffset = self._fh.tell()
            self._fh.write(struct.pack(self.tiff.offsetformat, 0))

        self._ome = None if ome is None else bool(ome)
        self._imagej = False if self._ome else bool(imagej)
        if self._imagej:
            self._ome = False
        if self._ome or self._imagej:
            self._tifffile = False
        else:
            self._tifffile = True if shaped is None else bool(shaped)

        if imagej and bigtiff:
            warnings.warn(
                f'{self!r} writing nonconformant BigTIFF ImageJ', UserWarning
            )

    def write(
        self,
        data: (
            ArrayLike | Iterator[NDArray[Any] | None] | Iterator[bytes] | None
        ) = None,
        *,
        shape: Sequence[int] | None = None,
        dtype: DTypeLike | None = None,
        photometric: PHOTOMETRIC | int | str | None = None,
        planarconfig: PLANARCONFIG | int | str | None = None,
        extrasamples: Sequence[EXTRASAMPLE | int | str] | None = None,
        volumetric: bool = False,
        tile: Sequence[int] | None = None,
        rowsperstrip: int | None = None,
        bitspersample: int | None = None,
        compression: COMPRESSION | int | str | bool | None = None,
        compressionargs: dict[str, Any] | None = None,
        predictor: PREDICTOR | int | str | bool | None = None,
        subsampling: tuple[int, int] | None = None,
        jpegtables: bytes | None = None,
        iccprofile: bytes | None = None,
        colormap: ArrayLike | None = None,
        description: str | bytes | None = None,
        datetime: str | bool | DateTime | None = None,
        resolution: (
            tuple[float | tuple[int, int], float | tuple[int, int]] | None
        ) = None,
        resolutionunit: RESUNIT | int | str | None = None,
        subfiletype: FILETYPE | int | None = None,
        software: str | bytes | bool | None = None,
        subifds: int | Sequence[int] | None = None,
        metadata: dict[str, Any] | None = {},
        extratags: Sequence[TagTuple] | None = None,
        contiguous: bool = False,
        truncate: bool = False,
        align: int | None = None,
        maxworkers: int | None = None,
        buffersize: int | None = None,
        returnoffset: bool = False,
    ) -> tuple[int, int] | None:
        r"""Write multi-dimensional image to series of TIFF pages.

        Metadata in JSON, ImageJ, or OME-XML format are written to the
        ImageDescription tag of the first page of a series by default,
        such that the image can later be read back as an array of the
        same shape.

        The values of the ImageWidth, ImageLength, ImageDepth, and
        SamplesPerPixel tags are inferred from the last dimensions of the
        data's shape.
        The value of the SampleFormat tag is inferred from the data's dtype.
        Image data are written uncompressed in one strip per plane by default.
        Dimensions higher than 2 to 4 (depending on photometric mode, planar
        configuration, and volumetric mode) are flattened and written as
        separate pages.
        If the data size is zero, write a single page with shape (0, 0).

        Parameters:
            data:
                Specifies image to write.
                If *None*, an empty image is written, which size and type must
                be specified using `shape` and `dtype` arguments.
                This option cannot be used with compression, predictors,
                packed integers, or bilevel images.
                A copy of array-like data is made if it is not a C-contiguous
                numpy or dask array with the same byteorder as the TIFF file.
                Iterators must yield ndarrays or bytes compatible with the
                file's byteorder as well as the `shape` and `dtype` arguments.
                Iterator bytes must be compatible with the `compression`,
                `predictor`, `subsampling`, and `jpegtables` arguments.
                If `tile` is specified, iterator items must match the tile
                shape. Incomplete tiles are zero-padded.
                Iterators of non-tiled images must yield ndarrays of
                `shape[1:]` or strips as bytes. Iterators of strip ndarrays
                are not supported.
                Writing dask arrays might be excruciatingly slow for arrays
                with many chunks or files with many segments.
                (https://github.com/dask/dask/issues/8570).
            shape:
                Shape of image to write.
                The default is inferred from the `data` argument if possible.
                A ValueError is raised if the value is incompatible with
                the `data` or other arguments.
            dtype:
                NumPy data type of image to write.
                The default is inferred from the `data` argument if possible.
                A ValueError is raised if the value is incompatible with
                the `data` argument.
            photometric:
                Color space of image.
                The default is inferred from the data shape, dtype, and the
                `colormap` argument.
                A UserWarning is logged if RGB color space is auto-detected.
                Specify this parameter to silence the warning and to avoid
                ambiguities.
                *MINISBLACK*: for bilevel and grayscale images, 0 is black.
                *MINISWHITE*: for bilevel and grayscale images, 0 is white.
                *RGB*: the image contains red, green and blue samples.
                *SEPARATED*: the image contains CMYK samples.
                *PALETTE*: the image is used as an index into a colormap.
                *CFA*: the image is a Color Filter Array. The
                CFARepeatPatternDim, CFAPattern, and other DNG or TIFF/EP tags
                must be specified in `extratags` to produce a valid file.
                The value is written to the PhotometricInterpretation tag.
            planarconfig:
                Specifies if samples are stored interleaved or in separate
                planes.
                *CONTIG*: the last dimension contains samples.
                *SEPARATE*: the 3rd or 4th last dimension contains samples.
                The default is inferred from the data shape and `photometric`
                mode.
                If this parameter is set, extra samples are used to store
                grayscale images.
                The value is written to the PlanarConfiguration tag.
            extrasamples:
                Interpretation of extra components in pixels.
                *UNSPECIFIED*: no transparency information (default).
                *ASSOCALPHA*: true transparency with premultiplied color.
                *UNASSALPHA*: independent transparency masks.
                The values are written to the ExtraSamples tag.
            volumetric:
                Write volumetric image to single page (instead of multiple
                pages) using SGI ImageDepth tag.
                The volumetric format is not part of the TIFF specification,
                and few software can read it.
                OME and ImageJ formats are not compatible with volumetric
                storage.
            tile:
                Shape ([depth,] length, width) of image tiles to write.
                By default, image data are written in strips.
                The tile length and width must be a multiple of 16.
                If a tile depth is provided, the SGI ImageDepth and TileDepth
                tags are used to write volumetric data.
                Tiles cannot be used to write contiguous series, except if
                the tile shape matches the data shape.
                The values are written to the TileWidth, TileLength, and
                TileDepth tags.
            rowsperstrip:
                Number of rows per strip.
                By default, strips are about 256 KB if `compression` is
                enabled, else rowsperstrip is set to the image length.
                The value is written to the RowsPerStrip tag.
            bitspersample:
                Number of bits per sample.
                The default is the number of bits of the data's dtype.
                Different values per samples are not supported.
                Unsigned integer data are packed into bytes as tightly as
                possible.
                Valid values are 1-8 for uint8, 9-16 for uint16, and 17-32
                for uint32.
                This setting cannot be used with compression, contiguous
                series, or empty files.
                The value is written to the BitsPerSample tag.
            compression:
                Compression scheme used on image data.
                By default, image data are written uncompressed.
                Compression cannot be used to write contiguous series.
                Compressors may require certain data shapes, types or value
                ranges. For example, JPEG compression requires grayscale or
                RGB(A), uint8 or 12-bit uint16.
                JPEG compression is experimental. JPEG markers and TIFF tags
                may not match.
                Only a limited set of compression schemes are implemented.
                'ZLIB' is short for ADOBE_DEFLATE.
                The value is written to the Compression tag.
            compressionargs:
                Extra arguments passed to compression codec, for example,
                compression level. Refer to the Imagecodecs implementation
                for supported arguments.
            predictor:
                Horizontal differencing operator applied to image data before
                compression.
                By default, no operator is applied.
                Predictors can only be used with certain compression schemes
                and data types.
                The value is written to the Predictor tag.
            subsampling:
                Horizontal and vertical subsampling factors used for the
                chrominance components of images: (1, 1), (2, 1), (2, 2), or
                (4, 1). The default is *(2, 2)*.
                Currently applies to JPEG compression of RGB images only.
                Images are stored in YCbCr color space, the value of the
                PhotometricInterpretation tag is *YCBCR*.
                Segment widths must be a multiple of 8 times the horizontal
                factor. Segment lengths and rowsperstrip must be a multiple
                of 8 times the vertical factor.
                The values are written to the YCbCrSubSampling tag.
            jpegtables:
                JPEG quantization and/or Huffman tables.
                Use for copying pre-compressed JPEG segments.
                The value is written to the JPEGTables tag.
            iccprofile:
                International Color Consortium (ICC) device profile
                characterizing image color space.
                The value is written verbatim to the InterColorProfile tag.
            colormap:
                RGB color values for corresponding data value.
                The colormap array must be of shape
                `(3, 2\*\*(data.itemsize*8))` (or `(3, 256)` for ImageJ)
                and dtype uint16.
                The image's data type must be uint8 or uint16 (or float32
                for ImageJ) and the values are indices into the last
                dimension of the colormap.
                The value is written to the ColorMap tag.
            description:
                Subject of image. Must be 7-bit ASCII.
                Cannot be used with the ImageJ or OME formats.
                The value is written to the ImageDescription tag of the
                first page of a series.
            datetime:
                Date and time of image creation in ``%Y:%m:%d %H:%M:%S``
                format or datetime object.
                If *True*, the current date and time is used.
                The value is written to the DateTime tag of the first page
                of a series.
            resolution:
                Number of pixels per `resolutionunit` in X and Y directions
                as float or rational numbers.
                The default is (1.0, 1.0).
                The values are written to the YResolution and XResolution tags.
            resolutionunit:
                Unit of measurement for `resolution` values.
                The default is *NONE* if `resolution` is not specified and
                for ImageJ format, else *INCH*.
                The value is written to the ResolutionUnit tags.
            subfiletype:
                Bitfield to indicate kind of image.
                Set bit 0 if the image is a reduced-resolution version of
                another image.
                Set bit 1 if the image is part of a multi-page image.
                Set bit 2 if the image is transparency mask for another
                image (photometric must be MASK, SamplesPerPixel and
                bitspersample must be 1).
            software:
                Name of software used to create file.
                Must be 7-bit ASCII. The default is 'tifffile.py'.
                Unless *False*, the value is written to the Software tag of
                the first page of a series.
            subifds:
                Number of child IFDs.
                If greater than 0, the following `subifds` number of series
                are written as child IFDs of the current series.
                The number of IFDs written for each SubIFD level must match
                the number of IFDs written for the current series.
                All pages written to a certain SubIFD level of the current
                series must have the same hash.
                SubIFDs cannot be used with truncated or ImageJ files.
                SubIFDs in OME-TIFF files must be sub-resolutions of the
                main IFDs.
            metadata:
                Additional metadata describing image, written along
                with shape information in JSON, OME-XML, or ImageJ formats
                in ImageDescription or IJMetadata tags.
                If *None*, or the `shaped` argument to :py:class:`TiffWriter`
                is *False*, no information in JSON format is written to
                the ImageDescription tag.
                The 'axes' item defines the character codes for dimensions in
                `data` or `shape`.
                Refer to :py:class:`OmeXml` for supported keys when writing
                OME-TIFF.
                Refer to :py:func:`imagej_description` and
                :py:func:`imagej_metadata_tag` for items supported
                by the ImageJ format. Items 'Info', 'Labels', 'Ranges',
                'LUTs', 'Plot', 'ROI', and 'Overlays' are written to the
                IJMetadata and IJMetadataByteCounts tags.
                Strings must be 7-bit ASCII.
                Written with the first page of a series only.
            extratags:
                Additional tags to write. A list of tuples with 5 items:

                0. code (int): Tag Id.
                1. dtype (:py:class:`DATATYPE`):
                   Data type of items in `value`.
                2. count (int): Number of data values.
                   Not used for string or bytes values.
                3. value (Sequence[Any]): `count` values compatible with
                   `dtype`. Bytes must contain count values of dtype packed
                   as binary data.
                4. writeonce (bool): If *True*, write tag to first page
                   of a series only.

                Duplicate and select tags in TIFF.TAG_FILTERED are not written
                if the extratag is specified by integer code.
                Extratags cannot be used to write IFD type tags.

            contiguous:
                If *False* (default), write data to a new series.
                If *True* and the data and arguments are compatible with
                previous written ones (same shape, no compression, etc.),
                the image data are stored contiguously after the previous one.
                In that case, `photometric`, `planarconfig`, and
                `rowsperstrip` are ignored.
                Metadata such as `description`, `metadata`, `datetime`,
                and `extratags` are written to the first page of a contiguous
                series only.
                Contiguous mode cannot be used with the OME or ImageJ formats.
            truncate:
                If *True*, only write first page of contiguous series
                if possible (uncompressed, contiguous, not tiled).
                Other TIFF readers will only be able to read part of the data.
                Cannot be used with the OME or ImageJ formats.
            align:
                Byte boundary on which to align image data in file.
                The default is 16.
                Use mmap.ALLOCATIONGRANULARITY for memory-mapped data.
                Following contiguous writes are not aligned.
            maxworkers:
                Maximum number of threads to concurrently compress tiles
                or strips.
                If *None* or *0*, use up to :py:attr:`_TIFF.MAXWORKERS` CPU
                cores for compressing large segments.
                Using multiple threads can significantly speed up this
                function if the bottleneck is encoding the data, for example,
                in case of large JPEG compressed tiles.
                If the bottleneck is I/O or pure Python code, using multiple
                threads might be detrimental.
            buffersize:
                Approximate number of bytes to compress in one pass.
                The default is :py:attr:`_TIFF.BUFFERSIZE` * 2.
            returnoffset:
                Return offset and number of bytes of memory-mappable image
                data in file.

        Returns:
            If `returnoffset` is *True* and the image data in the file are
            memory-mappable, return the offset and number of bytes of the
            image data in the file.

        """
        # TODO: refactor this function

        fh: FileHandle
        storedshape: StoredShape = StoredShape(frames=-1)
        byteorder: Literal['>', '<']
        inputshape: tuple[int, ...]
        datashape: tuple[int, ...]
        dataarray: NDArray[Any] | None = None
        dataiter: Iterator[NDArray[Any] | bytes | None] | None = None
        dataoffsetsoffset: tuple[int, int | None] | None = None
        databytecountsoffset: tuple[int, int | None] | None = None
        subifdsoffsets: tuple[int, int | None] | None = None
        datadtype: numpy.dtype[Any]
        bilevel: bool
        tiles: tuple[int, ...]
        ifdpos: int
        photometricsamples: int
        pos: int | None = None
        predictortag: int
        predictorfunc: Callable[..., Any] | None = None
        compressiontag: int
        compressionfunc: Callable[..., Any] | None = None
        tags: list[tuple[int, bytes, bytes | None, bool]]
        numtiles: int
        numstrips: int

        fh = self._fh
        byteorder = self.tiff.byteorder

        if data is None:
            # empty
            if shape is None or dtype is None:
                raise ValueError(
                    "missing required 'shape' or 'dtype' arguments"
                )
            dataarray = None
            dataiter = None
            datashape = tuple(shape)
            datadtype = numpy.dtype(dtype).newbyteorder(byteorder)

        elif hasattr(data, '__next__'):
            # iterator/generator
            if shape is None or dtype is None:
                raise ValueError(
                    "missing required 'shape' or 'dtype' arguments"
                )
            dataiter = data  # type: ignore[assignment]
            datashape = tuple(shape)
            datadtype = numpy.dtype(dtype).newbyteorder(byteorder)

        elif hasattr(data, 'dtype'):
            # numpy, zarr, or dask array
            data = cast(numpy.ndarray, data)  # type: ignore[type-arg]
            dataarray = data
            datadtype = numpy.dtype(data.dtype).newbyteorder(byteorder)
            if not hasattr(data, 'reshape'):
                # zarr array cannot be shape-normalized
                dataarray = numpy.asarray(data, datadtype, 'C')
            else:
                try:
                    # numpy array must be C contiguous
                    if data.flags.f_contiguous:
                        dataarray = numpy.asarray(data, datadtype, 'C')
                except AttributeError:
                    # not a numpy array
                    pass
            datashape = dataarray.shape
            dataiter = None
            if dtype is not None and numpy.dtype(dtype) != datadtype:
                raise ValueError(
                    f'dtype argument {dtype!r} does not match '
                    f'data dtype {datadtype}'
                )
            if shape is not None and shape != dataarray.shape:
                raise ValueError(
                    f'shape argument {shape!r} does not match '
                    f'data shape {dataarray.shape}'
                )

        else:
            # scalar, list, tuple, etc
            # if dtype is not specified, default to float64
            datadtype = numpy.dtype(dtype).newbyteorder(byteorder)
            dataarray = numpy.asarray(data, datadtype, 'C')
            datashape = dataarray.shape
            dataiter = None

        del data

        if any(size >= 4294967296 for size in datashape):
            raise ValueError('invalid data shape')

        bilevel = datadtype.char == '?'
        if bilevel:
            index = -1 if datashape[-1] > 1 else -2
            datasize = product(datashape[:index])
            if datashape[index] % 8:
                datasize *= datashape[index] // 8 + 1
            else:
                datasize *= datashape[index] // 8
        else:
            datasize = product(datashape) * datadtype.itemsize

        if datasize == 0:
            dataarray = None
            compression = False
            bitspersample = None
            if metadata is not None:
                truncate = True

        if (
            not compression
            or (
                not isinstance(compression, bool)  # because True == 1
                and compression in ('NONE', 'None', 'none', 1)
            )
            or (
                isinstance(compression, (tuple, list))
                and compression[0] in (None, 0, 1, 'NONE', 'None', 'none')
            )
        ):
            compression = False

        if not predictor or (
            not isinstance(predictor, bool)  # because True == 1
            and predictor in {'NONE', 'None', 'none', 1}
        ):
            predictor = False

        inputshape = datashape

        packints = (
            bitspersample is not None
            and bitspersample != datadtype.itemsize * 8
        )

        # just append contiguous data if possible
        if self._datashape is not None and self._datadtype is not None:
            if colormap is not None:
                colormap = numpy.asarray(colormap, dtype=byteorder + 'H')
            if (
                not contiguous
                or self._datashape[1:] != datashape
                or self._datadtype != datadtype
                or (colormap is None and self._colormap is not None)
                or (self._colormap is None and colormap is not None)
                or not numpy.array_equal(
                    colormap, self._colormap  # type: ignore[arg-type]
                )
            ):
                # incompatible shape, dtype, or colormap
                self._write_remaining_pages()

                if self._imagej:
                    raise ValueError(
                        'the ImageJ format does not support '
                        'non-contiguous series'
                    )
                if self._omexml is not None:
                    if self._subifdslevel < 0:
                        # add image to OME-XML
                        assert self._storedshape is not None
                        assert self._metadata is not None
                        self._omexml.addimage(
                            dtype=self._datadtype,
                            shape=self._datashape[
                                0 if self._datashape[0] != 1 else 1 :
                            ],
                            storedshape=self._storedshape.shape,
                            **self._metadata,
                        )
                elif metadata is not None:
                    self._write_image_description()
                    # description might have been appended to file
                    fh.seek(0, os.SEEK_END)

                if self._subifds:
                    if self._truncate or truncate:
                        raise ValueError(
                            'SubIFDs cannot be used with truncated series'
                        )
                    self._subifdslevel += 1
                    if self._subifdslevel == self._subifds:
                        # done with writing SubIFDs
                        self._nextifdoffsets = []
                        self._subifdsoffsets = []
                        self._subifdslevel = -1
                        self._subifds = 0
                        self._ifdindex = 0
                    elif subifds:
                        raise ValueError(
                            'SubIFDs in SubIFDs are not supported'
                        )

                self._datashape = None
                self._colormap = None

            elif compression or packints or tile:
                raise ValueError(
                    'contiguous mode cannot be used with compression or tiles'
                )

            else:
                # consecutive mode
                # write all data, write IFDs/tags later
                self._datashape = (self._datashape[0] + 1,) + datashape
                offset = fh.tell()
                if dataarray is None:
                    fh.write_empty(datasize)
                else:
                    fh.write_array(dataarray, datadtype)
                if returnoffset:
                    return offset, datasize
                return None

        if self._ome is None:
            if description is None:
                self._ome = '.ome.' in fh.extension
            else:
                self._ome = False

        if self._tifffile or self._imagej:
            self._truncate = bool(truncate)
        elif truncate:
            raise ValueError(
                'truncate can only be used with imagej or shaped formats'
            )
        else:
            self._truncate = False

        if self._truncate and (compression or packints or tile):
            raise ValueError(
                'truncate cannot be used with compression, packints, or tiles'
            )

        if datasize == 0:
            # write single placeholder TiffPage for arrays with size=0
            datashape = (0, 0)
            warnings.warn(
                f'{self!r} writing zero-size array to nonconformant TIFF',
                UserWarning,
            )
            # TODO: reconsider this
            # raise ValueError('cannot save zero size array')

        tagnoformat = self.tiff.tagnoformat
        offsetformat = self.tiff.offsetformat
        offsetsize = self.tiff.offsetsize
        tagsize = self.tiff.tagsize

        MINISBLACK = PHOTOMETRIC.MINISBLACK
        MINISWHITE = PHOTOMETRIC.MINISWHITE
        RGB = PHOTOMETRIC.RGB
        YCBCR = PHOTOMETRIC.YCBCR
        PALETTE = PHOTOMETRIC.PALETTE
        CONTIG = PLANARCONFIG.CONTIG
        SEPARATE = PLANARCONFIG.SEPARATE

        # parse input
        if photometric is not None:
            photometric = enumarg(PHOTOMETRIC, photometric)
        if planarconfig:
            planarconfig = enumarg(PLANARCONFIG, planarconfig)
        if extrasamples is not None:
            # TODO: deprecate non-sequence extrasamples
            extrasamples = tuple(
                int(enumarg(EXTRASAMPLE, x)) for x in sequence(extrasamples)
            )

        if compressionargs is None:
            compressionargs = {}

        if compression:
            if isinstance(compression, (tuple, list)):
                # TODO: unreachable
                raise TypeError(
                    "passing multiple values to the 'compression' "
                    "parameter was deprecated in 2022.7.28. "
                    "Use 'compressionargs' to pass extra arguments to the "
                    "compression codec.",
                )
            if isinstance(compression, str):
                compression = compression.upper()
                if compression == 'ZLIB':
                    compression = 8  # ADOBE_DEFLATE
            elif isinstance(compression, bool):
                compression = 8  # ADOBE_DEFLATE
            compressiontag = enumarg(COMPRESSION, compression).value
            compression = True
        else:
            compressiontag = 1
            compression = False

        if compressiontag == 1:
            compressionargs = {}
        elif compressiontag in {33003, 33004, 33005, 34712}:
            # JPEG2000: use J2K instead of JP2
            compressionargs['codecformat'] = 0  # OPJ_CODEC_J2K

        assert compressionargs is not None

        if predictor:
            if not compression:
                raise ValueError('cannot use predictor without compression')
            if compressiontag in TIFF.IMAGE_COMPRESSIONS:
                # don't use predictor with JPEG, JPEG2000, WEBP, PNG, ...
                raise ValueError(
                    'cannot use predictor with '
                    f'{COMPRESSION(compressiontag)!r}'
                )
            if isinstance(predictor, bool):
                if datadtype.kind == 'f':
                    predictortag = 3
                elif datadtype.kind in 'iu' and datadtype.itemsize <= 4:
                    predictortag = 2
                else:
                    raise ValueError(
                        f'cannot use predictor with {datadtype!r}'
                    )
            else:
                predictor = enumarg(PREDICTOR, predictor)
                if (
                    datadtype.kind in 'iu'
                    and predictor.value not in {2, 34892, 34893}
                    and datadtype.itemsize <= 4
                ) or (
                    datadtype.kind == 'f'
                    and predictor.value not in {3, 34894, 34895}
                ):
                    raise ValueError(
                        f'cannot use {predictor!r} with {datadtype!r}'
                    )
                predictortag = predictor.value
        else:
            predictortag = 1

        del predictor
        predictorfunc = TIFF.PREDICTORS[predictortag]

        if self._ome:
            if description is not None:
                warnings.warn(
                    f'{self!r} not writing description to OME-TIFF',
                    UserWarning,
                )
                description = None
            if self._omexml is None:
                if metadata is None:
                    self._omexml = OmeXml()
                else:
                    self._omexml = OmeXml(**metadata)
            if volumetric or (tile and len(tile) > 2):
                raise ValueError('OME-TIFF does not support ImageDepth')
            volumetric = False

        elif self._imagej:
            # if tile is not None or predictor or compression:
            #     warnings.warn(
            #         f'{self!r} the ImageJ format does not support '
            #         'tiles, predictors, compression'
            #     )
            if description is not None:
                warnings.warn(
                    f'{self!r} not writing description to ImageJ file',
                    UserWarning,
                )
                description = None
            if datadtype.char not in 'BHhf':
                raise ValueError(
                    'the ImageJ format does not support data type '
                    f'{datadtype.char!r}'
                )
            if volumetric or (tile and len(tile) > 2):
                raise ValueError(
                    'the ImageJ format does not support ImageDepth'
                )
            volumetric = False
            ijrgb = photometric == RGB if photometric else None
            if datadtype.char != 'B':
                if photometric == RGB:
                    raise ValueError(
                        'the ImageJ format does not support '
                        f'data type {datadtype!r} for RGB'
                    )
                ijrgb = False
            if colormap is not None:
                ijrgb = False
            if metadata is None:
                axes = None
            else:
                axes = metadata.get('axes', None)
            ijshape = imagej_shape(datashape, rgb=ijrgb, axes=axes)
            if planarconfig == SEPARATE:
                raise ValueError(
                    'the ImageJ format does not support planar samples'
                )
            if ijshape[-1] in {3, 4}:
                photometric = RGB
            elif photometric is None:
                if colormap is not None and datadtype.char == 'B':
                    photometric = PALETTE
                else:
                    photometric = MINISBLACK
                planarconfig = None
            planarconfig = CONTIG if ijrgb else None

        # verify colormap and indices
        if colormap is not None:
            colormap = numpy.asarray(colormap, dtype=byteorder + 'H')
            self._colormap = colormap
            if self._imagej:
                if colormap.shape != (3, 256):
                    raise ValueError('invalid colormap shape for ImageJ')
                if datadtype.char == 'B' and photometric in {
                    MINISBLACK,
                    MINISWHITE,
                }:
                    photometric = PALETTE
                elif not (
                    (datadtype.char == 'B' and photometric == PALETTE)
                    or (
                        datadtype.char in 'Hf'
                        and photometric in {MINISBLACK, MINISWHITE}
                    )
                ):
                    warnings.warn(
                        f'{self!r} not writing colormap to ImageJ image with '
                        f'dtype={datadtype} and {photometric=}',
                        UserWarning,
                    )
                    colormap = None
            elif photometric is None and datadtype.char in 'BH':
                photometric = PALETTE
                planarconfig = None
                if colormap.shape != (3, 2 ** (datadtype.itemsize * 8)):
                    raise ValueError('invalid colormap shape')
            elif photometric == PALETTE:
                planarconfig = None
                if datadtype.char not in 'BH':
                    raise ValueError('invalid data dtype for palette-image')
                if colormap.shape != (3, 2 ** (datadtype.itemsize * 8)):
                    raise ValueError('invalid colormap shape')
            else:
                warnings.warn(
                    f'{self!r} not writing colormap with image of '
                    f'dtype={datadtype} and {photometric=}',
                    UserWarning,
                )
                colormap = None

        if tile:
            # verify tile shape

            if (
                not 1 < len(tile) < 4
                or tile[-1] % 16
                or tile[-2] % 16
                or any(i < 1 for i in tile)
            ):
                raise ValueError(f'invalid tile shape {tile}')
            tile = tuple(int(i) for i in tile)
            if volumetric and len(tile) == 2:
                tile = (1,) + tile
            volumetric = len(tile) == 3
        else:
            tile = ()
            volumetric = bool(volumetric)
        assert isinstance(tile, tuple)  # for mypy

        # normalize data shape to 5D or 6D, depending on volume:
        #   (pages, separate_samples, [depth,] length, width, contig_samples)
        shape = reshape_nd(
            datashape,
            TIFF.PHOTOMETRIC_SAMPLES.get(
                photometric, 2  # type: ignore[arg-type]
            ),
        )
        ndim = len(shape)

        if volumetric and ndim < 3:
            volumetric = False

        if photometric is None:
            deprecate = False
            photometric = MINISBLACK
            if bilevel:
                photometric = MINISWHITE
            elif planarconfig == CONTIG:
                if ndim > 2 and shape[-1] in {3, 4}:
                    photometric = RGB
                    deprecate = datadtype.char not in 'BH'
            elif planarconfig == SEPARATE:
                if volumetric and ndim > 3 and shape[-4] in {3, 4}:
                    photometric = RGB
                    deprecate = True
                elif ndim > 2 and shape[-3] in {3, 4}:
                    photometric = RGB
                    deprecate = True
            elif ndim > 2 and shape[-1] in {3, 4}:
                photometric = RGB
                planarconfig = CONTIG
                deprecate = datadtype.char not in 'BH'
            elif self._imagej or self._ome:
                photometric = MINISBLACK
                planarconfig = None
            elif volumetric and ndim > 3 and shape[-4] in {3, 4}:
                photometric = RGB
                planarconfig = SEPARATE
                deprecate = True
            elif ndim > 2 and shape[-3] in {3, 4}:
                photometric = RGB
                planarconfig = SEPARATE
                deprecate = True

            if deprecate:
                if planarconfig == CONTIG:
                    msg = 'contiguous samples', 'parameter is'
                else:
                    msg = (
                        'separate component planes',
                        "and 'planarconfig' parameters are",
                    )
                warnings.warn(
                    f"<tifffile.TiffWriter.write> data with shape {datashape} "
                    f"and dtype '{datadtype}' are stored as RGB with {msg[0]}."
                    " Future versions will store such data as MINISBLACK in "
                    "separate pages by default, unless the 'photometric' "
                    f"{msg[1]} specified.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                del msg
            del deprecate

        del datashape
        assert photometric is not None
        photometricsamples = TIFF.PHOTOMETRIC_SAMPLES[photometric]

        if planarconfig and len(shape) <= (3 if volumetric else 2):
            # TODO: raise error?
            planarconfig = None
            if photometricsamples > 1:
                photometric = MINISBLACK

        if photometricsamples > 1:
            if len(shape) < 3:
                raise ValueError(f'not a {photometric!r} image')
            if len(shape) < 4:
                volumetric = False
            if planarconfig is None:
                if photometric == RGB:
                    samples_set = {photometricsamples, 4}  # allow common alpha
                else:
                    samples_set = {photometricsamples}
                if shape[-1] in samples_set:
                    planarconfig = CONTIG
                elif shape[-4 if volumetric else -3] in samples_set:
                    planarconfig = SEPARATE
                elif shape[-1] > shape[-4 if volumetric else -3]:
                    # TODO: deprecated this?
                    planarconfig = SEPARATE
                else:
                    planarconfig = CONTIG
            if planarconfig == CONTIG:
                storedshape.contig_samples = shape[-1]
                storedshape.width = shape[-2]
                storedshape.length = shape[-3]
                if volumetric:
                    storedshape.depth = shape[-4]
            else:
                storedshape.width = shape[-1]
                storedshape.length = shape[-2]
                if volumetric:
                    storedshape.depth = shape[-3]
                    storedshape.separate_samples = shape[-4]
                else:
                    storedshape.separate_samples = shape[-3]
            if storedshape.samples > photometricsamples:
                storedshape.extrasamples = (
                    storedshape.samples - photometricsamples
                )

        elif photometric == PHOTOMETRIC.CFA:
            if len(shape) != 2:
                raise ValueError('invalid CFA image')
            volumetric = False
            planarconfig = None
            storedshape.width = shape[-1]
            storedshape.length = shape[-2]
            # if all(et[0] != 50706 for et in extratags):
            #     raise ValueError('must specify DNG tags for CFA image')

        elif planarconfig and len(shape) > (3 if volumetric else 2):
            if planarconfig == CONTIG:
                if extrasamples is None or len(extrasamples) > 0:
                    # use extrasamples
                    storedshape.contig_samples = shape[-1]
                    storedshape.width = shape[-2]
                    storedshape.length = shape[-3]
                    if volumetric:
                        storedshape.depth = shape[-4]
                else:
                    planarconfig = None
                    storedshape.contig_samples = 1
                    storedshape.width = shape[-1]
                    storedshape.length = shape[-2]
                    if volumetric:
                        storedshape.depth = shape[-3]
            else:
                storedshape.width = shape[-1]
                storedshape.length = shape[-2]
                if extrasamples is None or len(extrasamples) > 0:
                    # use extrasamples
                    if volumetric:
                        storedshape.depth = shape[-3]
                        storedshape.separate_samples = shape[-4]
                    else:
                        storedshape.separate_samples = shape[-3]
                else:
                    planarconfig = None
                    storedshape.separate_samples = 1
                    if volumetric:
                        storedshape.depth = shape[-3]
            storedshape.extrasamples = storedshape.samples - 1

        else:
            # photometricsamples == 1
            planarconfig = None
            if self._tifffile and (metadata or metadata == {}):
                # remove trailing 1s in shaped series
                while len(shape) > 2 and shape[-1] == 1:
                    shape = shape[:-1]
            elif self._imagej and len(shape) > 2 and shape[-1] == 1:
                # TODO: remove this and sync with ImageJ shape
                shape = shape[:-1]
            if len(shape) < 3:
                volumetric = False
            if not extrasamples:
                storedshape.width = shape[-1]
                storedshape.length = shape[-2]
                if volumetric:
                    storedshape.depth = shape[-3]
            else:
                storedshape.contig_samples = shape[-1]
                storedshape.width = shape[-2]
                storedshape.length = shape[-3]
                if volumetric:
                    storedshape.depth = shape[-4]
                storedshape.extrasamples = storedshape.samples - 1

        if not volumetric and tile and len(tile) == 3 and tile[0] > 1:
            raise ValueError(
                f'<tifffile.TiffWriter.write> cannot write {storedshape!r} '
                f'using volumetric tiles {tile}'
            )

        if subfiletype is not None and subfiletype & 0b100:
            # FILETYPE_MASK
            if not (
                bilevel
                and storedshape.samples == 1
                and photometric in {0, 1, 4}
            ):
                raise ValueError('invalid SubfileType MASK')
            photometric = PHOTOMETRIC.MASK

        packints = False
        if bilevel:
            if bitspersample is not None and bitspersample != 1:
                raise ValueError(f'{bitspersample=} must be 1 for bilevel')
            bitspersample = 1
        elif compressiontag in {6, 7, 34892, 33007}:
            # JPEG
            # TODO: add bitspersample to compressionargs?
            if bitspersample is None:
                if 'bitspersample' in compressionargs:
                    bitspersample = compressionargs['bitspersample']
                else:
                    bitspersample = 12 if datadtype == 'uint16' else 8
            if not 2 <= bitspersample <= 16:
                raise ValueError(
                    f'{bitspersample=} invalid for JPEG compression'
                )
        elif compressiontag in {33003, 33004, 33005, 34712, 50002, 52546}:
            # JPEG2K, JPEGXL
            # TODO: unify with JPEG?
            if bitspersample is None:
                if 'bitspersample' in compressionargs:
                    bitspersample = compressionargs['bitspersample']
                else:
                    bitspersample = datadtype.itemsize * 8
            if not (
                bitspersample > {1: 0, 2: 8, 4: 16}[datadtype.itemsize]
                and bitspersample <= datadtype.itemsize * 8
            ):
                raise ValueError(
                    f'{bitspersample=} out of range of {datadtype=}'
                )
        elif bitspersample is None:
            bitspersample = datadtype.itemsize * 8
        elif (
            datadtype.kind != 'u' or datadtype.itemsize > 4
        ) and bitspersample != datadtype.itemsize * 8:
            raise ValueError(f'{bitspersample=} does not match {datadtype=}')
        elif not (
            bitspersample > {1: 0, 2: 8, 4: 16}[datadtype.itemsize]
            and bitspersample <= datadtype.itemsize * 8
        ):
            raise ValueError(f'{bitspersample=} out of range of {datadtype=}')
        elif compression:
            if bitspersample != datadtype.itemsize * 8:
                raise ValueError(
                    f'{bitspersample=} cannot be used with compression'
                )
        elif bitspersample != datadtype.itemsize * 8:
            packints = True

        if storedshape.frames == -1:
            s0 = storedshape.page_size
            storedshape.frames = 1 if s0 == 0 else product(inputshape) // s0

        if datasize > 0 and not storedshape.is_valid:
            raise RuntimeError(f'invalid {storedshape!r}')

        if photometric == PALETTE:
            if storedshape.samples != 1 or storedshape.extrasamples > 0:
                raise ValueError(f'invalid {storedshape!r} for palette mode')
        elif storedshape.samples < photometricsamples:
            raise ValueError(
                f'not enough samples for {photometric!r}: '
                f'expected {photometricsamples}, got {storedshape.samples}'
            )

        if (
            planarconfig is not None
            and storedshape.planarconfig != planarconfig
        ):
            raise ValueError(
                f'{planarconfig!r} does not match {storedshape!r}'
            )
        del planarconfig

        if dataarray is not None:
            dataarray = dataarray.reshape(storedshape.shape)

        tags = []  # list of (code, ifdentry, ifdvalue, writeonce)

        if tile:
            tagbytecounts = 325  # TileByteCounts
            tagoffsets = 324  # TileOffsets
        else:
            tagbytecounts = 279  # StripByteCounts
            tagoffsets = 273  # StripOffsets
        self._dataoffsetstag = tagoffsets

        pack = self._pack
        addtag = self._addtag

        if extratags is None:
            extratags = ()

        if description is not None:
            # ImageDescription: user provided description
            addtag(tags, 270, 2, 0, description, True)

        # write shape and metadata to ImageDescription
        self._metadata = {} if not metadata else metadata.copy()
        if self._omexml is not None:
            if len(self._omexml.images) == 0:
                # rewritten later at end of file
                description = '\x00\x00\x00\x00'
            else:
                description = None
        elif self._imagej:
            ijmetadata = parse_kwargs(
                self._metadata,
                'Info',
                'Labels',
                'Ranges',
                'LUTs',
                'Plot',
                'ROI',
                'Overlays',
                'Properties',
                'info',
                'labels',
                'ranges',
                'luts',
                'plot',
                'roi',
                'overlays',
                'prop',
            )

            for t in imagej_metadata_tag(ijmetadata, byteorder):
                addtag(tags, *t)
            description = imagej_description(
                inputshape,
                rgb=storedshape.contig_samples in {3, 4},
                colormaped=self._colormap is not None,
                **self._metadata,
            )
            description += '\x00' * 64  # add buffer for in-place update
        elif self._tifffile and (metadata or metadata == {}):
            if self._truncate:
                self._metadata.update(truncated=True)
            description = shaped_description(inputshape, **self._metadata)
            description += '\x00' * 16  # add buffer for in-place update
        # elif metadata is None and self._truncate:
        #     raise ValueError('cannot truncate without writing metadata')
        elif description is not None:
            if not isinstance(description, bytes):
                description = description.encode('ascii')
            self._descriptiontag = TiffTag(
                self, 0, 270, 2, len(description), description, 0
            )
            description = None

        if description is None:
            # disable shaped format if user disabled metadata
            self._tifffile = False
        else:
            description = description.encode('ascii')
            addtag(tags, 270, 2, 0, description, True)
            self._descriptiontag = TiffTag(
                self, 0, 270, 2, len(description), description, 0
            )
        del description

        if software is None:
            software = 'tifffile.py'
        if software:
            addtag(tags, 305, 2, 0, software, True)
        if datetime:
            if isinstance(datetime, str):
                if len(datetime) != 19 or datetime[16] != ':':
                    raise ValueError('invalid datetime string')
            elif isinstance(datetime, DateTime):
                datetime = datetime.strftime('%Y:%m:%d %H:%M:%S')
            else:
                datetime = DateTime.now().strftime('%Y:%m:%d %H:%M:%S')
            addtag(tags, 306, 2, 0, datetime, True)
        addtag(tags, 259, 3, 1, compressiontag)  # Compression
        if compressiontag == 34887:
            # LERC
            if 'compression' not in compressionargs:
                lerc_compression = 0
            elif compressionargs['compression'] is None:
                lerc_compression = 0
            elif compressionargs['compression'] == 'deflate':
                lerc_compression = 1
            elif compressionargs['compression'] == 'zstd':
                lerc_compression = 2
            else:
                raise ValueError(
                    'invalid LERC compression '
                    f'{compressionargs["compression"]!r}'
                )
            addtag(tags, 50674, 4, 2, (4, lerc_compression))
            del lerc_compression
        if predictortag != 1:
            addtag(tags, 317, 3, 1, predictortag)
        addtag(tags, 256, 4, 1, storedshape.width)  # ImageWidth
        addtag(tags, 257, 4, 1, storedshape.length)  # ImageLength
        if tile:
            addtag(tags, 322, 4, 1, tile[-1])  # TileWidth
            addtag(tags, 323, 4, 1, tile[-2])  # TileLength
        if volumetric:
            addtag(tags, 32997, 4, 1, storedshape.depth)  # ImageDepth
            if tile:
                addtag(tags, 32998, 4, 1, tile[0])  # TileDepth
        if subfiletype is not None:
            addtag(tags, 254, 4, 1, subfiletype)  # NewSubfileType
        if (subifds or self._subifds) and self._subifdslevel < 0:
            if self._subifds:
                subifds = self._subifds
            elif hasattr(subifds, '__len__'):
                # allow TiffPage.subifds tuple
                subifds = len(subifds)  # type: ignore[arg-type]
            else:
                subifds = int(subifds)  # type: ignore[arg-type]
            self._subifds = subifds
            addtag(
                tags, 330, 18 if offsetsize > 4 else 13, subifds, [0] * subifds
            )
        if not bilevel and not datadtype.kind == 'u':
            # SampleFormat
            sampleformat = {'u': 1, 'i': 2, 'f': 3, 'c': 6}[datadtype.kind]
            addtag(
                tags,
                339,
                3,
                storedshape.samples,
                (sampleformat,) * storedshape.samples,
            )
        if colormap is not None:
            addtag(tags, 320, 3, colormap.size, colormap)
        if iccprofile is not None:
            addtag(tags, 34675, 7, len(iccprofile), iccprofile)
        addtag(tags, 277, 3, 1, storedshape.samples)
        if bilevel:
            # PlanarConfiguration
            if storedshape.samples > 1:
                addtag(tags, 284, 3, 1, storedshape.planarconfig)
        elif storedshape.samples > 1:
            # PlanarConfiguration
            addtag(tags, 284, 3, 1, storedshape.planarconfig)
            # BitsPerSample
            addtag(
                tags,
                258,
                3,
                storedshape.samples,
                (bitspersample,) * storedshape.samples,
            )
        else:
            addtag(tags, 258, 3, 1, bitspersample)
        if storedshape.extrasamples > 0:
            if extrasamples is not None:
                if storedshape.extrasamples != len(extrasamples):
                    raise ValueError(
                        'wrong number of extrasamples '
                        f'{storedshape.extrasamples} != {len(extrasamples)}'
                    )
                addtag(tags, 338, 3, len(extrasamples), extrasamples)
            elif photometric == RGB and storedshape.extrasamples == 1:
                # Unassociated alpha channel
                addtag(tags, 338, 3, 1, 2)
            else:
                # Unspecified alpha channel
                addtag(
                    tags,
                    338,
                    3,
                    storedshape.extrasamples,
                    (0,) * storedshape.extrasamples,
                )

        if jpegtables is not None:
            addtag(tags, 347, 7, len(jpegtables), jpegtables)

        if (
            compressiontag == 7
            and storedshape.planarconfig == 1
            and photometric in {RGB, YCBCR}
        ):
            # JPEG compression with subsampling
            # TODO: use JPEGTables for multiple tiles or strips
            if subsampling is None:
                subsampling = (2, 2)
            elif subsampling not in {(1, 1), (2, 1), (2, 2), (4, 1)}:
                raise ValueError(
                    f'invalid subsampling factors {subsampling!r}'
                )
            maxsampling = max(subsampling) * 8
            if tile and (tile[-1] % maxsampling or tile[-2] % maxsampling):
                raise ValueError(f'tile shape not a multiple of {maxsampling}')
            if storedshape.extrasamples > 1:
                raise ValueError('JPEG subsampling requires RGB(A) images')
            addtag(tags, 530, 3, 2, subsampling)  # YCbCrSubSampling
            # use PhotometricInterpretation YCBCR by default
            outcolorspace = enumarg(
                PHOTOMETRIC, compressionargs.get('outcolorspace', 6)
            )
            compressionargs['subsampling'] = subsampling
            compressionargs['colorspace'] = photometric.name
            compressionargs['outcolorspace'] = outcolorspace.name
            addtag(tags, 262, 3, 1, outcolorspace)
            if outcolorspace == YCBCR:
                # ReferenceBlackWhite is required for YCBCR
                if all(et[0] != 532 for et in extratags):
                    addtag(
                        tags,
                        532,
                        5,
                        6,
                        (0, 1, 255, 1, 128, 1, 255, 1, 128, 1, 255, 1),
                    )
        else:
            if subsampling not in {None, (1, 1)}:
                logger().warning(
                    f'{self!r} cannot apply subsampling {subsampling!r}'
                )
            subsampling = None
            maxsampling = 1
            addtag(
                tags, 262, 3, 1, photometric.value
            )  # PhotometricInterpretation
            if photometric == YCBCR:
                # YCbCrSubSampling and ReferenceBlackWhite
                addtag(tags, 530, 3, 2, (1, 1))
                if all(et[0] != 532 for et in extratags):
                    addtag(
                        tags,
                        532,
                        5,
                        6,
                        (0, 1, 255, 1, 128, 1, 255, 1, 128, 1, 255, 1),
                    )

        if resolutionunit is not None:
            resolutionunit = enumarg(RESUNIT, resolutionunit)
        elif self._imagej or resolution is None:
            resolutionunit = RESUNIT.NONE
        else:
            resolutionunit = RESUNIT.INCH

        if resolution is not None:
            addtag(tags, 282, 5, 1, rational(resolution[0]))  # XResolution
            addtag(tags, 283, 5, 1, rational(resolution[1]))  # YResolution
            if len(resolution) > 2:
                # TODO: unreachable
                raise ValueError(
                    "passing a unit along with the 'resolution' parameter "
                    "was deprecated in 2022.7.28. "
                    "Use the 'resolutionunit' parameter.",
                )
            addtag(tags, 296, 3, 1, resolutionunit)  # ResolutionUnit
        else:
            addtag(tags, 282, 5, 1, (1, 1))  # XResolution
            addtag(tags, 283, 5, 1, (1, 1))  # YResolution
            addtag(tags, 296, 3, 1, resolutionunit)  # ResolutionUnit

        # can save data array contiguous
        contiguous = not (compression or packints or bilevel)
        if tile:
            # one chunk per tile per plane
            if len(tile) == 2:
                tiles = (
                    (storedshape.length + tile[0] - 1) // tile[0],
                    (storedshape.width + tile[1] - 1) // tile[1],
                )
                contiguous = (
                    contiguous
                    and storedshape.length == tile[0]
                    and storedshape.width == tile[1]
                )
            else:
                tiles = (
                    (storedshape.depth + tile[0] - 1) // tile[0],
                    (storedshape.length + tile[1] - 1) // tile[1],
                    (storedshape.width + tile[2] - 1) // tile[2],
                )
                contiguous = (
                    contiguous
                    and storedshape.depth == tile[0]
                    and storedshape.length == tile[1]
                    and storedshape.width == tile[2]
                )
            numtiles = product(tiles) * storedshape.separate_samples
            databytecounts = [
                product(tile) * storedshape.contig_samples * datadtype.itemsize
            ] * numtiles
            bytecountformat = self._bytecount_format(
                databytecounts, compressiontag
            )
            addtag(
                tags, tagbytecounts, bytecountformat, numtiles, databytecounts
            )
            addtag(tags, tagoffsets, offsetformat, numtiles, [0] * numtiles)
            bytecountformat = f'{numtiles}{bytecountformat}'
            if not contiguous:
                if dataarray is not None:
                    dataiter = iter_tiles(dataarray, tile, tiles)
                elif dataiter is None and not (
                    compression or packints or bilevel
                ):

                    def dataiter_(
                        numtiles: int = numtiles * storedshape.frames,
                        bytecount: int = databytecounts[0],
                    ) -> Iterator[bytes]:
                        # yield empty tiles
                        chunk = bytes(bytecount)
                        for _ in range(numtiles):
                            yield chunk

                    dataiter = dataiter_()

            rowsperstrip = 0

        elif contiguous and (
            rowsperstrip is None or rowsperstrip >= storedshape.length
        ):
            count = storedshape.separate_samples * storedshape.depth
            databytecounts = [
                storedshape.length
                * storedshape.width
                * storedshape.contig_samples
                * datadtype.itemsize
            ] * count
            bytecountformat = self._bytecount_format(
                databytecounts, compressiontag
            )
            addtag(tags, tagbytecounts, bytecountformat, count, databytecounts)
            addtag(tags, tagoffsets, offsetformat, count, [0] * count)
            addtag(tags, 278, 4, 1, storedshape.length)  # RowsPerStrip
            bytecountformat = f'{count}{bytecountformat}'
            rowsperstrip = storedshape.length
            numstrips = count

        else:
            # use rowsperstrip
            rowsize = (
                storedshape.width
                * storedshape.contig_samples
                * datadtype.itemsize
            )
            if compressiontag == 48124:
                # Jetraw works on whole camera frame
                rowsperstrip = storedshape.length
            if rowsperstrip is None:
                # compress ~256 KB chunks by default
                # TIFF-EP requires <= 64 KB
                if compression:
                    rowsperstrip = 262144 // rowsize
                else:
                    rowsperstrip = storedshape.length
            if rowsperstrip < 1:
                rowsperstrip = maxsampling
            elif rowsperstrip > storedshape.length:
                rowsperstrip = storedshape.length
            elif subsampling and rowsperstrip % maxsampling:
                rowsperstrip = (
                    math.ceil(rowsperstrip / maxsampling) * maxsampling
                )
            assert rowsperstrip is not None
            addtag(tags, 278, 4, 1, rowsperstrip)  # RowsPerStrip

            numstrips1 = (
                storedshape.length + rowsperstrip - 1
            ) // rowsperstrip
            numstrips = (
                numstrips1 * storedshape.separate_samples * storedshape.depth
            )
            # TODO: save bilevel data with rowsperstrip
            stripsize = rowsperstrip * rowsize
            databytecounts = [stripsize] * numstrips
            laststripsize = stripsize - rowsize * (
                numstrips1 * rowsperstrip - storedshape.length
            )
            for i in range(numstrips1 - 1, numstrips, numstrips1):
                databytecounts[i] = laststripsize
            bytecountformat = self._bytecount_format(
                databytecounts, compressiontag
            )
            addtag(
                tags, tagbytecounts, bytecountformat, numstrips, databytecounts
            )
            addtag(tags, tagoffsets, offsetformat, numstrips, [0] * numstrips)
            bytecountformat = bytecountformat * numstrips

            if dataarray is not None and not contiguous:
                dataiter = iter_images(dataarray)

        if dataiter is None and not contiguous:
            raise ValueError('cannot write non-contiguous empty file')

        # add extra tags from user; filter duplicate and select tags
        extratag: TagTuple
        tagset = {t[0] for t in tags}
        tagset.update(TIFF.TAG_FILTERED)
        for extratag in extratags:
            if extratag[0] in tagset:
                logger().warning(
                    f'{self!r} not writing extratag {extratag[0]}'
                )
            else:
                addtag(tags, *extratag)
        del tagset
        del extratags

        # TODO: check TIFFReadDirectoryCheckOrder warning in files containing
        #   multiple tags of same code
        # the entries in an IFD must be sorted in ascending order by tag code
        tags = sorted(tags, key=lambda x: x[0])

        # define compress function
        compressionaxis: int = -2
        bytesiter: bool = False

        iteritem: NDArray[Any] | bytes | None
        if dataiter is not None:
            iteritem, dataiter = peek_iterator(dataiter)
            bytesiter = isinstance(iteritem, bytes)
            if not bytesiter:
                iteritem = numpy.asarray(iteritem)
                if (
                    tile
                    and storedshape.contig_samples == 1
                    and iteritem.shape[-1] != 1
                ):
                    # issue 185
                    compressionaxis = -1
                if iteritem.dtype.char != datadtype.char:
                    raise ValueError(
                        f'dtype of iterator {iteritem.dtype!r} '
                        f'does not match dtype {datadtype!r}'
                    )
        else:
            iteritem = None

        if bilevel:
            if compressiontag == 1:

                def compressionfunc1(
                    data: Any, axis: int = compressionaxis
                ) -> bytes:
                    return numpy.packbits(data, axis=axis).tobytes()

                compressionfunc = compressionfunc1

            elif compressiontag in {5, 32773, 8, 32946, 50013, 34925, 50000}:
                # LZW, PackBits, deflate, LZMA, ZSTD
                def compressionfunc2(
                    data: Any,
                    compressor: Any = TIFF.COMPRESSORS[compressiontag],
                    axis: int = compressionaxis,
                    kwargs: Any = compressionargs,
                ) -> bytes:
                    data = numpy.packbits(data, axis=axis).tobytes()
                    return compressor(data, **kwargs)

                compressionfunc = compressionfunc2

            else:
                raise NotImplementedError('cannot compress bilevel image')

        elif compression:
            compressor = TIFF.COMPRESSORS[compressiontag]

            if compressiontag == 32773:
                # PackBits
                compressionargs['axis'] = compressionaxis

            # elif compressiontag == 48124:
            #     # Jetraw
            #     imagecodecs.jetraw_init(
            #         parameters=compressionargs.pop('parameters', None),
            #         verbose=compressionargs.pop('verbose', None),
            #     )
            #     if not 'identifier' in compressionargs:
            #         raise ValueError(
            #             "jetraw_encode() missing argument: 'identifier'"
            #         )

            if subsampling:
                # JPEG with subsampling
                def compressionfunc(
                    data: Any,
                    compressor: Any = compressor,
                    kwargs: Any = compressionargs,
                ) -> bytes:
                    return compressor(data, **kwargs)

            elif predictorfunc is not None:

                def compressionfunc(
                    data: Any,
                    predictorfunc: Any = predictorfunc,
                    compressor: Any = compressor,
                    axis: int = compressionaxis,
                    kwargs: Any = compressionargs,
                ) -> bytes:
                    data = predictorfunc(data, axis=axis)
                    return compressor(data, **kwargs)

            elif compressionargs:

                def compressionfunc(
                    data: Any,
                    compressor: Any = compressor,
                    kwargs: Any = compressionargs,
                ) -> bytes:
                    return compressor(data, **kwargs)

            elif compressiontag > 1:
                compressionfunc = compressor

            else:
                compressionfunc = None

        elif packints:

            def compressionfunc(
                data: Any,
                bps: Any = bitspersample,
                axis: int = compressionaxis,
            ) -> bytes:
                return imagecodecs.packints_encode(data, bps, axis=axis)

        else:
            compressionfunc = None

        del compression
        if not contiguous and not bytesiter and compressionfunc is not None:
            # create iterator of encoded tiles or strips
            bytesiter = True
            if tile:
                # dataiter yields tiles
                tileshape = tile + (storedshape.contig_samples,)
                tilesize = product(tileshape) * datadtype.itemsize
                maxworkers = TiffWriter._maxworkers(
                    maxworkers,
                    numtiles * storedshape.frames,
                    tilesize,
                    compressiontag,
                )
                # yield encoded tiles
                dataiter = encode_chunks(
                    numtiles * storedshape.frames,
                    dataiter,  # type: ignore[arg-type]
                    compressionfunc,
                    tileshape,
                    datadtype,
                    maxworkers,
                    buffersize,
                    True,
                )
            else:
                # dataiter yields frames
                maxworkers = TiffWriter._maxworkers(
                    maxworkers,
                    numstrips * storedshape.frames,
                    stripsize,
                    compressiontag,
                )
                # yield strips
                dataiter = iter_strips(
                    dataiter,  # type: ignore[arg-type]
                    storedshape.page_shape,
                    datadtype,
                    rowsperstrip,
                )
                # yield encoded strips
                dataiter = encode_chunks(
                    numstrips * storedshape.frames,
                    dataiter,
                    compressionfunc,
                    (
                        rowsperstrip,
                        storedshape.width,
                        storedshape.contig_samples,
                    ),
                    datadtype,
                    maxworkers,
                    buffersize,
                    False,
                )

        fhpos = fh.tell()
        # commented out to allow image data beyond 4GB in classic TIFF
        # if (
        #     not (
        #         offsetsize > 4
        #         or self._imagej or compressionfunc is not None
        #     )
        #     and fhpos + datasize > 2**32 - 1
        # ):
        #     raise ValueError('data too large for classic TIFF format')

        dataoffset: int = 0

        # if not compressed or multi-tiled, write the first IFD and then
        # all data contiguously; else, write all IFDs and data interleaved
        for pageindex in range(1 if contiguous else storedshape.frames):
            ifdpos = fhpos
            if ifdpos % 2:
                # position of IFD must begin on a word boundary
                fh.write(b'\x00')
                ifdpos += 1

            if self._subifdslevel < 0:
                # update pointer at ifdoffset
                fh.seek(self._ifdoffset)
                fh.write(pack(offsetformat, ifdpos))

            fh.seek(ifdpos)

            # create IFD in memory
            if pageindex < 2:
                subifdsoffsets = None
                ifd = io.BytesIO()
                ifd.write(pack(tagnoformat, len(tags)))
                tagoffset = ifd.tell()
                ifd.write(b''.join(t[1] for t in tags))
                ifdoffset = ifd.tell()
                ifd.write(pack(offsetformat, 0))  # offset to next IFD
                # write tag values and patch offsets in ifdentries
                for tagindex, tag in enumerate(tags):
                    offset = tagoffset + tagindex * tagsize + 4 + offsetsize
                    code = tag[0]
                    value = tag[2]
                    if value:
                        pos = ifd.tell()
                        if pos % 2:
                            # tag value is expected to begin on word boundary
                            ifd.write(b'\x00')
                            pos += 1
                        ifd.seek(offset)
                        ifd.write(pack(offsetformat, ifdpos + pos))
                        ifd.seek(pos)
                        ifd.write(value)
                        if code == tagoffsets:
                            dataoffsetsoffset = offset, pos
                        elif code == tagbytecounts:
                            databytecountsoffset = offset, pos
                        elif code == 270:
                            if (
                                self._descriptiontag is not None
                                and self._descriptiontag.offset == 0
                                and value.startswith(
                                    self._descriptiontag.value
                                )
                            ):
                                self._descriptiontag.offset = (
                                    ifdpos + tagoffset + tagindex * tagsize
                                )
                                self._descriptiontag.valueoffset = ifdpos + pos
                        elif code == 330:
                            subifdsoffsets = offset, pos
                    elif code == tagoffsets:
                        dataoffsetsoffset = offset, None
                    elif code == tagbytecounts:
                        databytecountsoffset = offset, None
                    elif code == 270:
                        if (
                            self._descriptiontag is not None
                            and self._descriptiontag.offset == 0
                            and self._descriptiontag.value in tag[1][-4:]
                        ):
                            self._descriptiontag.offset = (
                                ifdpos + tagoffset + tagindex * tagsize
                            )
                            self._descriptiontag.valueoffset = (
                                self._descriptiontag.offset + offsetsize + 4
                            )
                    elif code == 330:
                        subifdsoffsets = offset, None
                ifdsize = ifd.tell()
                if ifdsize % 2:
                    ifd.write(b'\x00')
                    ifdsize += 1

            # write IFD later when strip/tile bytecounts and offsets are known
            fh.seek(ifdsize, os.SEEK_CUR)

            # write image data
            dataoffset = fh.tell()
            if align is None:
                align = 16
            skip = (align - (dataoffset % align)) % align
            fh.seek(skip, os.SEEK_CUR)
            dataoffset += skip

            if contiguous:
                # write all image data contiguously
                if dataiter is not None:
                    byteswritten = 0
                    if bytesiter:
                        for iteritem in dataiter:
                            # assert isinstance(iteritem, bytes)
                            byteswritten += fh.write(
                                iteritem  # type: ignore[arg-type]
                            )
                            del iteritem
                    else:
                        pagesize = storedshape.page_size * datadtype.itemsize
                        for iteritem in dataiter:
                            if iteritem is None:
                                byteswritten += fh.write_empty(pagesize)
                            else:
                                # assert isinstance(iteritem, numpy.ndarray)
                                byteswritten += fh.write_array(
                                    iteritem,  # type: ignore[arg-type]
                                    datadtype,
                                )
                            del iteritem
                    if byteswritten != datasize:
                        raise ValueError(
                            'iterator contains wrong number of bytes '
                            f'{byteswritten} != {datasize}'
                        )
                elif dataarray is None:
                    fh.write_empty(datasize)
                else:
                    fh.write_array(dataarray, datadtype)

            elif bytesiter:
                # write tiles or strips
                assert dataiter is not None
                for chunkindex in range(numtiles if tile else numstrips):
                    iteritem = cast(bytes, next(dataiter))
                    # assert isinstance(iteritem, bytes)
                    databytecounts[chunkindex] = len(iteritem)
                    fh.write(iteritem)
                    del iteritem

            elif tile:
                # write uncompressed tiles
                assert dataiter is not None
                tileshape = tile + (storedshape.contig_samples,)
                tilesize = product(tileshape) * datadtype.itemsize
                for tileindex in range(numtiles):
                    iteritem = next(dataiter)
                    if iteritem is None:
                        databytecounts[tileindex] = 0
                        # fh.write_empty(tilesize)
                        continue
                    # assert not isinstance(iteritem, bytes)
                    iteritem = numpy.ascontiguousarray(iteritem, datadtype)
                    if iteritem.nbytes != tilesize:
                        # if iteritem.dtype != datadtype:
                        #     raise ValueError(
                        #         'dtype of tile does not match data'
                        #     )
                        if iteritem.nbytes > tilesize:
                            raise ValueError('tile is too large')
                        pad = tuple(
                            (0, i - j)
                            for i, j in zip(tileshape, iteritem.shape)
                        )
                        iteritem = numpy.pad(iteritem, pad)
                    fh.write_array(iteritem)
                    del iteritem

            else:
                raise RuntimeError('unreachable code')

            # update strip/tile offsets
            assert dataoffsetsoffset is not None
            offset, pos = dataoffsetsoffset
            ifd.seek(offset)
            if pos is not None:
                ifd.write(pack(offsetformat, ifdpos + pos))
                ifd.seek(pos)
                offset = dataoffset
                for size in databytecounts:
                    ifd.write(pack(offsetformat, offset if size > 0 else 0))
                    offset += size
            else:
                ifd.write(pack(offsetformat, dataoffset))

            if compressionfunc is not None or (tile and dataarray is None):
                # update strip/tile bytecounts
                assert databytecountsoffset is not None
                offset, pos = databytecountsoffset
                ifd.seek(offset)
                if pos is not None:
                    ifd.write(pack(offsetformat, ifdpos + pos))
                    ifd.seek(pos)
                ifd.write(pack(bytecountformat, *databytecounts))

            if subifdsoffsets is not None:
                # update and save pointer to SubIFDs tag values if necessary
                offset, pos = subifdsoffsets
                if pos is not None:
                    ifd.seek(offset)
                    ifd.write(pack(offsetformat, ifdpos + pos))
                    self._subifdsoffsets.append(ifdpos + pos)
                else:
                    self._subifdsoffsets.append(ifdpos + offset)

            fhpos = fh.tell()
            fh.seek(ifdpos)
            fh.write(ifd.getbuffer())
            fh.flush()

            if self._subifdslevel < 0:
                self._ifdoffset = ifdpos + ifdoffset
            else:
                # update SubIFDs tag values
                fh.seek(
                    self._subifdsoffsets[self._ifdindex]
                    + self._subifdslevel * offsetsize
                )
                fh.write(pack(offsetformat, ifdpos))

                # update SubIFD chain offsets
                if self._subifdslevel == 0:
                    self._nextifdoffsets.append(ifdpos + ifdoffset)
                else:
                    fh.seek(self._nextifdoffsets[self._ifdindex])
                    fh.write(pack(offsetformat, ifdpos))
                    self._nextifdoffsets[self._ifdindex] = ifdpos + ifdoffset
                self._ifdindex += 1
                self._ifdindex %= len(self._subifdsoffsets)

            fh.seek(fhpos)

            # remove tags that should be written only once
            if pageindex == 0:
                tags = [tag for tag in tags if not tag[-1]]

        assert dataoffset > 0

        self._datashape = (1,) + inputshape
        self._datadtype = datadtype
        self._dataoffset = dataoffset
        self._databytecounts = databytecounts
        self._storedshape = storedshape

        if contiguous:
            # write remaining IFDs/tags later
            self._tags = tags
            # return offset and size of image data
            if returnoffset:
                return dataoffset, sum(databytecounts)
        return None

    def overwrite_description(self, description: str, /) -> None:
        """Overwrite value of last ImageDescription tag.

        Can be used to write OME-XML after writing images.
        Ends a contiguous series.

        """
        if self._descriptiontag is None:
            raise ValueError('no ImageDescription tag found')
        self._write_remaining_pages()
        self._descriptiontag.overwrite(description, erase=False)
        self._descriptiontag = None

    def close(self) -> None:
        """Write remaining pages and close file handle."""
        try:
            if not self._truncate:
                self._write_remaining_pages()
            self._write_image_description()
        finally:
            try:
                self._fh.close()
            except Exception:
                pass

    @property
    def filehandle(self) -> FileHandle:
        """File handle to write file."""
        return self._fh

    def _write_remaining_pages(self) -> None:
        """Write outstanding IFDs and tags to file."""
        if not self._tags or self._truncate or self._datashape is None:
            return

        assert self._storedshape is not None
        assert self._databytecounts is not None
        assert self._dataoffset is not None

        pageno: int = self._storedshape.frames * self._datashape[0] - 1
        if pageno < 1:
            self._tags = None
            self._dataoffset = None
            self._databytecounts = None
            return

        fh = self._fh
        fhpos: int = fh.tell()
        if fhpos % 2:
            fh.write(b'\x00')
            fhpos += 1

        pack = struct.pack
        offsetformat: str = self.tiff.offsetformat
        offsetsize: int = self.tiff.offsetsize
        tagnoformat: str = self.tiff.tagnoformat
        tagsize: int = self.tiff.tagsize
        dataoffset: int = self._dataoffset
        pagedatasize: int = sum(self._databytecounts)
        subifdsoffsets: tuple[int, int | None] | None = None
        dataoffsetsoffset: tuple[int, int | None]
        pos: int | None
        offset: int

        # construct template IFD in memory
        # must patch offsets to next IFD and data before writing to file
        ifd = io.BytesIO()
        ifd.write(pack(tagnoformat, len(self._tags)))
        tagoffset = ifd.tell()
        ifd.write(b''.join(t[1] for t in self._tags))
        ifdoffset = ifd.tell()
        ifd.write(pack(offsetformat, 0))  # offset to next IFD
        # tag values
        for tagindex, tag in enumerate(self._tags):
            offset = tagoffset + tagindex * tagsize + offsetsize + 4
            code = tag[0]
            value = tag[2]
            if value:
                pos = ifd.tell()
                if pos % 2:
                    # tag value is expected to begin on word boundary
                    ifd.write(b'\x00')
                    pos += 1
                ifd.seek(offset)
                try:
                    ifd.write(pack(offsetformat, fhpos + pos))
                except Exception as exc:  # struct.error
                    if self._imagej:
                        warnings.warn(
                            f'{self!r} truncating ImageJ file', UserWarning
                        )
                        self._truncate = True
                        return
                    raise ValueError(
                        'data too large for non-BigTIFF file'
                    ) from exc
                ifd.seek(pos)
                ifd.write(value)
                if code == self._dataoffsetstag:
                    # save strip/tile offsets for later updates
                    dataoffsetsoffset = offset, pos
                elif code == 330:
                    # save subifds offsets for later updates
                    subifdsoffsets = offset, pos
            elif code == self._dataoffsetstag:
                dataoffsetsoffset = offset, None
            elif code == 330:
                subifdsoffsets = offset, None

        ifdsize = ifd.tell()
        if ifdsize % 2:
            ifd.write(b'\x00')
            ifdsize += 1

        # check if all IFDs fit in file
        if offsetsize < 8 and fhpos + ifdsize * pageno > 2**32 - 32:
            if self._imagej:
                warnings.warn(f'{self!r} truncating ImageJ file', UserWarning)
                self._truncate = True
                return
            raise ValueError('data too large for non-BigTIFF file')

        # assemble IFD chain in memory from IFD template
        ifds = io.BytesIO(bytes(ifdsize * pageno))
        ifdpos = fhpos
        for _ in range(pageno):
            # update strip/tile offsets in IFD
            dataoffset += pagedatasize  # offset to image data
            offset, pos = dataoffsetsoffset
            ifd.seek(offset)
            if pos is not None:
                ifd.write(pack(offsetformat, ifdpos + pos))
                ifd.seek(pos)
                offset = dataoffset
                for size in self._databytecounts:
                    ifd.write(pack(offsetformat, offset))
                    offset += size
            else:
                ifd.write(pack(offsetformat, dataoffset))

            if subifdsoffsets is not None:
                offset, pos = subifdsoffsets
                self._subifdsoffsets.append(
                    ifdpos + (pos if pos is not None else offset)
                )

            if self._subifdslevel < 0:
                if subifdsoffsets is not None:
                    # update pointer to SubIFDs tag values if necessary
                    offset, pos = subifdsoffsets
                    if pos is not None:
                        ifd.seek(offset)
                        ifd.write(pack(offsetformat, ifdpos + pos))

                # update pointer at ifdoffset to point to next IFD in file
                ifdpos += ifdsize
                ifd.seek(ifdoffset)
                ifd.write(pack(offsetformat, ifdpos))

            else:
                # update SubIFDs tag values in file
                fh.seek(
                    self._subifdsoffsets[self._ifdindex]
                    + self._subifdslevel * offsetsize
                )
                fh.write(pack(offsetformat, ifdpos))

                # update SubIFD chain
                if self._subifdslevel == 0:
                    self._nextifdoffsets.append(ifdpos + ifdoffset)
                else:
                    fh.seek(self._nextifdoffsets[self._ifdindex])
                    fh.write(pack(offsetformat, ifdpos))
                    self._nextifdoffsets[self._ifdindex] = ifdpos + ifdoffset
                self._ifdindex += 1
                self._ifdindex %= len(self._subifdsoffsets)
                ifdpos += ifdsize

            # write IFD entry
            ifds.write(ifd.getbuffer())

        # terminate IFD chain
        ifdoffset += ifdsize * (pageno - 1)
        ifds.seek(ifdoffset)
        ifds.write(pack(offsetformat, 0))
        # write IFD chain to file
        fh.seek(fhpos)
        fh.write(ifds.getbuffer())

        if self._subifdslevel < 0:
            # update file to point to new IFD chain
            pos = fh.tell()
            fh.seek(self._ifdoffset)
            fh.write(pack(offsetformat, fhpos))
            fh.flush()
            fh.seek(pos)
            self._ifdoffset = fhpos + ifdoffset

        self._tags = None
        self._dataoffset = None
        self._databytecounts = None
        # do not reset _storedshape, _datashape, _datadtype

    def _write_image_description(self) -> None:
        """Write metadata to ImageDescription tag."""
        if self._datashape is None or self._descriptiontag is None:
            self._descriptiontag = None
            return

        assert self._storedshape is not None
        assert self._datadtype is not None

        if self._omexml is not None:
            if self._subifdslevel < 0:
                assert self._metadata is not None
                self._omexml.addimage(
                    dtype=self._datadtype,
                    shape=self._datashape[
                        0 if self._datashape[0] != 1 else 1 :
                    ],
                    storedshape=self._storedshape.shape,
                    **self._metadata,
                )
            description = self._omexml.tostring(declaration=True)
        elif self._datashape[0] == 1:
            # description already up-to-date
            self._descriptiontag = None
            return
        # elif self._subifdslevel >= 0:
        #     # don't write metadata to SubIFDs
        #     return
        elif self._imagej:
            assert self._metadata is not None
            colormapped = self._colormap is not None
            isrgb = self._storedshape.samples in {3, 4}
            description = imagej_description(
                self._datashape,
                rgb=isrgb,
                colormaped=colormapped,
                **self._metadata,
            )
        elif not self._tifffile:
            self._descriptiontag = None
            return
        else:
            assert self._metadata is not None
            description = shaped_description(self._datashape, **self._metadata)

        self._descriptiontag.overwrite(description.encode(), erase=False)
        self._descriptiontag = None

    def _addtag(
        self,
        tags: list[tuple[int, bytes, bytes | None, bool]],
        code: int | str,
        dtype: int | str,
        count: int | None,
        value: Any,
        writeonce: bool = False,
        /,
    ) -> None:
        """Append (code, ifdentry, ifdvalue, writeonce) to tags list.

        Compute ifdentry and ifdvalue bytes from code, dtype, count, value.

        """
        pack = self._pack

        if not isinstance(code, int):
            code = TIFF.TAGS[code]
        try:
            datatype = cast(int, dtype)
            dataformat = TIFF.DATA_FORMATS[datatype][-1]
        except KeyError as exc:
            try:
                dataformat = cast(str, dtype)
                if dataformat[0] in '<>':
                    dataformat = dataformat[1:]
                datatype = TIFF.DATA_DTYPES[dataformat]
            except (KeyError, TypeError):
                raise ValueError(f'unknown dtype {dtype}') from exc
        del dtype

        rawcount = count
        if datatype == 2:
            # string
            if isinstance(value, str):
                # enforce 7-bit ASCII on Unicode strings
                try:
                    value = value.encode('ascii')
                except UnicodeEncodeError as exc:
                    raise ValueError(
                        'TIFF strings must be 7-bit ASCII'
                    ) from exc
            elif not isinstance(value, bytes):
                raise ValueError('TIFF strings must be 7-bit ASCII')

            if len(value) == 0 or value[-1:] != b'\x00':
                value += b'\x00'
            count = len(value)
            if code == 270:
                rawcount = int(value.find(b'\x00\x00'))
                if rawcount < 0:
                    rawcount = count
                else:
                    # length of string without buffer
                    rawcount = max(self.tiff.offsetsize + 1, rawcount + 1)
                    rawcount = min(count, rawcount)
            else:
                rawcount = count
            value = (value,)

        elif isinstance(value, bytes):
            # packed binary data
            itemsize = struct.calcsize(dataformat)
            if len(value) % itemsize:
                raise ValueError('invalid packed binary data')
            count = len(value) // itemsize
            rawcount = count

        elif count is None:
            raise ValueError('invalid count')
        else:
            count = int(count)

        if datatype in {5, 10}:  # rational
            count *= 2
            dataformat = dataformat[-1]

        ifdentry = [
            pack('HH', code, datatype),
            pack(self.tiff.offsetformat, rawcount),
        ]

        ifdvalue = None
        if struct.calcsize(dataformat) * count <= self.tiff.offsetsize:
            # value(s) can be written directly
            valueformat = f'{self.tiff.offsetsize}s'
            if isinstance(value, bytes):
                ifdentry.append(pack(valueformat, value))
            elif count == 1:
                if isinstance(value, (tuple, list, numpy.ndarray)):
                    value = value[0]
                ifdentry.append(pack(valueformat, pack(dataformat, value)))
            else:
                ifdentry.append(
                    pack(valueformat, pack(f'{count}{dataformat}', *value))
                )
        else:
            # use offset to value(s)
            ifdentry.append(pack(self.tiff.offsetformat, 0))
            if isinstance(value, bytes):
                ifdvalue = value
            elif isinstance(value, numpy.ndarray):
                if value.size != count:
                    raise RuntimeError('value.size != count')
                if value.dtype.char != dataformat:
                    raise RuntimeError('value.dtype.char != dtype')
                ifdvalue = value.tobytes()
            elif isinstance(value, (tuple, list)):
                ifdvalue = pack(f'{count}{dataformat}', *value)
            else:
                ifdvalue = pack(dataformat, value)
        tags.append((code, b''.join(ifdentry), ifdvalue, writeonce))

    def _pack(self, fmt: str, *val: Any) -> bytes:
        """Return values packed to bytes according to format."""
        if fmt[0] not in '<>':
            fmt = self.tiff.byteorder + fmt
        return struct.pack(fmt, *val)

    def _bytecount_format(
        self, bytecounts: Sequence[int], compression: int, /
    ) -> str:
        """Return small bytecount format."""
        if len(bytecounts) == 1:
            return self.tiff.offsetformat[1]
        bytecount = bytecounts[0]
        if compression > 1:
            bytecount = bytecount * 10
        if bytecount < 2**16:
            return 'H'
        if bytecount < 2**32:
            return 'I'
        return self.tiff.offsetformat[1]

    @staticmethod
    def _maxworkers(
        maxworkers: int | None,
        numchunks: int,
        chunksize: int,
        compression: int,
    ) -> int:
        """Return number of threads to encode segments."""
        if maxworkers is not None:
            return maxworkers
        if (
            # imagecodecs is None or
            compression <= 1
            or numchunks < 2
            or chunksize < 1024
            or compression == 48124  # Jetraw is not thread-safe?
        ):
            return 1
        # the following is based on benchmarking RGB tile sizes vs maxworkers
        # using a (8228, 11500, 3) uint8 WSI slide:
        if chunksize < 131072 and compression in {
            7,  # JPEG
            33007,  # ALT_JPG
            34892,  # JPEG_LOSSY
            32773,  # PackBits
            34887,  # LERC
        }:
            return 1
        if chunksize < 32768 and compression in {
            5,  # LZW
            8,  # zlib
            32946,  # zlib
            50000,  # zstd
            50013,  # zlib/pixtiff
        }:
            # zlib,
            return 1
        if chunksize < 8192 and compression in {
            34934,  # JPEG XR
            22610,  # JPEG XR
            34933,  # PNG
        }:
            return 1
        if chunksize < 2048 and compression in {
            33003,  # JPEG2000
            33004,  # JPEG2000
            33005,  # JPEG2000
            34712,  # JPEG2000
            50002,  # JPEG XL
            52546,  # JPEG XL DNG
        }:
            return 1
        if chunksize < 1024 and compression in {
            34925,  # LZMA
            50001,  # WebP
        }:
            return 1
        if compression == 34887:  # LERC
            # limit to 4 threads
            return min(numchunks, 4)
        return min(numchunks, TIFF.MAXWORKERS)

    def __enter__(self) -> TiffWriter:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        return f'<tifffile.TiffWriter {snipstr(self.filehandle.name, 32)!r}>'


@final
class TiffFile:
    """Read image and metadata from TIFF file.

    TiffFile instances must be closed with :py:meth:`TiffFile.close`, which
    is automatically called when using the 'with' context manager.

    TiffFile instances are not thread-safe. All attributes are read-only.

    Parameters:
        file:
            Specifies TIFF file to read.
            Open file objects must be positioned at the TIFF header.
        mode:
            File open mode if `file` is file name. The default is 'rb'.
        name:
            Name of file if `file` is file handle.
        offset:
            Start position of embedded file.
            The default is the current file position.
        size:
            Size of embedded file. The default is the number of bytes
            from the `offset` to the end of the file.
        omexml:
            OME metadata in XML format, for example, from external companion
            file or sanitized XML overriding XML in file.
        _multifile, _useframes, _parent:
            Internal use.
        **is_flags:
            Override `TiffFile.is_` flags, for example:

            ``is_ome=False``: disable processing of OME-XML metadata.
            ``is_lsm=False``: disable special handling of LSM files.
            ``is_ndpi=True``: force file to be NDPI format.

    Raises:
        TiffFileError: Invalid TIFF structure.

    """

    tiff: TiffFormat
    """Properties of TIFF file format."""

    pages: TiffPages
    """Sequence of pages in TIFF file."""

    _fh: FileHandle
    _multifile: bool
    _parent: TiffFile  # OME master file
    _files: dict[str | None, TiffFile]  # cache of TiffFile instances
    _omexml: str | None  # external OME-XML
    _decoders: dict[  # cache of TiffPage.decode functions
        int,
        Callable[
            ...,
            tuple[
                NDArray[Any] | None,
                tuple[int, int, int, int, int],
                tuple[int, int, int, int],
            ],
        ],
    ]

    def __init__(
        self,
        file: str | os.PathLike[Any] | FileHandle | IO[bytes],
        /,
        *,
        mode: Literal['r', 'r+'] | None = None,
        name: str | None = None,
        offset: int | None = None,
        size: int | None = None,
        omexml: str | None = None,
        _multifile: bool | None = None,
        _useframes: bool | None = None,
        _parent: TiffFile | None = None,
        **is_flags: bool | None,
    ) -> None:
        for key, value in is_flags.items():
            if key[:3] == 'is_' and key[3:] in TIFF.FILE_FLAGS:
                if value is not None:
                    setattr(self, key, bool(value))
            else:
                raise TypeError(f'unexpected keyword argument: {key}')

        if mode not in {None, 'r', 'r+', 'rb', 'r+b'}:
            raise ValueError(f'invalid mode {mode!r}')

        self._omexml = None
        if omexml:
            if omexml.strip()[-4:] != 'OME>':
                raise ValueError('invalid OME-XML')
            self._omexml = omexml
            self.is_ome = True

        fh = FileHandle(file, mode=mode, name=name, offset=offset, size=size)
        self._fh = fh
        self._multifile = True if _multifile is None else bool(_multifile)
        self._files = {fh.name: self}
        self._decoders = {}
        self._parent = self if _parent is None else _parent

        try:
            fh.seek(0)
            header = fh.read(4)
            try:
                byteorder = {b'II': '<', b'MM': '>', b'EP': '<'}[header[:2]]
            except KeyError as exc:
                raise TiffFileError(f'not a TIFF file {header!r}') from exc

            version = struct.unpack(byteorder + 'H', header[2:4])[0]
            if version == 43:
                # BigTiff
                offsetsize, zero = struct.unpack(byteorder + 'HH', fh.read(4))
                if zero != 0 or offsetsize != 8:
                    raise TiffFileError(
                        f'invalid BigTIFF offset size {(offsetsize, zero)}'
                    )
                if byteorder == '>':
                    self.tiff = TIFF.BIG_BE
                else:
                    self.tiff = TIFF.BIG_LE
            elif version == 42:
                # Classic TIFF
                if byteorder == '>':
                    self.tiff = TIFF.CLASSIC_BE
                elif is_flags.get('is_ndpi', fh.extension == '.ndpi'):
                    # NDPI uses 64 bit IFD offsets
                    if is_flags.get('is_ndpi', True):
                        self.tiff = TIFF.NDPI_LE
                    else:
                        self.tiff = TIFF.CLASSIC_LE
                else:
                    self.tiff = TIFF.CLASSIC_LE
            elif version == 0x4E31:
                # NIFF
                if byteorder == '>':
                    raise TiffFileError('invalid NIFF file')
                logger().error(f'{self!r} NIFF format not supported')
                self.tiff = TIFF.CLASSIC_LE
            elif version in {0x55, 0x4F52, 0x5352}:
                # Panasonic or Olympus RAW
                logger().error(
                    f'{self!r} RAW format 0x{version:04X} not supported'
                )
                if byteorder == '>':
                    self.tiff = TIFF.CLASSIC_BE
                else:
                    self.tiff = TIFF.CLASSIC_LE
            else:
                raise TiffFileError(f'invalid TIFF version {version}')

            # file handle is at offset to offset to first page
            self.pages = TiffPages(self)

            if self.is_lsm and (
                self.filehandle.size >= 2**32
                or self.pages[0].compression != 1
                or self.pages[1].compression != 1
            ):
                self._lsm_load_pages()

            elif self.is_scanimage and not self.is_bigtiff:
                # ScanImage <= 2015
                try:
                    self.pages._load_virtual_frames()
                except Exception as exc:
                    logger().error(
                        f'{self!r} <TiffPages._load_virtual_frames> '
                        f'raised {exc!r:.128}'
                    )

            elif self.is_ndpi:
                try:
                    self._ndpi_load_pages()
                except Exception as exc:
                    logger().error(
                        f'{self!r} <_ndpi_load_pages> raised {exc!r:.128}'
                    )

            elif _useframes:
                self.pages.useframes = True

        except Exception:
            fh.close()
            raise

    @property
    def byteorder(self) -> Literal['>', '<']:
        """Byteorder of TIFF file."""
        return self.tiff.byteorder

    @property
    def filehandle(self) -> FileHandle:
        """File handle."""
        return self._fh

    @property
    def filename(self) -> str:
        """Name of file handle."""
        return self._fh.name

    @cached_property
    def fstat(self) -> Any:
        """Status of file handle's descriptor, if any."""
        try:
            return os.fstat(self._fh.fileno())
        except Exception:  # io.UnsupportedOperation
            return None

    def close(self) -> None:
        """Close open file handle(s)."""
        for tif in self._files.values():
            tif.filehandle.close()

    def asarray(
        self,
        key: int | slice | Iterable[int] | None = None,
        *,
        series: int | TiffPageSeries | None = None,
        level: int | None = None,
        squeeze: bool | None = None,
        out: OutputType = None,
        maxworkers: int | None = None,
        buffersize: int | None = None,
    ) -> NDArray[Any]:
        """Return images from select pages as NumPy array.

        By default, the image array from the first level of the first series
        is returned.

        Parameters:
            key:
                Specifies which pages to return as array.
                By default, the image of the specified `series` and `level`
                is returned.
                If not *None*, the images from the specified pages in the
                whole file (if `series` is *None*) or a specified series are
                returned as a stacked array.
                Requesting an array from multiple pages that are not
                compatible wrt. shape, dtype, compression etc. is undefined,
                that is, it may crash or return incorrect values.
            series:
                Specifies which series of pages to return as array.
                The default is 0.
            level:
                Specifies which level of multi-resolution series to return
                as array. The default is 0.
            squeeze:
                If *True*, remove all length-1 dimensions (except X and Y)
                from array.
                If *False*, single pages are returned as 5D array of shape
                :py:attr:`TiffPage.shaped`.
                For series, the shape of the returned array also includes
                singlet dimensions specified in some file formats.
                For example, ImageJ series and most commonly also OME series,
                are returned in TZCYXS order.
                By default, all but `"shaped"` series are squeezed.
            out:
                Specifies how image array is returned.
                By default, a new NumPy array is created.
                If a *numpy.ndarray*, a writable array to which the image
                is copied.
                If *'memmap'*, directly memory-map the image data in the
                file if possible; else create a memory-mapped array in a
                temporary file.
                If a *string* or *open file*, the file used to create a
                memory-mapped array.
            maxworkers:
                Maximum number of threads to concurrently decode data from
                multiple pages or compressed segments.
                If *None* or *0*, use up to :py:attr:`_TIFF.MAXWORKERS`
                threads. Reading data from file is limited to the main thread.
                Using multiple threads can significantly speed up this
                function if the bottleneck is decoding compressed data,
                for example, in case of large LZW compressed LSM files or
                JPEG compressed tiled slides.
                If the bottleneck is I/O or pure Python code, using multiple
                threads might be detrimental.
            buffersize:
                Approximate number of bytes to read from file in one pass.
                The default is :py:attr:`_TIFF.BUFFERSIZE`.

        Returns:
            Images from specified pages. See `TiffPage.asarray`
            for operations that are applied (or not) to the image data
            stored in the file.

        """
        if not self.pages:
            return numpy.array([])
        if key is None and series is None:
            series = 0

        pages: Any  # TiffPages | TiffPageSeries | list[TiffPage | TiffFrame]
        page0: TiffPage | TiffFrame | None

        if series is None:
            pages = self.pages
        else:
            if not isinstance(series, TiffPageSeries):
                series = self.series[series]
            if level is not None:
                series = series.levels[level]
            pages = series

        if key is None:
            pass
        elif series is None:
            pages = pages._getlist(key)
        elif isinstance(key, (int, numpy.integer)):
            pages = [pages[int(key)]]
        elif isinstance(key, slice):
            pages = pages[key]
        elif isinstance(key, Iterable) and not isinstance(key, str):
            pages = [pages[k] for k in key]
        else:
            raise TypeError(
                f'key must be an integer, slice, or sequence, not {type(key)}'
            )

        if pages is None or len(pages) == 0:
            raise ValueError('no pages selected')

        if (
            key is None
            and series is not None
            and series.dataoffset is not None
        ):
            typecode = self.byteorder + series.dtype.char
            if (
                series.keyframe.is_memmappable
                and isinstance(out, str)
                and out == 'memmap'
            ):
                # direct mapping
                shape = series.get_shape(squeeze)
                result = self.filehandle.memmap_array(
                    typecode, shape, series.dataoffset
                )
            else:
                # read into output
                shape = series.get_shape(squeeze)
                if out is not None:
                    out = create_output(out, shape, series.dtype)
                result = self.filehandle.read_array(
                    typecode,
                    series.size,
                    series.dataoffset,
                    out=out,
                )
        elif len(pages) == 1:
            page0 = pages[0]
            if page0 is None:
                raise ValueError('page is None')
            result = page0.asarray(
                out=out, maxworkers=maxworkers, buffersize=buffersize
            )
        else:
            result = stack_pages(
                pages, out=out, maxworkers=maxworkers, buffersize=buffersize
            )

        assert result is not None

        if key is None:
            assert series is not None  # TODO: ?
            shape = series.get_shape(squeeze)
            try:
                result.shape = shape
            except ValueError as exc:
                try:
                    logger().warning(
                        f'{self!r} <asarray> failed to reshape '
                        f'{result.shape} to {shape}, raised {exc!r:.128}'
                    )
                    # try series of expected shapes
                    result.shape = (-1,) + shape
                except ValueError:
                    # revert to generic shape
                    result.shape = (-1,) + series.keyframe.shape
        elif len(pages) == 1:
            if squeeze is None:
                squeeze = True
            page0 = pages[0]
            if page0 is None:
                raise ValueError('page is None')
            result.shape = page0.shape if squeeze else page0.shaped
        else:
            if squeeze is None:
                squeeze = True
            try:
                page0 = next(p for p in pages if p is not None)
            except StopIteration as exc:
                raise ValueError('pages are all None') from exc
            assert page0 is not None
            result.shape = (-1,) + (page0.shape if squeeze else page0.shaped)
        return result

    def aszarr(
        self,
        key: int | None = None,
        *,
        series: int | TiffPageSeries | None = None,
        level: int | None = None,
        **kwargs: Any,
    ) -> ZarrTiffStore:
        """Return images from select pages as Zarr 2 store.

        By default, the images from the first series, including all levels,
        are wrapped as a Zarr 2 store.

        Parameters:
            key:
                Index of page in file (if `series` is None) or series to wrap
                as Zarr 2 store.
                By default, a series is wrapped.
            series:
                Index of series to wrap as Zarr 2 store.
                The default is 0 (if `key` is None).
            level:
                Index of pyramid level in series to wrap as Zarr 2 store.
                By default, all levels are included as a multi-scale group.
            **kwargs:
                Additional arguments passed to :py:meth:`TiffPage.aszarr`
                or :py:meth:`TiffPageSeries.aszarr`.

        """
        if not self.pages:
            raise NotImplementedError('empty Zarr arrays not supported')
        if key is None and series is None:
            return self.series[0].aszarr(level=level, **kwargs)

        pages: Any
        if series is None:
            pages = self.pages
        else:
            if not isinstance(series, TiffPageSeries):
                series = self.series[series]
            if key is None:
                return series.aszarr(level=level, **kwargs)
            if level is not None:
                series = series.levels[level]
            pages = series

        if isinstance(key, (int, numpy.integer)):
            page: TiffPage | TiffFrame = pages[key]
            return page.aszarr(**kwargs)
        raise TypeError('key must be an integer index')

    @cached_property
    def series(self) -> list[TiffPageSeries]:
        """Series of pages with compatible shape and data type.

        Side effect: after accessing this property, `TiffFile.pages` might
        contain `TiffPage` and `TiffFrame` instead of only `TiffPage`
        instances.

        """
        if not self.pages:
            return []
        assert self.pages.keyframe is not None
        useframes = self.pages.useframes
        keyframe = self.pages.keyframe.index
        series: list[TiffPageSeries] | None = None
        for kind in (
            'shaped',
            'lsm',
            'mmstack',
            'ome',
            'imagej',
            'ndtiff',
            'fluoview',
            'stk',
            'sis',
            'svs',
            'scn',
            'qpi',
            'ndpi',
            'bif',
            'avs',
            'philips',
            'scanimage',
            # 'indica',  # TODO: rewrite _series_indica()
            'nih',
            'mdgel',  # adds second page to cache
            'uniform',
        ):
            if getattr(self, 'is_' + kind, False):
                series = getattr(self, '_series_' + kind)()
                if not series:
                    if kind == 'ome' and self.is_imagej:
                        # try ImageJ series if OME series fails.
                        # clear pages cache since _series_ome() might leave
                        # some frames without keyframe
                        self.pages._clear()
                        continue
                    if kind == 'mmstack':
                        # try OME, ImageJ, uniform
                        continue
                break
        if not series:
            series = self._series_generic()

        self.pages.useframes = useframes
        self.pages.set_keyframe(keyframe)

        # remove empty series, for example, in MD Gel files
        # series = [s for s in series if product(s.shape) > 0]
        assert series is not None
        for i, s in enumerate(series):
            s._index = i
        return series

    def _series_uniform(self) -> list[TiffPageSeries] | None:
        """Return all images in file as single series."""
        self.pages.useframes = True
        self.pages.set_keyframe(0)
        page = self.pages.first
        validate = not (page.is_scanimage or page.is_nih)
        pages = self.pages._getlist(validate=validate)
        if len(pages) == 1:
            shape = page.shape
            axes = page.axes
        else:
            shape = (len(pages),) + page.shape
            axes = 'I' + page.axes
        dtype = page.dtype
        return [TiffPageSeries(pages, shape, dtype, axes, kind='uniform')]

    def _series_generic(self) -> list[TiffPageSeries] | None:
        """Return image series in file.

        A series is a sequence of TiffPages with the same hash.

        """
        pages = self.pages
        pages._clear(False)
        pages.useframes = False
        if pages.cache:
            pages._load()

        series = []
        keys = []
        seriesdict: dict[int, list[TiffPage | TiffFrame]] = {}

        def addpage(page: TiffPage | TiffFrame, /) -> None:
            # add page to seriesdict
            if not page.shape:  # or product(page.shape) == 0:
                return
            key = page.hash
            if key in seriesdict:
                for p in seriesdict[key]:
                    if p.offset == page.offset:
                        break  # remove duplicate page
                else:
                    seriesdict[key].append(page)
            else:
                keys.append(key)
                seriesdict[key] = [page]

        for page in pages:
            addpage(page)
            if page.subifds is not None:
                for i, offset in enumerate(page.subifds):
                    if offset < 8:
                        continue
                    try:
                        self._fh.seek(offset)
                        subifd = TiffPage(self, (page.index, i))
                    except Exception as exc:
                        logger().warning(
                            f'{self!r} generic series raised {exc!r:.128}'
                        )
                    else:
                        addpage(subifd)

        for key in keys:
            pagelist = seriesdict[key]
            page = pagelist[0]
            shape = (len(pagelist),) + page.shape
            axes = 'I' + page.axes
            if 'S' not in axes:
                shape += (1,)
                axes += 'S'
            series.append(
                TiffPageSeries(
                    pagelist, shape, page.dtype, axes, kind='generic'
                )
            )

        self.is_uniform = len(series) == 1  # replaces is_uniform method
        if not self.is_agilent:
            pyramidize_series(series)
        return series

    def _series_shaped(self) -> list[TiffPageSeries] | None:
        """Return image series in tifffile "shaped" formatted file."""
        # TODO: all series need to have JSON metadata for this to succeed

        def append(
            series: list[TiffPageSeries],
            pages: list[TiffPage | TiffFrame | None],
            axes: str | None,
            shape: tuple[int, ...] | None,
            reshape: tuple[int, ...],
            name: str,
            truncated: bool | None,
            /,
        ) -> None:
            # append TiffPageSeries to series
            assert isinstance(pages[0], TiffPage)
            page = pages[0]
            if not check_shape(page.shape, reshape):
                logger().warning(
                    f'{self!r} shaped series metadata does not match '
                    f'page shape {page.shape} != {tuple(reshape)}'
                )
                failed = True
            else:
                failed = False
            if failed or axes is None or shape is None:
                shape = page.shape
                axes = page.axes
                if len(pages) > 1:
                    shape = (len(pages),) + shape
                    axes = 'Q' + axes
                if failed:
                    reshape = shape
            size = product(shape)
            resize = product(reshape)
            if page.is_contiguous and resize > size and resize % size == 0:
                if truncated is None:
                    truncated = True
                axes = 'Q' + axes
                shape = (resize // size,) + shape
            try:
                axes = reshape_axes(axes, shape, reshape)
                shape = reshape
            except ValueError as exc:
                logger().error(
                    f'{self!r} shaped series failed to reshape, '
                    f'raised {exc!r:.128}'
                )
            series.append(
                TiffPageSeries(
                    pages,
                    shape,
                    page.dtype,
                    axes,
                    name=name,
                    kind='shaped',
                    truncated=bool(truncated),
                    squeeze=False,
                )
            )

        def detect_series(
            pages: TiffPages | list[TiffPage | TiffFrame | None],
            series: list[TiffPageSeries],
            /,
        ) -> list[TiffPageSeries] | None:
            shape: tuple[int, ...] | None
            reshape: tuple[int, ...]
            page: TiffPage | TiffFrame | None
            keyframe: TiffPage
            subifds: list[TiffPage | TiffFrame | None] = []
            subifd: TiffPage | TiffFrame
            keysubifd: TiffPage
            axes: str | None
            name: str

            lenpages = len(pages)
            index = 0
            while True:
                if index >= lenpages:
                    break

                if isinstance(pages, TiffPages):
                    # new keyframe; start of new series
                    pages.set_keyframe(index)
                    keyframe = cast(TiffPage, pages.keyframe)
                else:
                    # pages is list of SubIFDs
                    keyframe = cast(TiffPage, pages[0])

                if keyframe.shaped_description is None:
                    logger().error(
                        f'{self!r} '
                        'invalid shaped series metadata or corrupted file'
                    )
                    return None
                # read metadata
                axes = None
                shape = None
                metadata = shaped_description_metadata(
                    keyframe.shaped_description
                )
                name = metadata.get('name', '')
                reshape = metadata['shape']
                truncated = None if keyframe.subifds is None else False
                truncated = metadata.get('truncated', truncated)
                if 'axes' in metadata:
                    axes = cast(str, metadata['axes'])
                    if len(axes) == len(reshape):
                        shape = reshape
                    else:
                        axes = ''
                        logger().error(
                            f'{self!r} shaped series axes do not match shape'
                        )
                # skip pages if possible
                spages: list[TiffPage | TiffFrame | None] = [keyframe]
                size = product(reshape)
                if size > 0:
                    npages, mod = divmod(size, product(keyframe.shape))
                else:
                    npages = 1
                    mod = 0
                if mod:
                    logger().error(
                        f'{self!r} '
                        'shaped series shape does not match page shape'
                    )
                    return None

                if 1 < npages <= lenpages - index:
                    assert keyframe._dtype is not None
                    size *= keyframe._dtype.itemsize
                    if truncated:
                        npages = 1
                    else:
                        page = pages[index + 1]
                        if (
                            keyframe.is_final
                            and page is not None
                            and keyframe.offset + size < page.offset
                            and keyframe.subifds is None
                        ):
                            truncated = False
                        else:
                            # must read all pages for series
                            truncated = False
                            for j in range(index + 1, index + npages):
                                page = pages[j]
                                assert page is not None
                                page.keyframe = keyframe
                                spages.append(page)
                append(series, spages, axes, shape, reshape, name, truncated)
                index += npages

                # create series from SubIFDs
                if keyframe.subifds:
                    subifds_size = len(keyframe.subifds)
                    for i, offset in enumerate(keyframe.subifds):
                        if offset < 8:
                            continue
                        subifds = []
                        for j, page in enumerate(spages):
                            # if page.subifds is not None:
                            try:
                                if (
                                    page is None
                                    or page.subifds is None
                                    or len(page.subifds) < subifds_size
                                ):
                                    raise ValueError(
                                        f'{page!r} contains invalid subifds'
                                    )
                                self._fh.seek(page.subifds[i])
                                if j == 0:
                                    subifd = TiffPage(self, (page.index, i))
                                    keysubifd = subifd
                                else:
                                    subifd = TiffFrame(
                                        self,
                                        (page.index, i),
                                        keyframe=keysubifd,
                                    )
                            except Exception as exc:
                                logger().error(
                                    f'{self!r} shaped series '
                                    f'raised {exc!r:.128}'
                                )
                                return None
                            subifds.append(subifd)
                        if subifds:
                            series_or_none = detect_series(subifds, series)
                            if series_or_none is None:
                                return None
                            series = series_or_none
            return series

        self.pages.useframes = True
        series = detect_series(self.pages, [])
        if series is None:
            return None
        self.is_uniform = len(series) == 1
        pyramidize_series(series, isreduced=True)
        return series

    def _series_imagej(self) -> list[TiffPageSeries] | None:
        """Return image series in ImageJ file."""
        # ImageJ's dimension order is TZCYXS
        # TODO: fix loading of color, composite, or palette images
        meta = self.imagej_metadata
        if meta is None:
            return None

        pages = self.pages
        pages.useframes = True
        pages.set_keyframe(0)
        page = self.pages.first

        order = meta.get('order', 'czt').lower()
        frames = meta.get('frames', 1)
        slices = meta.get('slices', 1)
        channels = meta.get('channels', 1)
        images = meta.get('images', 1)  # not reliable

        if images < 1 or frames < 1 or slices < 1 or channels < 1:
            logger().warning(
                f'{self!r} ImageJ series metadata invalid or corrupted file'
            )
            return None

        if channels == 1:
            images = frames * slices
        elif page.shaped[0] > 1 and page.shaped[0] == channels:
            # Bio-Formats declares separate samples as channels
            images = frames * slices
        elif images == frames * slices and page.shaped[4] == channels:
            # RGB contig samples declared as channel
            channels = 1
        else:
            images = frames * slices * channels

        if images == 1 and pages.is_multipage:
            images = len(pages)

        nbytes = images * page.nbytes

        # ImageJ virtual hyperstacks store all image metadata in the first
        # page and image data are stored contiguously before the second
        # page, if any
        if not page.is_final:
            isvirtual = False
        elif page.dataoffsets[0] + nbytes > self.filehandle.size:
            logger().error(
                f'{self!r} ImageJ series metadata invalid or corrupted file'
            )
            return None
        elif images <= 1:
            isvirtual = True
        elif (
            pages.is_multipage
            and page.dataoffsets[0] + nbytes > pages[1].offset
        ):
            # next page is not stored after data
            isvirtual = False
        else:
            isvirtual = True

        page_list: list[TiffPage | TiffFrame]
        if isvirtual:
            # no need to read other pages
            page_list = [page]
        else:
            page_list = pages[:]

        shape: tuple[int, ...]
        axes: str

        if order in {'czt', 'default'}:
            axes = 'TZC'
            shape = (frames, slices, channels)
        elif order == 'ctz':
            axes = 'ZTC'
            shape = (slices, frames, channels)
        elif order == 'zct':
            axes = 'TCZ'
            shape = (frames, channels, slices)
        elif order == 'ztc':
            axes = 'CTZ'
            shape = (channels, frames, slices)
        elif order == 'tcz':
            axes = 'ZCT'
            shape = (slices, channels, frames)
        elif order == 'tzc':
            axes = 'CZT'
            shape = (channels, slices, frames)
        else:
            axes = 'TZC'
            shape = (frames, slices, channels)
            logger().warning(
                f'{self!r} ImageJ series of unknown order {order!r}'
            )

        remain = images // product(shape)
        if remain > 1:
            logger().debug(
                f'{self!r} ImageJ series contains unidentified dimension'
            )
            shape = (remain,) + shape
            axes = 'I' + axes

        if page.shaped[0] > 1:
            # Bio-Formats declares separate samples as channels
            assert axes[-1] == 'C'
            shape = shape[:-1] + page.shape
            axes += page.axes[1:]
        else:
            shape += page.shape
            axes += page.axes

        if 'S' not in axes:
            shape += (1,)
            axes += 'S'
        # assert axes.endswith('TZCYXS'), axes

        truncated = (
            isvirtual and not pages.is_multipage and page.nbytes != nbytes
        )

        self.is_uniform = True
        return [
            TiffPageSeries(
                page_list,
                shape,
                page.dtype,
                axes,
                kind='imagej',
                truncated=truncated,
            )
        ]

    def _series_nih(self) -> list[TiffPageSeries] | None:
        """Return all images in NIH Image file as single series."""
        series = self._series_uniform()
        if series is not None:
            for s in series:
                s.kind = 'nih'
        return series

    def _series_scanimage(self) -> list[TiffPageSeries] | None:
        """Return image series in ScanImage file."""
        pages = self.pages._getlist(validate=False)
        page = self.pages.first
        dtype = page.dtype
        shape = None

        meta = self.scanimage_metadata
        if meta is None:
            framedata = {}
        else:
            framedata = meta.get('FrameData', {})
        if 'SI.hChannels.channelSave' in framedata:
            try:
                channels = framedata['SI.hChannels.channelSave']
                try:
                    # channelSave is a list of channel IDs
                    channels = len(channels)
                except TypeError:
                    # channelSave is a single channel ID
                    channels = 1
                # slices = framedata.get(
                #    'SI.hStackManager.actualNumSlices',
                #     framedata.get('SI.hStackManager.numSlices', None),
                # )
                # if slices is None:
                #     raise ValueError('unable to determine numSlices')
                slices = None
                try:
                    frames = int(framedata['SI.hStackManager.framesPerSlice'])
                except Exception as exc:
                    # framesPerSlice is inf
                    slices = 1
                    if len(pages) % channels:
                        raise ValueError(
                            'unable to determine framesPerSlice'
                        ) from exc
                    frames = len(pages) // channels
                if slices is None:
                    slices = max(len(pages) // (frames * channels), 1)
                shape = (slices, frames, channels) + page.shape
                axes = 'ZTC' + page.axes
            except Exception as exc:
                logger().warning(
                    f'{self!r} ScanImage series raised {exc!r:.128}'
                )

        # TODO: older versions of ScanImage store non-varying frame data in
        # the ImageDescription tag. Candidates are scanimage.SI5.channelsSave,
        # scanimage.SI5.stackNumSlices, scanimage.SI5.acqNumFrames
        # scanimage.SI4., state.acq.numberOfFrames, state.acq.numberOfFrames...

        if shape is None:
            shape = (len(pages),) + page.shape
            axes = 'I' + page.axes

        return [TiffPageSeries(pages, shape, dtype, axes, kind='scanimage')]

    def _series_fluoview(self) -> list[TiffPageSeries] | None:
        """Return image series in FluoView file."""
        meta = self.fluoview_metadata
        if meta is None:
            return None
        pages = self.pages._getlist(validate=False)
        mmhd = list(reversed(meta['Dimensions']))
        axes = ''.join(TIFF.MM_DIMENSIONS.get(i[0].upper(), 'Q') for i in mmhd)
        shape = tuple(int(i[1]) for i in mmhd)
        self.is_uniform = True
        return [
            TiffPageSeries(
                pages,
                shape,
                pages[0].dtype,
                axes,
                name=meta['ImageName'],
                kind='fluoview',
            )
        ]

    def _series_mdgel(self) -> list[TiffPageSeries] | None:
        """Return image series in MD Gel file."""
        # only a single page, scaled according to metadata in second page
        meta = self.mdgel_metadata
        if meta is None:
            return None
        transform: Callable[[NDArray[Any]], NDArray[Any]] | None
        self.pages.useframes = False
        self.pages.set_keyframe(0)

        if meta['FileTag'] in {2, 128}:
            dtype = numpy.dtype(numpy.float32)
            scale = meta['ScalePixel']
            scale = scale[0] / scale[1]  # rational
            if meta['FileTag'] == 2:
                # squary root data format
                def transform(a: NDArray[Any], /) -> NDArray[Any]:
                    return a.astype(numpy.float32) ** 2 * scale

            else:

                def transform(a: NDArray[Any], /) -> NDArray[Any]:
                    return a.astype(numpy.float32) * scale

        else:
            transform = None
        page = self.pages.first
        self.is_uniform = False
        return [
            TiffPageSeries(
                [page],
                page.shape,
                dtype,
                page.axes,
                transform=transform,
                kind='mdgel',
            )
        ]

    def _series_ndpi(self) -> list[TiffPageSeries] | None:
        """Return pyramidal image series in NDPI file."""
        series = self._series_generic()
        if series is None:
            return None
        for s in series:
            s.kind = 'ndpi'
            if s.axes[0] == 'I':
                s._set_dimensions(s.shape, 'Z' + s.axes[1:], None, True)
            if s.is_pyramidal:
                name = s.keyframe.tags.valueof(65427)
                s.name = 'Baseline' if name is None else name
                continue
            mag = s.keyframe.tags.valueof(65421)
            if mag is not None:
                if mag == -1.0:
                    s.name = 'Macro'
                    # s.kind += '_macro'
                elif mag == -2.0:
                    s.name = 'Map'
                    # s.kind += '_map'
        self.is_uniform = False
        return series

    def _series_avs(self) -> list[TiffPageSeries] | None:
        """Return pyramidal image series in AVS file."""
        series = self._series_generic()
        if series is None:
            return None
        if len(series) != 3:
            logger().warning(
                f'{self!r} AVS series expected 3 series, got {len(series)}'
            )
        s = series[0]
        s.kind = 'avs'
        if s.axes[0] == 'I':
            s._set_dimensions(s.shape, 'Z' + s.axes[1:], None, True)
        if s.is_pyramidal:
            s.name = 'Baseline'
        if len(series) == 3:
            series[1].name = 'Map'
            series[1].kind = 'avs'
            series[2].name = 'Macro'
            series[2].kind = 'avs'
        self.is_uniform = False
        return series

    def _series_philips(self) -> list[TiffPageSeries] | None:
        """Return pyramidal image series in Philips DP file."""
        from xml.etree import ElementTree as etree

        series = []
        pages = self.pages
        pages.cache = False
        pages.useframes = False
        pages.set_keyframe(0)
        pages._load()

        meta = self.philips_metadata
        assert meta is not None

        try:
            tree = etree.fromstring(meta)
        except etree.ParseError as exc:
            logger().error(f'{self!r} Philips series raised {exc!r:.128}')
            return None

        pixel_spacing = [
            tuple(float(v) for v in elem.text.replace('"', '').split())
            for elem in tree.findall(
                './/*'
                '/DataObject[@ObjectType="PixelDataRepresentation"]'
                '/Attribute[@Name="DICOM_PIXEL_SPACING"]'
            )
            if elem.text is not None
        ]
        if len(pixel_spacing) < 2:
            logger().error(
                f'{self!r} Philips series {len(pixel_spacing)=} < 2'
            )
            return None

        series_dict: dict[str, list[TiffPage]] = {}
        series_dict['Level'] = []
        series_dict['Other'] = []
        for page in pages:
            assert isinstance(page, TiffPage)
            if page.description.startswith('Macro'):
                series_dict['Macro'] = [page]
            elif page.description.startswith('Label'):
                series_dict['Label'] = [page]
            elif not page.is_tiled:
                series_dict['Other'].append(page)
            else:
                series_dict['Level'].append(page)

        levels = series_dict.pop('Level')
        if len(levels) != len(pixel_spacing):
            logger().error(
                f'{self!r} Philips series '
                f'{len(levels)=} != {len(pixel_spacing)=}'
            )
            return None

        # fix padding of sublevels
        imagewidth0 = levels[0].imagewidth
        imagelength0 = levels[0].imagelength
        h0, w0 = pixel_spacing[0]
        for serie, (h, w) in zip(levels[1:], pixel_spacing[1:]):
            page = serie.keyframe
            # if page.dtype.itemsize == 1:
            #     page.nodata = 255

            imagewidth = imagewidth0 // int(round(w / w0))
            imagelength = imagelength0 // int(round(h / h0))

            if page.imagewidth - page.tilewidth >= imagewidth:
                logger().warning(
                    f'{self!r} Philips series {page.index=} '
                    f'{page.imagewidth=}-{page.tilewidth=} >= {imagewidth=}'
                )
                page.imagewidth -= page.tilewidth - 1
            elif page.imagewidth < imagewidth:
                logger().warning(
                    f'{self!r} Philips series {page.index=} '
                    f'{page.imagewidth=} < {imagewidth=}'
                )
            else:
                page.imagewidth = imagewidth
            imagewidth = page.imagewidth

            if page.imagelength - page.tilelength >= imagelength:
                logger().warning(
                    f'{self!r} Philips series {page.index=} '
                    f'{page.imagelength=}-{page.tilelength=} >= {imagelength=}'
                )
                page.imagelength -= page.tilelength - 1
            # elif page.imagelength < imagelength:
            #    # in this case image is padded with zero
            else:
                page.imagelength = imagelength
            imagelength = page.imagelength

            if page.shaped[-1] > 1:
                page.shape = (imagelength, imagewidth, page.shape[-1])
            elif page.shaped[0] > 1:
                page.shape = (page.shape[0], imagelength, imagewidth)
            else:
                page.shape = (imagelength, imagewidth)
            page.shaped = (
                page.shaped[:2] + (imagelength, imagewidth) + page.shaped[-1:]
            )

        series = [TiffPageSeries([levels[0]], name='Baseline', kind='philips')]
        for i, page in enumerate(levels[1:]):
            series[0].levels.append(
                TiffPageSeries([page], name=f'Level{i + 1}', kind='philips')
            )
        for key, value in series_dict.items():
            for page in value:
                series.append(TiffPageSeries([page], name=key, kind='philips'))

        self.is_uniform = False
        return series

    def _series_indica(self) -> list[TiffPageSeries] | None:
        """Return pyramidal image series in IndicaLabs file."""
        # TODO: need more IndicaLabs sample files
        # TODO: parse indica series from XML
        # TODO: alpha channels in SubIFDs or main IFDs

        from xml.etree import ElementTree as etree

        series = self._series_generic()
        if series is None or len(series) != 1:
            return series

        try:
            tree = etree.fromstring(self.pages.first.description)
        except etree.ParseError as exc:
            logger().error(f'{self!r} Indica series raised {exc!r:.128}')
            return series

        channel_names = [
            channel.attrib['name'] for channel in tree.iter('channel')
        ]
        for s in series:
            s.kind = 'indica'
            # TODO: identify other dimensions
            if s.axes[0] == 'I' and s.shape[0] == len(channel_names):
                s._set_dimensions(s.shape, 'C' + s.axes[1:], None, True)
            if s.is_pyramidal:
                s.name = 'Baseline'
        self.is_uniform = False
        return series

    def _series_sis(self) -> list[TiffPageSeries] | None:
        """Return image series in Olympus SIS file."""
        meta = self.sis_metadata
        if meta is None:
            return None
        pages = self.pages._getlist(validate=False)  # TODO: this fails for VSI
        page = pages[0]
        lenpages = len(pages)

        if 'shape' in meta and 'axes' in meta:
            shape = meta['shape'] + page.shape
            axes = meta['axes'] + page.axes
        else:
            shape = (lenpages,) + page.shape
            axes = 'I' + page.axes
        self.is_uniform = True
        return [TiffPageSeries(pages, shape, page.dtype, axes, kind='sis')]

    def _series_qpi(self) -> list[TiffPageSeries] | None:
        """Return image series in PerkinElmer QPI file."""
        series = []
        pages = self.pages
        pages.cache = True
        pages.useframes = False
        pages.set_keyframe(0)
        pages._load()
        page0 = self.pages.first

        # Baseline
        # TODO: get name from ImageDescription XML
        ifds = []
        index = 0
        axes = 'C' + page0.axes
        dtype = page0.dtype
        pshape = page0.shape
        while index < len(pages):
            page = pages[index]
            if page.shape != pshape:
                break
            ifds.append(page)
            index += 1
        shape = (len(ifds),) + pshape
        series.append(
            TiffPageSeries(
                ifds, shape, dtype, axes, name='Baseline', kind='qpi'
            )
        )

        if index < len(pages):
            # Thumbnail
            page = pages[index]
            series.append(
                TiffPageSeries(
                    [page],
                    page.shape,
                    page.dtype,
                    page.axes,
                    name='Thumbnail',
                    kind='qpi',
                )
            )
            index += 1

        if page0.is_tiled:
            # Resolutions
            while index < len(pages):
                pshape = (pshape[0] // 2, pshape[1] // 2) + pshape[2:]
                ifds = []
                while index < len(pages):
                    page = pages[index]
                    if page.shape != pshape:
                        break
                    ifds.append(page)
                    index += 1
                if len(ifds) != len(series[0].pages):
                    break
                shape = (len(ifds),) + pshape
                series[0].levels.append(
                    TiffPageSeries(
                        ifds, shape, dtype, axes, name='Resolution', kind='qpi'
                    )
                )

        if series[0].is_pyramidal and index < len(pages):
            # Macro
            page = pages[index]
            series.append(
                TiffPageSeries(
                    [page],
                    page.shape,
                    page.dtype,
                    page.axes,
                    name='Macro',
                    kind='qpi',
                )
            )
            index += 1
            # Label
            if index < len(pages):
                page = pages[index]
                series.append(
                    TiffPageSeries(
                        [page],
                        page.shape,
                        page.dtype,
                        page.axes,
                        name='Label',
                        kind='qpi',
                    )
                )

        self.is_uniform = False
        return series

    def _series_svs(self) -> list[TiffPageSeries] | None:
        """Return image series in Aperio SVS file."""
        if not self.pages.first.is_tiled:
            return None

        series = []
        self.pages.cache = True
        self.pages.useframes = False
        self.pages.set_keyframe(0)
        self.pages._load()

        # baseline
        firstpage = self.pages.first
        if len(self.pages) == 1:
            self.is_uniform = False
            return [
                TiffPageSeries(
                    [firstpage],
                    firstpage.shape,
                    firstpage.dtype,
                    firstpage.axes,
                    name='Baseline',
                    kind='svs',
                )
            ]

        # thumbnail
        page = self.pages[1]
        thumnail = TiffPageSeries(
            [page],
            page.shape,
            page.dtype,
            page.axes,
            name='Thumbnail',
            kind='svs',
        )

        # resolutions and focal planes
        levels = {firstpage.shape: [firstpage]}
        index = 2
        while index < len(self.pages):
            page = cast(TiffPage, self.pages[index])
            if not page.is_tiled or page.is_reduced:
                break
            if page.shape in levels:
                levels[page.shape].append(page)
            else:
                levels[page.shape] = [page]
            index += 1

        zsize = len(levels[firstpage.shape])
        if not all(len(level) == zsize for level in levels.values()):
            logger().warning(f'{self!r} SVS series focal planes do not match')
            zsize = 1
        baseline = TiffPageSeries(
            levels[firstpage.shape],
            (zsize,) + firstpage.shape,
            firstpage.dtype,
            'Z' + firstpage.axes,
            name='Baseline',
            kind='svs',
        )
        for shape, level in levels.items():
            if shape == firstpage.shape:
                continue
            page = level[0]
            baseline.levels.append(
                TiffPageSeries(
                    level,
                    (zsize,) + page.shape,
                    page.dtype,
                    'Z' + page.axes,
                    name='Resolution',
                    kind='svs',
                )
            )
        series.append(baseline)
        series.append(thumnail)

        # Label, Macro; subfiletype 1, 9
        for _ in range(2):
            if index == len(self.pages):
                break
            page = self.pages[index]
            assert isinstance(page, TiffPage)
            if page.subfiletype == 9:
                name = 'Macro'
            else:
                name = 'Label'
            series.append(
                TiffPageSeries(
                    [page],
                    page.shape,
                    page.dtype,
                    page.axes,
                    name=name,
                    kind='svs',
                )
            )
            index += 1
        self.is_uniform = False
        return series

    def _series_scn(self) -> list[TiffPageSeries] | None:
        """Return pyramidal image series in Leica SCN file."""
        # TODO: support collections
        from xml.etree import ElementTree as etree

        scnxml = self.pages.first.description
        root = etree.fromstring(scnxml)

        series = []
        self.pages.cache = True
        self.pages.useframes = False
        self.pages.set_keyframe(0)
        self.pages._load()

        for collection in root:
            if not collection.tag.endswith('collection'):
                continue
            for image in collection:
                if not image.tag.endswith('image'):
                    continue
                name = image.attrib.get('name', 'Unknown')
                for pixels in image:
                    if not pixels.tag.endswith('pixels'):
                        continue
                    resolutions: dict[int, dict[str, Any]] = {}
                    for dimension in pixels:
                        if not dimension.tag.endswith('dimension'):
                            continue
                        if int(image.attrib.get('sizeZ', 1)) > 1:
                            raise NotImplementedError(
                                'SCN series: Z-Stacks not supported. '
                                'Please submit a sample file.'
                            )
                        sizex = int(dimension.attrib['sizeX'])
                        sizey = int(dimension.attrib['sizeY'])
                        c = int(dimension.attrib.get('c', 0))
                        z = int(dimension.attrib.get('z', 0))
                        r = int(dimension.attrib.get('r', 0))
                        ifd = int(dimension.attrib['ifd'])
                        if r in resolutions:
                            level = resolutions[r]
                            if c > level['channels']:
                                level['channels'] = c
                            if z > level['sizez']:
                                level['sizez'] = z
                            level['ifds'][(c, z)] = ifd
                        else:
                            resolutions[r] = {
                                'size': [sizey, sizex],
                                'channels': c,
                                'sizez': z,
                                'ifds': {(c, z): ifd},
                            }
                    if not resolutions:
                        continue
                    levels = []
                    for r, level in sorted(resolutions.items()):
                        shape: tuple[int, ...] = (
                            level['channels'] + 1,
                            level['sizez'] + 1,
                        )
                        axes = 'CZ'

                        ifds: list[TiffPage | TiffFrame | None] = [
                            None
                        ] * product(shape)
                        for (c, z), ifd in sorted(level['ifds'].items()):
                            ifds[c * shape[1] + z] = self.pages[ifd]

                        assert ifds[0] is not None
                        axes += ifds[0].axes
                        shape += ifds[0].shape
                        dtype = ifds[0].dtype

                        levels.append(
                            TiffPageSeries(
                                ifds,
                                shape,
                                dtype,
                                axes,
                                parent=self,
                                name=name,
                                kind='scn',
                            )
                        )
                    levels[0].levels.extend(levels[1:])
                    series.append(levels[0])

        self.is_uniform = False
        return series

    def _series_bif(self) -> list[TiffPageSeries] | None:
        """Return image series in Ventana/Roche BIF file."""
        series = []
        baseline: TiffPageSeries | None = None
        self.pages.cache = True
        self.pages.useframes = False
        self.pages.set_keyframe(0)
        self.pages._load()

        for page in self.pages:
            page = cast(TiffPage, page)
            if page.description[:5] == 'Label':
                series.append(
                    TiffPageSeries(
                        [page],
                        page.shape,
                        page.dtype,
                        page.axes,
                        name='Label',
                        kind='bif',
                    )
                )
            elif (
                page.description == 'Thumbnail'
                or page.description[:11] == 'Probability'
            ):
                series.append(
                    TiffPageSeries(
                        [page],
                        page.shape,
                        page.dtype,
                        page.axes,
                        name='Thumbnail',
                        kind='bif',
                    )
                )
            elif 'level' not in page.description:
                # TODO: is this necessary?
                series.append(
                    TiffPageSeries(
                        [page],
                        page.shape,
                        page.dtype,
                        page.axes,
                        name='Unknown',
                        kind='bif',
                    )
                )
            elif baseline is None:
                baseline = TiffPageSeries(
                    [page],
                    page.shape,
                    page.dtype,
                    page.axes,
                    name='Baseline',
                    kind='bif',
                )
                series.insert(0, baseline)
            else:
                baseline.levels.append(
                    TiffPageSeries(
                        [page],
                        page.shape,
                        page.dtype,
                        page.axes,
                        name='Resolution',
                        kind='bif',
                    )
                )

        logger().warning(f'{self!r} BIF series tiles are not stiched')
        self.is_uniform = False
        return series

    def _series_ome(self) -> list[TiffPageSeries] | None:
        """Return image series in OME-TIFF file(s)."""
        # xml.etree found to be faster than lxml
        from xml.etree import ElementTree as etree

        omexml = self.ome_metadata
        if omexml is None:
            return None
        try:
            root = etree.fromstring(omexml)
        except etree.ParseError as exc:
            # TODO: test badly encoded OME-XML
            logger().error(f'{self!r} OME series raised {exc!r:.128}')
            return None

        keyframe: TiffPage
        ifds: list[TiffPage | TiffFrame | None]
        size: int = -1

        def load_pages(tif: TiffFile, /) -> None:
            tif.pages.cache = True
            tif.pages.useframes = True
            tif.pages.set_keyframe(0)
            tif.pages._load(None)

        load_pages(self)

        root_uuid = root.attrib.get('UUID', None)
        self._files = {root_uuid: self}
        dirname = self._fh.dirname
        files_missing = 0
        moduloref = []
        modulo: dict[str, dict[str, tuple[str, int]]] = {}
        series: list[TiffPageSeries] = []
        for element in root:
            if element.tag.endswith('BinaryOnly'):
                # TODO: load OME-XML from master or companion file
                logger().debug(
                    f'{self!r} OME series is BinaryOnly, '
                    'not an OME-TIFF master file'
                )
                break
            if element.tag.endswith('StructuredAnnotations'):
                for annot in element:
                    if not annot.attrib.get('Namespace', '').endswith(
                        'modulo'
                    ):
                        continue
                    modulo[annot.attrib['ID']] = mod = {}
                    for value in annot:
                        for modulo_ns in value:
                            for along in modulo_ns:
                                if not along.tag[:-1].endswith('Along'):
                                    continue
                                axis = along.tag[-1]
                                newaxis = along.attrib.get('Type', 'other')
                                newaxis = TIFF.AXES_CODES[newaxis]
                                if 'Start' in along.attrib:
                                    step = float(along.attrib.get('Step', 1))
                                    start = float(along.attrib['Start'])
                                    stop = float(along.attrib['End']) + step
                                    labels = len(
                                        numpy.arange(start, stop, step)
                                    )
                                else:
                                    labels = len(
                                        [
                                            label
                                            for label in along
                                            if label.tag.endswith('Label')
                                        ]
                                    )
                                mod[axis] = (newaxis, labels)

            if not element.tag.endswith('Image'):
                continue

            for annot in element:
                if annot.tag.endswith('AnnotationRef'):
                    annotationref = annot.attrib['ID']
                    break
            else:
                annotationref = None

            attr = element.attrib
            name = attr.get('Name', None)

            for pixels in element:
                if not pixels.tag.endswith('Pixels'):
                    continue
                attr = pixels.attrib
                # dtype = attr.get('PixelType', None)
                axes = ''.join(reversed(attr['DimensionOrder']))
                shape = [int(attr['Size' + ax]) for ax in axes]
                ifds = []
                spp = 1  # samples per pixel
                first = True

                for data in pixels:
                    if data.tag.endswith('Channel'):
                        attr = data.attrib
                        if first:
                            first = False
                            spp = int(attr.get('SamplesPerPixel', spp))
                            if spp > 1:
                                # correct channel dimension for spp
                                shape = [
                                    shape[i] // spp if ax == 'C' else shape[i]
                                    for i, ax in enumerate(axes)
                                ]
                        elif int(attr.get('SamplesPerPixel', 1)) != spp:
                            raise ValueError(
                                'OME series cannot handle differing '
                                'SamplesPerPixel'
                            )
                        continue

                    if not data.tag.endswith('TiffData'):
                        continue

                    attr = data.attrib
                    ifd_index = int(attr.get('IFD', 0))
                    num = int(attr.get('NumPlanes', 1 if 'IFD' in attr else 0))
                    num = int(attr.get('PlaneCount', num))
                    idxs = [int(attr.get('First' + ax, 0)) for ax in axes[:-2]]
                    try:
                        idx = int(numpy.ravel_multi_index(idxs, shape[:-2]))
                    except ValueError as exc:
                        # ImageJ produces invalid ome-xml when cropping
                        logger().warning(
                            f'{self!r} '
                            'OME series contains invalid TiffData index, '
                            f'raised {exc!r:.128}',
                        )
                        continue
                    for uuid in data:
                        if not uuid.tag.endswith('UUID'):
                            continue
                        if (
                            root_uuid is None
                            and uuid.text is not None
                            and (
                                uuid.attrib.get('FileName', '').lower()
                                == self.filename.lower()
                            )
                        ):
                            # no global UUID, use this file
                            root_uuid = uuid.text
                            self._files[root_uuid] = self._files[None]
                            del self._files[None]
                        elif uuid.text not in self._files:
                            if not self._multifile:
                                # abort reading multifile OME series
                                # and fall back to generic series
                                return []
                            fname = uuid.attrib['FileName']
                            try:
                                if not self.filehandle.is_file:
                                    raise ValueError
                                tif = TiffFile(
                                    os.path.join(dirname, fname), _parent=self
                                )
                                load_pages(tif)
                            except (
                                OSError,
                                FileNotFoundError,
                                ValueError,
                            ) as exc:
                                if files_missing == 0:
                                    logger().warning(
                                        f'{self!r} OME series failed to read '
                                        f'{fname!r}, raised {exc!r:.128}. '
                                        'Missing data are zeroed'
                                    )
                                files_missing += 1
                                # assume that size is same as in previous file
                                # if no NumPlanes or PlaneCount are given
                                if num:
                                    size = num
                                elif size == -1:
                                    raise ValueError(
                                        'OME series missing '
                                        'NumPlanes or PlaneCount'
                                    ) from exc
                                ifds.extend([None] * (size + idx - len(ifds)))
                                break
                            self._files[uuid.text] = tif
                            tif.close()
                        pages = self._files[uuid.text].pages
                        try:
                            size = num if num else len(pages)
                            ifds.extend([None] * (size + idx - len(ifds)))
                            for i in range(size):
                                ifds[idx + i] = pages[ifd_index + i]
                        except IndexError as exc:
                            logger().warning(
                                f'{self!r} '
                                'OME series contains index out of range, '
                                f'raised {exc!r:.128}'
                            )
                        # only process first UUID
                        break
                    else:
                        # no uuid found
                        pages = self.pages
                        try:
                            size = num if num else len(pages)
                            ifds.extend([None] * (size + idx - len(ifds)))
                            for i in range(size):
                                ifds[idx + i] = pages[ifd_index + i]
                        except IndexError as exc:
                            logger().warning(
                                f'{self!r} '
                                'OME series contains index out of range, '
                                f'raised {exc!r:.128}'
                            )

                if not ifds or all(i is None for i in ifds):
                    # skip images without data
                    continue

                # find a keyframe
                for ifd in ifds:
                    # try find a TiffPage
                    if ifd is not None and ifd == ifd.keyframe:
                        keyframe = cast(TiffPage, ifd)
                        break
                else:
                    # reload a TiffPage from file
                    for i, ifd in enumerate(ifds):
                        if ifd is not None:
                            isclosed = ifd.parent.filehandle.closed
                            if isclosed:
                                ifd.parent.filehandle.open()
                            ifd.parent.pages.set_keyframe(ifd.index)
                            keyframe = cast(
                                TiffPage, ifd.parent.pages[ifd.index]
                            )
                            ifds[i] = keyframe
                            if isclosed:
                                keyframe.parent.filehandle.close()
                            break

                # does the series spawn multiple files
                multifile = False
                for ifd in ifds:
                    if ifd and ifd.parent != keyframe.parent:
                        multifile = True
                        break

                if spp > 1:
                    if keyframe.planarconfig == 1:
                        shape += [spp]
                        axes += 'S'
                    else:
                        shape = shape[:-2] + [spp] + shape[-2:]
                        axes = axes[:-2] + 'S' + axes[-2:]
                if 'S' not in axes:
                    shape += [1]
                    axes += 'S'

                # number of pages in the file might mismatch XML metadata, for
                # example Nikon-cell011.ome.tif or stack_t24_y2048_x2448.tiff
                size = max(product(shape) // keyframe.size, 1)
                if size < len(ifds):
                    logger().warning(
                        f'{self!r} '
                        f'OME series expected {size} frames, got {len(ifds)}'
                    )
                    ifds = ifds[:size]
                elif size > len(ifds):
                    logger().warning(
                        f'{self!r} '
                        f'OME series is missing {size - len(ifds)} frames.'
                        ' Missing data are zeroed'
                    )
                    ifds.extend([None] * (size - len(ifds)))

                # FIXME: this implementation assumes the last dimensions are
                # stored in TIFF pages. Apparently that is not always the case.
                # For example, TCX (20000, 2, 500) is stored in 2 pages of
                # (20000, 500) in 'Image 7.ome_h00.tiff'.
                # For now, verify that shapes of keyframe and series match.
                # If not, skip series.
                squeezed = squeeze_axes(shape, axes)[0]
                if keyframe.shape != tuple(squeezed[-len(keyframe.shape) :]):
                    logger().warning(
                        f'{self!r} OME series cannot handle discontiguous '
                        f'storage ({keyframe.shape} != '
                        f'{tuple(squeezed[-len(keyframe.shape) :])})',
                    )
                    del ifds
                    continue

                # set keyframe on all IFDs
                # each series must contain a TiffPage used as keyframe
                keyframes: dict[str, TiffPage] = {
                    keyframe.parent.filehandle.name: keyframe
                }
                for i, page in enumerate(ifds):
                    if page is None:
                        continue
                    fh = page.parent.filehandle
                    if fh.name not in keyframes:
                        if page.keyframe != page:
                            # reload TiffPage from file
                            isclosed = fh.closed
                            if isclosed:
                                fh.open()
                            page.parent.pages.set_keyframe(page.index)
                            page = page.parent.pages[page.index]
                            ifds[i] = page
                            if isclosed:
                                fh.close()
                        keyframes[fh.name] = cast(TiffPage, page)
                    if page.keyframe != page:
                        page.keyframe = keyframes[fh.name]

                moduloref.append(annotationref)
                series.append(
                    TiffPageSeries(
                        ifds,
                        shape,
                        keyframe.dtype,
                        axes,
                        parent=self,
                        name=name,
                        multifile=multifile,
                        kind='ome',
                    )
                )
                del ifds

        if files_missing > 1:
            logger().warning(
                f'{self!r} OME series failed to read {files_missing} files'
            )

        # apply modulo according to AnnotationRef
        for aseries, annotationref in zip(series, moduloref):
            if annotationref not in modulo:
                continue
            shape = list(aseries.get_shape(False))
            axes = aseries.get_axes(False)
            for axis, (newaxis, size) in modulo[annotationref].items():
                i = axes.index(axis)
                if shape[i] == size:
                    axes = axes.replace(axis, newaxis, 1)
                else:
                    shape[i] //= size
                    shape.insert(i + 1, size)
                    axes = axes.replace(axis, axis + newaxis, 1)
            aseries._set_dimensions(shape, axes, None)

        # pyramids
        for aseries in series:
            keyframe = aseries.keyframe
            if keyframe.subifds is None:
                continue
            if len(self._files) > 1:
                # TODO: support multi-file pyramids; must re-open/close
                logger().warning(
                    f'{self!r} OME series cannot read multi-file pyramids'
                )
                break
            for level in range(len(keyframe.subifds)):
                found_keyframe = False
                ifds = []
                for page in aseries.pages:
                    if (
                        page is None
                        or page.subifds is None
                        or page.subifds[level] < 8
                    ):
                        ifds.append(None)
                        continue
                    page.parent.filehandle.seek(page.subifds[level])
                    if page.keyframe == page:
                        ifd = keyframe = TiffPage(
                            self, (page.index, level + 1)
                        )
                        found_keyframe = True
                    elif not found_keyframe:
                        raise RuntimeError('no keyframe found')
                    else:
                        ifd = TiffFrame(
                            self, (page.index, level + 1), keyframe=keyframe
                        )
                    ifds.append(ifd)
                if all(ifd_or_none is None for ifd_or_none in ifds):
                    logger().warning(
                        f'{self!r} OME series level {level + 1} is empty'
                    )
                    break
                # fix shape
                shape = list(aseries.get_shape(False))
                axes = aseries.get_axes(False)
                for i, ax in enumerate(axes):
                    if ax == 'X':
                        shape[i] = keyframe.imagewidth
                    elif ax == 'Y':
                        shape[i] = keyframe.imagelength
                # add series
                aseries.levels.append(
                    TiffPageSeries(
                        ifds,
                        tuple(shape),
                        keyframe.dtype,
                        axes,
                        parent=self,
                        name=f'level {level + 1}',
                        kind='ome',
                    )
                )

        self.is_uniform = len(series) == 1 and len(series[0].levels) == 1

        return series

    def _series_mmstack(self) -> list[TiffPageSeries] | None:
        """Return series in Micro-Manager stack file(s)."""
        settings = self.micromanager_metadata
        if (
            settings is None
            or 'Summary' not in settings
            or 'IndexMap' not in settings
        ):
            return None

        pages: list[TiffPage | TiffFrame | None]
        page_count: int

        summary = settings['Summary']
        indexmap = settings['IndexMap']
        indexmap = indexmap[indexmap[:, 4].argsort()]

        if 'MicroManagerVersion' not in summary or 'Frames' not in summary:
            # TODO: handle MagellanStack?
            return None

        # determine CZTR shape from indexmap; TODO: is this necessary?
        indexmap_shape = (numpy.max(indexmap[:, :4], axis=0) + 1).tolist()
        indexmap_index = {'C': 0, 'Z': 1, 'T': 2, 'R': 3}

        # TODO: activate this?
        # if 'AxisOrder' in summary:
        #     axesorder = summary['AxisOrder']
        #     keys = {
        #         'channel': 'C',
        #         'z': 'Z',
        #         'slice': 'Z',
        #         'position': 'R',
        #         'time': 'T',
        #     }
        #     axes = ''.join(keys[ax] for ax in reversed(axesorder))

        axes = 'TR' if summary.get('TimeFirst', True) else 'RT'
        axes += 'ZC' if summary.get('SlicesFirst', True) else 'CZ'

        keys = {
            'C': 'Channels',
            'Z': 'Slices',
            'R': 'Positions',
            'T': 'Frames',
        }
        shape = tuple(
            max(
                indexmap_shape[indexmap_index[ax]],
                int(summary.get(keys[ax], 1)),
            )
            for ax in axes
        )
        size = product(shape)

        indexmap_order = tuple(indexmap_index[ax] for ax in axes)

        def add_file(tif: TiffFile, indexmap: NDArray[Any]) -> int:
            # add virtual TiffFrames to pages list
            page_count = 0
            offsets: list[int]
            offsets = indexmap[:, 4].tolist()
            indices = numpy.ravel_multi_index(
                # type: ignore[call-overload]
                indexmap[:, indexmap_order].T,
                shape,
            ).tolist()
            keyframe = tif.pages.first
            filesize = tif.filehandle.size - keyframe.databytecounts[0] - 162
            index: int
            offset: int
            for index, offset in zip(indices, offsets):
                if offset == keyframe.offset:
                    pages[index] = keyframe
                    page_count += 1
                    continue
                if 0 < offset <= filesize:
                    dataoffsets = (offset + 162,)
                    databytecounts = keyframe.databytecounts
                    page_count += 1
                else:
                    # assume file is truncated
                    dataoffsets = databytecounts = (0,)
                    offset = 0
                pages[index] = TiffFrame(
                    tif,
                    index=index,
                    offset=offset,
                    dataoffsets=dataoffsets,
                    databytecounts=databytecounts,
                    keyframe=keyframe,
                )
            return page_count

        multifile = size > indexmap.shape[0]
        if multifile:
            # get multifile prefix
            if not self.filehandle.is_file:
                logger().warning(
                    f'{self!r} MMStack multi-file series cannot be read from '
                    f'{self.filehandle._fh!r}'
                )
                multifile = False
            elif '_MMStack' not in self.filename:
                logger().warning(f'{self!r} MMStack file name is invalid')
                multifile = False
            elif 'Prefix' in summary:
                prefix = summary['Prefix']
                if not self.filename.startswith(prefix):
                    logger().warning(f'{self!r} MMStack file name is invalid')
                    multifile = False
            else:
                prefix = self.filename.split('_MMStack')[0]

        if multifile:
            # read other files
            pattern = os.path.join(
                self.filehandle.dirname, prefix + '_MMStack*.tif'
            )
            filenames = glob.glob(pattern)
            if len(filenames) == 1:
                multifile = False
            else:
                pages = [None] * size
                page_count = add_file(self, indexmap)
                for fname in filenames:
                    if self.filename == os.path.split(fname)[-1]:
                        continue
                    with TiffFile(fname) as tif:
                        indexmap = read_micromanager_metadata(
                            tif.filehandle, {'IndexMap'}
                        )['IndexMap']
                        indexmap = indexmap[indexmap[:, 4].argsort()]
                        page_count += add_file(tif, indexmap)

        if multifile:
            pass
        elif size > indexmap.shape[0]:
            # other files missing: squeeze shape
            old_shape = shape
            min_index = numpy.min(indexmap[:, :4], axis=0)
            max_index = numpy.max(indexmap[:, :4], axis=0)
            indexmap = indexmap.copy()
            indexmap[:, :4] -= min_index
            shape = tuple(
                j - i + 1
                for i, j in zip(min_index.tolist(), max_index.tolist())
            )
            shape = tuple(shape[i] for i in indexmap_order)
            size = product(shape)
            pages = [None] * size
            page_count = add_file(self, indexmap)
            logger().warning(
                f'{self!r} MMStack series is missing files. '
                f'Returning subset {shape!r} of {old_shape!r}'
            )
        else:
            # single file
            pages = [None] * size
            page_count = add_file(self, indexmap)

        if page_count != size:
            logger().warning(
                f'{self!r} MMStack is missing {size - page_count} pages.'
                ' Missing data are zeroed'
            )

        keyframe = self.pages.first
        return [
            TiffPageSeries(
                pages,
                shape=shape + keyframe.shape,
                dtype=keyframe.dtype,
                axes=axes + keyframe.axes,
                # axestiled=axestiled,
                # axesoverlap=axesoverlap,
                # coords=coords,
                parent=self,
                kind='mmstack',
                multifile=multifile,
                squeeze=True,
            )
        ]

    def _series_ndtiff(self) -> list[TiffPageSeries] | None:
        """Return series in NDTiff v2 and v3 files."""
        # TODO: implement fallback for missing index file, versions 0 and 1
        if not self.filehandle.is_file:
            logger().warning(
                f'{self!r} NDTiff.index not found for {self.filehandle._fh!r}'
            )
            return None

        indexfile = os.path.join(self.filehandle.dirname, 'NDTiff.index')
        if not os.path.exists(indexfile):
            logger().warning(f'{self!r} NDTiff.index not found')
            return None

        keyframes: dict[str, TiffPage] = {}
        shape: tuple[int, ...]
        dims: tuple[str, ...]
        page: TiffPage | TiffFrame
        pageindex = 0
        pixel_types = {
            0: ('uint8', 8),  # 8bit monochrome
            1: ('uint16', 16),  # 16bit monochrome
            2: ('uint8', 8),  # 8bit RGB
            3: ('uint16', 10),  # 10bit monochrome
            4: ('uint16', 12),  # 12bit monochrome
            5: ('uint16', 14),  # 14bit monochrome
            6: ('uint16', 11),  # 11bit monochrome
        }

        indices: dict[tuple[int, ...], TiffPage | TiffFrame] = {}
        categories: dict[str, dict[str, int]] = {}
        first = True

        for (
            axes_dict,
            filename,
            dataoffset,
            width,
            height,
            pixeltype,
            compression,
            metaoffset,
            metabytecount,
            metacompression,
        ) in read_ndtiff_index(indexfile):
            if filename in keyframes:
                # create virtual frame from index
                pageindex += 1  # TODO
                keyframe = keyframes[filename]
                page = TiffFrame(
                    keyframe.parent,
                    pageindex,
                    offset=None,  # virtual frame
                    keyframe=keyframe,
                    dataoffsets=(dataoffset,),
                    databytecounts=keyframe.databytecounts,
                )
                if page.shape[:2] != (height, width):
                    raise ValueError(
                        'NDTiff.index does not match TIFF shape '
                        f'{page.shape[:2]} != {(height, width)}'
                    )
                if compression != 0:
                    raise ValueError(
                        'NDTiff.index compression {compression} not supported'
                    )
                if page.compression != 1:
                    raise ValueError(
                        'NDTiff.index does not match TIFF compression '
                        f'{page.compression!r}'
                    )
                if pixeltype not in pixel_types:
                    raise ValueError(
                        f'NDTiff.index unknown pixel type {pixeltype}'
                    )
                dtype, _ = pixel_types[pixeltype]
                if page.dtype != dtype:
                    raise ValueError(
                        'NDTiff.index pixeltype does not match TIFF dtype '
                        f'{page.dtype} != {dtype}'
                    )
            elif filename == self.filename:
                # use first page as keyframe
                pageindex = 0
                page = self.pages.first
                keyframes[filename] = page
            else:
                # read keyframe from file
                pageindex = 0
                with TiffFile(
                    os.path.join(self.filehandle.dirname, filename)
                ) as tif:
                    page = tif.pages.first
                keyframes[filename] = page

            # replace string with integer indices
            index: int | str
            if first:
                for axis, index in axes_dict.items():
                    if isinstance(index, str):
                        categories[axis] = {index: 0}
                        axes_dict[axis] = 0
                first = False
            elif categories:
                for axis, values in categories.items():
                    index = axes_dict[axis]
                    assert isinstance(index, str)
                    if index not in values:
                        values[index] = max(values.values()) + 1
                    axes_dict[axis] = values[index]

            indices[tuple(axes_dict.values())] = page  # type: ignore[arg-type]
            dims = tuple(axes_dict.keys())

        # indices may be negative or missing
        indices_array = numpy.array(list(indices.keys()), dtype=numpy.int32)
        min_index = numpy.min(indices_array, axis=0).tolist()
        max_index = numpy.max(indices_array, axis=0).tolist()
        shape = tuple(j - i + 1 for i, j in zip(min_index, max_index))

        # change axes to match storage order
        order = order_axes(indices_array, squeeze=False)
        shape = tuple(shape[i] for i in order)
        dims = tuple(dims[i] for i in order)
        indices = {
            tuple(index[i] - min_index[i] for i in order): value
            for index, value in indices.items()
        }

        pages: list[TiffPage | TiffFrame | None] = []
        for idx in numpy.ndindex(shape):
            pages.append(indices.get(idx, None))

        keyframe = next(i for i in keyframes.values())
        shape += keyframe.shape
        dims += keyframe.dims
        axes = ''.join(TIFF.AXES_CODES.get(i.lower(), 'Q') for i in dims)

        # TODO: support tiled axes and overlap
        # meta: Any = self.micromanager_metadata
        # if meta is None:
        #     meta = {}
        # elif 'Summary' in meta:
        #     meta = meta['Summary']

        # # map axes column->x, row->y
        # axestiled: dict[int, int] = {}
        # axesoverlap: dict[int, int] = {}
        # if 'column' in dims:
        #     key = dims.index('column')
        #     axestiled[key] = keyframe.axes.index('X')
        #     axesoverlap[key] = meta.get('GridPixelOverlapX', 0)
        # if 'row' in dims:
        #     key = dims.index('row')
        #     axestiled[key] = keyframe.axes.index('Y')
        #     axesoverlap[key] = meta.get('GridPixelOverlapY', 0)

        # if all(i == 0 for i in axesoverlap.values()):
        #     axesoverlap = {}

        self.is_uniform = True
        return [
            TiffPageSeries(
                pages,
                shape=shape,
                dtype=keyframe.dtype,
                axes=axes,
                # axestiled=axestiled,
                # axesoverlap=axesoverlap,
                # coords=coords,
                parent=self,
                kind='ndtiff',
                multifile=len(keyframes) > 1,
                squeeze=True,
            )
        ]

    def _series_stk(self) -> list[TiffPageSeries] | None:
        """Return series in STK file."""
        meta = self.stk_metadata
        if meta is None:
            return None
        page = self.pages.first
        planes = meta['NumberPlanes']
        name = meta.get('Name', '')
        if planes == 1:
            shape = (1,) + page.shape
            axes = 'I' + page.axes
        elif numpy.all(meta['ZDistance'] != 0):
            shape = (planes,) + page.shape
            axes = 'Z' + page.axes
        elif numpy.all(numpy.diff(meta['TimeCreated']) != 0):
            shape = (planes,) + page.shape
            axes = 'T' + page.axes
        else:
            # TODO: determine other/combinations of dimensions
            shape = (planes,) + page.shape
            axes = 'I' + page.axes
        self.is_uniform = True
        series = TiffPageSeries(
            [page],
            shape,
            page.dtype,
            axes,
            name=name,
            truncated=planes > 1,
            kind='stk',
        )
        return [series]

    def _series_lsm(self) -> list[TiffPageSeries] | None:
        """Return main and thumbnail series in LSM file."""
        lsmi = self.lsm_metadata
        if lsmi is None:
            return None
        axes = TIFF.CZ_LSMINFO_SCANTYPE[lsmi['ScanType']]
        if self.pages.first.planarconfig == 1:
            axes = axes.replace('C', '').replace('X', 'XC')
        elif self.pages.first.planarconfig == 2:
            # keep axis for `get_shape(False)`
            pass
        elif self.pages.first.samplesperpixel == 1:
            axes = axes.replace('C', '')
        if lsmi.get('DimensionP', 0) > 0:
            axes = 'P' + axes
        if lsmi.get('DimensionM', 0) > 0:
            axes = 'M' + axes
        shape = tuple(int(lsmi[TIFF.CZ_LSMINFO_DIMENSIONS[i]]) for i in axes)

        name = lsmi.get('Name', '')
        pages = self.pages._getlist(slice(0, None, 2), validate=False)
        dtype = pages[0].dtype
        series = [
            TiffPageSeries(pages, shape, dtype, axes, name=name, kind='lsm')
        ]

        page = cast(TiffPage, self.pages[1])
        if page.is_reduced:
            pages = self.pages._getlist(slice(1, None, 2), validate=False)
            dtype = page.dtype
            cp = 1
            i = 0
            while cp < len(pages) and i < len(shape) - 2:
                cp *= shape[i]
                i += 1
            shape = shape[:i] + page.shape
            axes = axes[:i] + page.axes
            series.append(
                TiffPageSeries(
                    pages, shape, dtype, axes, name=name, kind='lsm'
                )
            )

        self.is_uniform = False
        return series

    def _lsm_load_pages(self) -> None:
        """Read and fix all pages from LSM file."""
        # cache all pages to preserve corrected values
        pages = self.pages
        pages.cache = True
        pages.useframes = True
        # use first and second page as keyframes
        pages.set_keyframe(1)
        pages.set_keyframe(0)
        # load remaining pages as frames
        pages._load(None)
        # fix offsets and bytecounts first
        # TODO: fix multiple conversions between lists and tuples
        self._lsm_fix_strip_offsets()
        self._lsm_fix_strip_bytecounts()
        # assign keyframes for data and thumbnail series
        keyframe = self.pages.first
        for page in pages._pages[::2]:
            page.keyframe = keyframe  # type: ignore[union-attr]
        keyframe = cast(TiffPage, pages[1])
        for page in pages._pages[1::2]:
            page.keyframe = keyframe  # type: ignore[union-attr]

    def _lsm_fix_strip_offsets(self) -> None:
        """Unwrap strip offsets for LSM files greater than 4 GB.

        Each series and position require separate unwrapping (undocumented).

        """
        if self.filehandle.size < 2**32:
            return

        indices: NDArray[Any]
        pages = self.pages
        npages = len(pages)
        series = self.series[0]
        axes = series.axes

        # find positions
        positions = 1
        for i in 0, 1:
            if series.axes[i] in 'PM':
                positions *= series.shape[i]

        # make time axis first
        if positions > 1:
            ntimes = 0
            for i in 1, 2:
                if axes[i] == 'T':
                    ntimes = series.shape[i]
                    break
            if ntimes:
                div, mod = divmod(npages, 2 * positions * ntimes)
                if mod != 0:
                    raise RuntimeError('mod != 0')
                shape = (positions, ntimes, div, 2)
                indices = numpy.arange(product(shape)).reshape(shape)
                indices = numpy.moveaxis(indices, 1, 0)
            else:
                indices = numpy.arange(npages).reshape(-1, 2)
        else:
            indices = numpy.arange(npages).reshape(-1, 2)

        # images of reduced page might be stored first
        if pages[0].dataoffsets[0] > pages[1].dataoffsets[0]:
            indices = indices[..., ::-1]

        # unwrap offsets
        wrap = 0
        previousoffset = 0
        for npi in indices.flat:
            page = pages[int(npi)]
            dataoffsets = []
            if all(i <= 0 for i in page.dataoffsets):
                logger().warning(
                    f'{self!r} LSM file incompletely written at {page}'
                )
                break
            for currentoffset in page.dataoffsets:
                if currentoffset < previousoffset:
                    wrap += 2**32
                dataoffsets.append(currentoffset + wrap)
                previousoffset = currentoffset
            page.dataoffsets = tuple(dataoffsets)

    def _lsm_fix_strip_bytecounts(self) -> None:
        """Set databytecounts to size of compressed data.

        The StripByteCounts tag in LSM files contains the number of bytes
        for the uncompressed data.

        """
        if self.pages.first.compression == 1:
            return
        # sort pages by first strip offset
        pages = sorted(self.pages, key=lambda p: p.dataoffsets[0])
        npages = len(pages) - 1
        for i, page in enumerate(pages):
            if page.index % 2:
                continue
            offsets = page.dataoffsets
            bytecounts = page.databytecounts
            if i < npages:
                lastoffset = pages[i + 1].dataoffsets[0]
            else:
                # LZW compressed strips might be longer than uncompressed
                lastoffset = min(
                    offsets[-1] + 2 * bytecounts[-1], self._fh.size
                )
            bytecount_list = list(bytecounts)
            for j in range(len(bytecounts) - 1):
                bytecount_list[j] = offsets[j + 1] - offsets[j]
            bytecount_list[-1] = lastoffset - offsets[-1]
            page.databytecounts = tuple(bytecount_list)

    def _ndpi_load_pages(self) -> None:
        """Read and fix pages from NDPI slide file if CaptureMode > 6.

        If the value of the CaptureMode tag is greater than 6, change the
        attributes of TiffPage instances that are part of the pyramid to
        match 16-bit grayscale data. TiffTag values are not corrected.

        """
        pages = self.pages
        capturemode = self.pages.first.tags.valueof(65441)
        if capturemode is None or capturemode < 6:
            return

        pages.cache = True
        pages.useframes = False
        pages._load()

        for page in pages:
            assert isinstance(page, TiffPage)
            mag = page.tags.valueof(65421)
            if mag is None or mag > 0:
                page.photometric = PHOTOMETRIC.MINISBLACK
                page.sampleformat = SAMPLEFORMAT.UINT
                page.samplesperpixel = 1
                page.bitspersample = 16
                page.dtype = page._dtype = numpy.dtype(numpy.uint16)
                if page.shaped[-1] > 1:
                    page.axes = page.axes[:-1]
                    page.shape = page.shape[:-1]
                    page.shaped = page.shaped[:-1] + (1,)

    def __getattr__(self, name: str, /) -> bool:
        """Return `is_flag` attributes from first page."""
        if name[3:] in TIFF.PAGE_FLAGS:
            if not self.pages:
                return False
            value = bool(getattr(self.pages.first, name))
            setattr(self, name, value)
            return value
        raise AttributeError(
            f'{self.__class__.__name__!r} object has no attribute {name!r}'
        )

    def __enter__(self) -> TiffFile:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        return f'<tifffile.TiffFile {snipstr(self._fh.name, 32)!r}>'

    def __str__(self) -> str:
        return self._str()

    def _str(self, detail: int = 0, width: int = 79) -> str:
        """Return string containing information about TiffFile.

        The `detail` parameter specifies the level of detail returned:

        0: file only.
        1: all series, first page of series and its tags.
        2: large tag values and file metadata.
        3: all pages.

        """
        info_list = [
            "TiffFile '{}'",
            format_size(self._fh.size),
            (
                ''
                if byteorder_isnative(self.byteorder)
                else {'<': 'little-endian', '>': 'big-endian'}[self.byteorder]
            ),
        ]
        if self.is_bigtiff:
            info_list.append('BigTiff')
        if len(self.pages) > 1:
            info_list.append(f'{len(self.pages)} Pages')
        if len(self.series) > 1:
            info_list.append(f'{len(self.series)} Series')
        if len(self._files) > 1:
            info_list.append(f'{len(self._files)} Files')
        flags = self.flags
        if 'uniform' in flags and len(self.pages) == 1:
            flags.discard('uniform')
        info_list.append('|'.join(f.lower() for f in sorted(flags)))
        info = '  '.join(info_list)
        info = info.replace('    ', '  ').replace('   ', '  ')
        info = info.format(
            snipstr(self._fh.name, max(12, width + 2 - len(info)))
        )
        if detail <= 0:
            return info
        info_list = [info]
        info_list.append('\n'.join(str(s) for s in self.series))
        if detail >= 3:
            for page in self.pages:
                if page is None:
                    continue
                info_list.append(page._str(detail=detail, width=width))
                if page.pages is not None:
                    for subifd in page.pages:
                        info_list.append(
                            subifd._str(detail=detail, width=width)
                        )
        elif self.series:
            info_list.extend(
                s.keyframe._str(detail=detail, width=width)
                for s in self.series
                if not s.keyframe.parent.filehandle.closed  # avoid warning
            )
        elif self.pages:  # and self.pages.first:
            info_list.append(self.pages.first._str(detail=detail, width=width))
        if detail >= 2:
            for name in sorted(self.flags):
                if hasattr(self, name + '_metadata'):
                    m = getattr(self, name + '_metadata')
                    if m:
                        info_list.append(
                            f'{name.upper()}_METADATA\n'
                            f'{pformat(m, width=width, height=detail * 24)}'
                        )
        return '\n\n'.join(info_list).replace('\n\n\n', '\n\n')

    @cached_property
    def flags(self) -> set[str]:
        """Set of file flags (a potentially expensive operation)."""
        return {
            name.lower()
            for name in TIFF.FILE_FLAGS
            if getattr(self, 'is_' + name)
        }

    @cached_property
    def is_uniform(self) -> bool:
        """File contains uniform series of pages."""
        # the hashes of IFDs 0, 7, and -1 are the same
        pages = self.pages
        try:
            page = self.pages.first
        except IndexError:
            return False
        if page.subifds:
            return False
        if page.is_scanimage or page.is_nih:
            return True
        i = 0
        useframes = pages.useframes
        try:
            pages.useframes = False
            h = page.hash
            for i in (1, 7, -1):
                if pages[i].aspage().hash != h:
                    return False
        except IndexError:
            return i == 1  # single page TIFF is uniform
        finally:
            pages.useframes = useframes
        return True

    @property
    def is_appendable(self) -> bool:
        """Pages can be appended to file without corrupting."""
        # TODO: check other formats
        return not (
            self.is_ome
            or self.is_lsm
            or self.is_stk
            or self.is_imagej
            or self.is_fluoview
            or self.is_micromanager
        )

    @property
    def is_bigtiff(self) -> bool:
        """File has BigTIFF format."""
        return self.tiff.is_bigtiff

    @cached_property
    def is_ndtiff(self) -> bool:
        """File has NDTiff format."""
        # file should be accompanied by NDTiff.index
        meta = self.micromanager_metadata
        if meta is not None and meta.get('MajorVersion', 0) >= 2:
            self.is_uniform = True
            return True
        return False

    @cached_property
    def is_mmstack(self) -> bool:
        """File has Micro-Manager stack format."""
        meta = self.micromanager_metadata
        if (
            meta is not None
            and 'Summary' in meta
            and 'IndexMap' in meta
            and meta.get('MajorVersion', 1) == 0
            # and 'MagellanStack' not in self.filename:
        ):
            self.is_uniform = True
            return True
        return False

    @cached_property
    def is_mdgel(self) -> bool:
        """File has MD Gel format."""
        # side effect: add second page, if exists, to cache
        try:
            ismdgel = (
                self.pages.first.is_mdgel
                or self.pages.get(1, cache=True).is_mdgel
            )
            if ismdgel:
                self.is_uniform = False
            return ismdgel
        except IndexError:
            return False

    @property
    def is_sis(self) -> bool:
        """File is Olympus SIS format."""
        try:
            return (
                self.pages.first.is_sis
                and not self.filename.lower().endswith('.vsi')
            )
        except IndexError:
            return False

    @cached_property
    def shaped_metadata(self) -> tuple[dict[str, Any], ...] | None:
        """Tifffile metadata from JSON formatted ImageDescription tags."""
        if not self.is_shaped:
            return None
        result = []
        for s in self.series:
            if s.kind.lower() != 'shaped':
                continue
            page = s.pages[0]
            if (
                not isinstance(page, TiffPage)
                or page.shaped_description is None
            ):
                continue
            result.append(shaped_description_metadata(page.shaped_description))
        return tuple(result)

    @property
    def ome_metadata(self) -> str | None:
        """OME XML metadata from ImageDescription tag."""
        if not self.is_ome:
            return None
        # return xml2dict(self.pages.first.description)['OME']
        if self._omexml:
            return self._omexml
        return self.pages.first.description

    @property
    def scn_metadata(self) -> str | None:
        """Leica SCN XML metadata from ImageDescription tag."""
        if not self.is_scn:
            return None
        return self.pages.first.description

    @property
    def philips_metadata(self) -> str | None:
        """Philips DP XML metadata from ImageDescription tag."""
        if not self.is_philips:
            return None
        return self.pages.first.description

    @property
    def indica_metadata(self) -> str | None:
        """IndicaLabs XML metadata from ImageDescription tag."""
        if not self.is_indica:
            return None
        return self.pages.first.description

    @property
    def avs_metadata(self) -> str | None:
        """Argos AVS XML metadata from tag 65000."""
        if not self.is_avs:
            return None
        return self.pages.first.tags.valueof(65000)

    @property
    def lsm_metadata(self) -> dict[str, Any] | None:
        """LSM metadata from CZ_LSMINFO tag."""
        if not self.is_lsm:
            return None
        return self.pages.first.tags.valueof(34412)  # CZ_LSMINFO

    @cached_property
    def stk_metadata(self) -> dict[str, Any] | None:
        """STK metadata from UIC tags."""
        if not self.is_stk:
            return None
        page = self.pages.first
        tags = page.tags
        result: dict[str, Any] = {}
        if page.description:
            result['PlaneDescriptions'] = page.description.split('\x00')
        tag = tags.get(33629)  # UIC2tag
        result['NumberPlanes'] = 1 if tag is None else tag.count
        value = tags.valueof(33628)  # UIC1tag
        if value is not None:
            result.update(value)
        value = tags.valueof(33630)  # UIC3tag
        if value is not None:
            result.update(value)  # wavelengths
        value = tags.valueof(33631)  # UIC4tag
        if value is not None:
            result.update(value)  # override UIC1 tags
        uic2tag = tags.valueof(33629)
        if uic2tag is not None:
            result['ZDistance'] = uic2tag['ZDistance']
            result['TimeCreated'] = uic2tag['TimeCreated']
            result['TimeModified'] = uic2tag['TimeModified']
            for key in ('Created', 'Modified'):
                try:
                    result['Datetime' + key] = numpy.array(
                        [
                            julian_datetime(*dt)
                            for dt in zip(
                                uic2tag['Date' + key], uic2tag['Time' + key]
                            )
                        ],
                        dtype='datetime64[ns]',
                    )
                except Exception as exc:
                    result['Datetime' + key] = None
                    logger().warning(
                        f'{self!r} STK Datetime{key} raised {exc!r:.128}'
                    )
        return result

    @cached_property
    def imagej_metadata(self) -> dict[str, Any] | None:
        """ImageJ metadata from ImageDescription and IJMetadata tags."""
        if not self.is_imagej:
            return None
        page = self.pages.first
        if page.imagej_description is None:
            return None
        result = imagej_description_metadata(page.imagej_description)
        value = page.tags.valueof(50839)  # IJMetadata
        if value is not None:
            try:
                result.update(value)
            except Exception:
                pass
        return result

    @cached_property
    def fluoview_metadata(self) -> dict[str, Any] | None:
        """FluoView metadata from MM_Header and MM_Stamp tags."""
        if not self.is_fluoview:
            return None
        result = {}
        page = self.pages.first
        value = page.tags.valueof(34361)  # MM_Header
        if value is not None:
            result.update(value)
        # TODO: read stamps from all pages
        value = page.tags.valueof(34362)  # MM_Stamp
        if value is not None:
            result['Stamp'] = value
        # skip parsing image description; not reliable
        # try:
        #     t = fluoview_description_metadata(page.image_description)
        #     if t is not None:
        #         result['ImageDescription'] = t
        # except Exception as exc:
        #     logger().warning(
        #         f'{self!r} <fluoview_description_metadata> '
        #         f'raised {exc!r:.128}'
        #     )
        return result

    @property
    def nih_metadata(self) -> dict[str, Any] | None:
        """NIHImage metadata from NIHImageHeader tag."""
        if not self.is_nih:
            return None
        return self.pages.first.tags.valueof(43314)  # NIHImageHeader

    @property
    def fei_metadata(self) -> dict[str, Any] | None:
        """FEI metadata from SFEG or HELIOS tags."""
        if not self.is_fei:
            return None
        tags = self.pages.first.tags
        result = {}
        try:
            result.update(tags.valueof(34680))  # FEI_SFEG
        except Exception:
            pass
        try:
            result.update(tags.valueof(34682))  # FEI_HELIOS
        except Exception:
            pass
        return result

    @property
    def sem_metadata(self) -> dict[str, Any] | None:
        """SEM metadata from CZ_SEM tag."""
        if not self.is_sem:
            return None
        return self.pages.first.tags.valueof(34118)

    @property
    def sis_metadata(self) -> dict[str, Any] | None:
        """Olympus SIS metadata from OlympusSIS and OlympusINI tags."""
        if not self.pages.first.is_sis:
            return None
        tags = self.pages.first.tags
        result = {}
        try:
            result.update(tags.valueof(33471))  # OlympusINI
        except Exception:
            pass
        try:
            result.update(tags.valueof(33560))  # OlympusSIS
        except Exception:
            pass
        return result if result else None

    @cached_property
    def mdgel_metadata(self) -> dict[str, Any] | None:
        """MD-GEL metadata from MDFileTag tags."""
        if not self.is_mdgel:
            return None
        if 33445 in self.pages.first.tags:
            tags = self.pages.first.tags
        else:
            page = cast(TiffPage, self.pages[1])
            if 33445 in page.tags:
                tags = page.tags
            else:
                return None
        result = {}
        for code in range(33445, 33453):
            if code not in tags:
                continue
            name = TIFF.TAGS[code]
            result[name[2:]] = tags.valueof(code)
        return result

    @property
    def andor_metadata(self) -> dict[str, Any] | None:
        """Andor metadata from Andor tags."""
        return self.pages.first.andor_tags

    @property
    def epics_metadata(self) -> dict[str, Any] | None:
        """EPICS metadata from areaDetector tags."""
        return self.pages.first.epics_tags

    @property
    def tvips_metadata(self) -> dict[str, Any] | None:
        """TVIPS metadata from tag."""
        if not self.is_tvips:
            return None
        return self.pages.first.tags.valueof(37706)

    @cached_property
    def metaseries_metadata(self) -> dict[str, Any] | None:
        """MetaSeries metadata from ImageDescription tag of first tag."""
        # TODO: remove this? It is a per page property
        if not self.is_metaseries:
            return None
        return metaseries_description_metadata(self.pages.first.description)

    @cached_property
    def pilatus_metadata(self) -> dict[str, Any] | None:
        """Pilatus metadata from ImageDescription tag."""
        if not self.is_pilatus:
            return None
        return pilatus_description_metadata(self.pages.first.description)

    @cached_property
    def micromanager_metadata(self) -> dict[str, Any] | None:
        """Non-TIFF Micro-Manager metadata."""
        if not self.is_micromanager:
            return None
        return read_micromanager_metadata(self._fh)

    @cached_property
    def gdal_structural_metadata(self) -> dict[str, Any] | None:
        """Non-TIFF GDAL structural metadata."""
        return read_gdal_structural_metadata(self._fh)

    @cached_property
    def scanimage_metadata(self) -> dict[str, Any] | None:
        """ScanImage non-varying frame and ROI metadata.

        The returned dict may contain 'FrameData', 'RoiGroups', and 'version'
        keys.

        Varying frame data can be found in the ImageDescription tags.

        """
        if not self.is_scanimage:
            return None
        result: dict[str, Any] = {}
        try:
            framedata, roidata, version = read_scanimage_metadata(self._fh)
            result['version'] = version
            result['FrameData'] = framedata
            result.update(roidata)
        except ValueError:
            pass
        return result

    @property
    def geotiff_metadata(self) -> dict[str, Any] | None:
        """GeoTIFF metadata from tags."""
        if not self.is_geotiff:
            return None
        return self.pages.first.geotiff_tags

    @property
    def gdal_metadata(self) -> dict[str, Any] | None:
        """GDAL XML metadata from GDAL_METADATA tag."""
        if not self.is_gdal:
            return None
        return self.pages.first.tags.valueof(42112)

    @cached_property
    def astrotiff_metadata(self) -> dict[str, Any] | None:
        """AstroTIFF metadata from ImageDescription tag."""
        if not self.is_astrotiff:
            return None
        return astrotiff_description_metadata(self.pages.first.description)

    @cached_property
    def streak_metadata(self) -> dict[str, Any] | None:
        """Hamamatsu streak metadata from ImageDescription tag."""
        if not self.is_streak:
            return None
        return streak_description_metadata(
            self.pages.first.description, self.filehandle
        )

    @property
    def eer_metadata(self) -> str | None:
        """EER AcquisitionMetadata XML from tag 65001."""
        if not self.is_eer:
            return None
        value = self.pages.first.tags.valueof(65001)
        return None if value is None else value.decode()


@final
class TiffFormat:
    """TIFF format properties."""

    __slots__ = (
        'version',
        'byteorder',
        'offsetsize',
        'offsetformat',
        'tagnosize',
        'tagnoformat',
        'tagsize',
        'tagformat1',
        'tagformat2',
        'tagoffsetthreshold',
        '_hash',
    )

    version: int
    """Version of TIFF header."""

    byteorder: Literal['>', '<']
    """Byteorder of TIFF header."""

    offsetsize: int
    """Size of offsets."""

    offsetformat: str
    """Struct format for offset values."""

    tagnosize: int
    """Size of `tagnoformat`."""

    tagnoformat: str
    """Struct format for number of TIFF tags."""

    tagsize: int
    """Size of `tagformat1` and `tagformat2`."""

    tagformat1: str
    """Struct format for code and dtype of TIFF tag."""

    tagformat2: str
    """Struct format for count and value of TIFF tag."""

    tagoffsetthreshold: int
    """Size of inline tag values."""

    _hash: int

    def __init__(
        self,
        version: int,
        byteorder: Literal['>', '<'],
        offsetsize: int,
        offsetformat: str,
        tagnosize: int,
        tagnoformat: str,
        tagsize: int,
        tagformat1: str,
        tagformat2: str,
        tagoffsetthreshold: int,
    ) -> None:
        self.version = version
        self.byteorder = byteorder
        self.offsetsize = offsetsize
        self.offsetformat = offsetformat
        self.tagnosize = tagnosize
        self.tagnoformat = tagnoformat
        self.tagsize = tagsize
        self.tagformat1 = tagformat1
        self.tagformat2 = tagformat2
        self.tagoffsetthreshold = tagoffsetthreshold
        self._hash = hash((version, byteorder, offsetsize))

    @property
    def is_bigtiff(self) -> bool:
        """Format is 64-bit BigTIFF."""
        return self.version == 43

    @property
    def is_ndpi(self) -> bool:
        """Format is 32-bit TIFF with 64-bit offsets used by NDPI."""
        return self.version == 42 and self.offsetsize == 8

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        bits = '32' if self.version == 42 else '64'
        endian = 'little' if self.byteorder == '<' else 'big'
        ndpi = ' with 64-bit offsets' if self.is_ndpi else ''
        return f'<tifffile.TiffFormat {bits}-bit {endian}-endian{ndpi}>'

    def __str__(self) -> str:
        return indent(
            repr(self),
            *(
                f'{attr}: {getattr(self, attr)!r}'
                for attr in TiffFormat.__slots__
            ),
        )


@final
class TiffPage:
    """TIFF image file directory (IFD).

    TiffPage instances are not thread-safe. All attributes are read-only.

    Parameters:
        parent:
            TiffFile instance to read page from.
            The file handle position must be at an offset to an IFD structure.
        index:
            Index of page in IFD tree.
        keyframe:
            Not used.

    Raises:
        TiffFileError: Invalid TIFF structure.

    """

    # instance attributes

    tags: TiffTags
    """Tags belonging to page."""

    parent: TiffFile
    """TiffFile instance page belongs to."""

    offset: int
    """Position of page in file."""

    shape: tuple[int, ...]
    """Shape of image array in page."""

    dtype: numpy.dtype[Any] | None
    """Data type of image array in page."""

    shaped: tuple[int, int, int, int, int]
    """Normalized 5-dimensional shape of image array in page:

        0. separate samplesperpixel or 1.
        1. imagedepth or 1.
        2. imagelength.
        3. imagewidth.
        4. contig samplesperpixel or 1.

    """

    axes: str
    """Character codes for dimensions in image array:
    'S' sample, 'X' width, 'Y' length, 'Z' depth.
    """

    dataoffsets: tuple[int, ...]
    """Positions of strips or tiles in file."""

    databytecounts: tuple[int, ...]
    """Size of strips or tiles in file."""

    _dtype: numpy.dtype[Any] | None
    _index: tuple[int, ...]  # index of page in IFD tree

    # default properties; might be updated from tags

    subfiletype: int = 0
    """:py:class:`FILETYPE` kind of image."""

    imagewidth: int = 0
    """Number of columns (pixels per row) in image."""

    imagelength: int = 0
    """Number of rows in image."""

    imagedepth: int = 1
    """Number of Z slices in image."""

    tilewidth: int = 0
    """Number of columns in each tile."""

    tilelength: int = 0
    """Number of rows in each tile."""

    tiledepth: int = 1
    """Number of Z slices in each tile."""

    samplesperpixel: int = 1
    """Number of components per pixel."""

    bitspersample: int = 1
    """Number of bits per pixel component."""

    sampleformat: int = 1
    """:py:class:`SAMPLEFORMAT` type of pixel components."""

    rowsperstrip: int = 2**32 - 1
    """Number of rows per strip."""

    compression: int = 1
    """:py:class:`COMPRESSION` scheme used on image data."""

    planarconfig: int = 1
    """:py:class:`PLANARCONFIG` type of storage of components in pixel."""

    fillorder: int = 1
    """Logical order of bits within byte of image data."""

    photometric: int = 0
    """:py:class:`PHOTOMETRIC` color space of image."""

    predictor: int = 1
    """:py:class:`PREDICTOR` applied to image data before compression."""

    extrasamples: tuple[int, ...] = ()
    """:py:class:`EXTRASAMPLE` interpretation of extra components in pixel."""

    subsampling: tuple[int, int] | None = None
    """Subsampling factors used for chrominance components."""

    subifds: tuple[int, ...] | None = None
    """Positions of SubIFDs in file."""

    jpegtables: bytes | None = None
    """JPEG quantization and Huffman tables."""

    jpegheader: bytes | None = None
    """JPEG header for NDPI."""

    software: str = ''
    """Software used to create image."""

    description: str = ''
    """Subject of image."""

    description1: str = ''
    """Value of second ImageDescription tag."""

    nodata: int | float = 0
    """Value used for missing data. The value of the GDAL_NODATA tag or 0."""

    def __init__(
        self,
        parent: TiffFile,
        /,
        index: int | Sequence[int],
        *,
        keyframe: TiffPage | None = None,
    ) -> None:
        tag: TiffTag | None
        tiff = parent.tiff

        self.parent = parent
        self.shape = ()
        self.shaped = (0, 0, 0, 0, 0)
        self.dtype = self._dtype = None
        self.axes = ''
        self.tags = tags = TiffTags()
        self.dataoffsets = ()
        self.databytecounts = ()
        if isinstance(index, int):
            self._index = (index,)
        else:
            self._index = tuple(index)

        # read IFD structure and its tags from file
        fh = parent.filehandle
        self.offset = fh.tell()  # offset to this IFD
        try:
            tagno: int = struct.unpack(
                tiff.tagnoformat, fh.read(tiff.tagnosize)
            )[0]
            if tagno > 4096:
                raise ValueError(f'suspicious number of tags {tagno}')
        except Exception as exc:
            raise TiffFileError(f'corrupted tag list @{self.offset}') from exc

        tagoffset = self.offset + tiff.tagnosize  # fh.tell()
        tagsize = tagsize_ = tiff.tagsize

        data = fh.read(tagsize * tagno)
        if len(data) != tagsize * tagno:
            raise TiffFileError('corrupted IFD structure')
        if tiff.is_ndpi:
            # patch offsets/values for 64-bit NDPI file
            tagsize = 16
            fh.seek(8, os.SEEK_CUR)
            ext = fh.read(4 * tagno)  # high bits
            data = b''.join(
                data[i * 12 : i * 12 + 12] + ext[i * 4 : i * 4 + 4]
                for i in range(tagno)
            )

        tagindex = -tagsize
        for i in range(tagno):
            tagindex += tagsize
            tagdata = data[tagindex : tagindex + tagsize]
            try:
                tag = TiffTag.fromfile(
                    parent, offset=tagoffset + i * tagsize_, header=tagdata
                )
            except TiffFileError as exc:
                logger().error(f'<TiffTag.fromfile> raised {exc!r:.128}')
                continue
            tags.add(tag)

        if not tags:
            return  # found in FIBICS

        for code, name in TIFF.TAG_ATTRIBUTES.items():
            value = tags.valueof(code)
            if value is None:
                continue
            if code in {270, 305} and not isinstance(value, str):
                # wrong string type for software or description
                continue
            setattr(self, name, value)

        value = tags.valueof(270, index=1)
        if isinstance(value, str):
            self.description1 = value

        if self.subfiletype == 0:
            value = tags.valueof(255)  # SubfileType
            if value == 2:
                self.subfiletype = 0b1  # reduced image
            elif value == 3:
                self.subfiletype = 0b10  # multi-page
        elif not isinstance(self.subfiletype, int):
            # files created by IDEAS
            logger().warning(f'{self!r} invalid {self.subfiletype=}')
            self.subfiletype = 0

        # consolidate private tags; remove them from self.tags
        # if self.is_andor:
        #     self.andor_tags
        # elif self.is_epics:
        #     self.epics_tags
        # elif self.is_ndpi:
        #     self.ndpi_tags
        # if self.is_sis and 34853 in tags:
        #     # TODO: cannot change tag.name
        #     tags[34853].name = 'OlympusSIS2'

        # dataoffsets and databytecounts
        # TileOffsets
        self.dataoffsets = tags.valueof(324)
        if self.dataoffsets is None:
            # StripOffsets
            self.dataoffsets = tags.valueof(273)
            if self.dataoffsets is None:
                # JPEGInterchangeFormat et al.
                self.dataoffsets = tags.valueof(513)
                if self.dataoffsets is None:
                    self.dataoffsets = ()
                    logger().error(f'{self!r} missing data offset tag')
        # TileByteCounts
        self.databytecounts = tags.valueof(325)
        if self.databytecounts is None:
            # StripByteCounts
            self.databytecounts = tags.valueof(279)
            if self.databytecounts is None:
                # JPEGInterchangeFormatLength et al.
                self.databytecounts = tags.valueof(514)

        if (
            self.imagewidth == 0
            and self.imagelength == 0
            and self.dataoffsets
            and self.databytecounts
        ):
            # dimensions may be missing in some RAW formats
            # read dimensions from assumed JPEG encoded segment
            try:
                fh.seek(self.dataoffsets[0])
                (
                    precision,
                    imagelength,
                    imagewidth,
                    samplesperpixel,
                ) = jpeg_shape(fh.read(min(self.databytecounts[0], 4096)))
            except Exception:
                pass
            else:
                self.imagelength = imagelength
                self.imagewidth = imagewidth
                self.samplesperpixel = samplesperpixel
                if 258 not in tags:
                    self.bitspersample = 8 if precision <= 8 else 16
                if 262 not in tags and samplesperpixel == 3:
                    self.photometric = PHOTOMETRIC.YCBCR
                if 259 not in tags:
                    self.compression = COMPRESSION.OJPEG
                if 278 not in tags:
                    self.rowsperstrip = imagelength

        elif self.compression == 6:
            # OJPEG hack. See libtiff v4.2.0 tif_dirread.c#L4082
            if 262 not in tags:
                # PhotometricInterpretation missing
                self.photometric = PHOTOMETRIC.YCBCR
            elif self.photometric == 2:
                # RGB -> YCbCr
                self.photometric = PHOTOMETRIC.YCBCR
            if 258 not in tags:
                # BitsPerSample missing
                self.bitspersample = 8
            if 277 not in tags:
                # SamplesPerPixel missing
                if self.photometric in {2, 6}:
                    self.samplesperpixel = 3
                elif self.photometric in {0, 1}:
                    self.samplesperpixel = 3

        elif self.is_lsm or (self.index != 0 and self.parent.is_lsm):
            # correct non standard LSM bitspersample tags
            tags[258]._fix_lsm_bitspersample()
            if self.compression == 1 and self.predictor != 1:
                # work around bug in LSM510 software
                self.predictor = PREDICTOR.NONE

        elif self.is_vista or (self.index != 0 and self.parent.is_vista):
            # ISS Vista writes wrong ImageDepth tag
            self.imagedepth = 1

        elif self.is_stk:
            # read UIC1tag again now that plane count is known
            tag = tags.get(33628)  # UIC1tag
            assert tag is not None
            fh.seek(tag.valueoffset)
            uic2tag = tags.get(33629)  # UIC2tag
            try:
                tag.value = read_uic1tag(
                    fh,
                    tiff.byteorder,
                    tag.dtype,
                    tag.count,
                    0,
                    planecount=uic2tag.count if uic2tag is not None else 1,
                )
            except Exception as exc:
                logger().warning(
                    f'{self!r} <tifffile.read_uic1tag> raised {exc!r:.128}'
                )

        tag = tags.get(50839)
        if tag is not None:
            # decode IJMetadata tag
            try:
                tag.value = imagej_metadata(
                    tag.value,
                    tags[50838].value,  # IJMetadataByteCounts
                    tiff.byteorder,
                )
            except Exception as exc:
                logger().warning(
                    f'{self!r} <tifffile.imagej_metadata> raised {exc!r:.128}'
                )

        # BitsPerSample
        value = tags.valueof(258)
        if value is not None:
            if self.bitspersample != 1:
                pass  # bitspersample was set by ojpeg hack
            elif tags[258].count == 1:
                self.bitspersample = int(value)
            else:
                # LSM might list more items than samplesperpixel
                value = value[: self.samplesperpixel]
                if any(v - value[0] for v in value):
                    self.bitspersample = value
                else:
                    self.bitspersample = int(value[0])

        # SampleFormat
        value = tags.valueof(339)
        if value is not None:
            if tags[339].count == 1:
                try:
                    self.sampleformat = SAMPLEFORMAT(value)
                except ValueError:
                    self.sampleformat = int(value)
            else:
                value = value[: self.samplesperpixel]
                if any(v - value[0] for v in value):
                    try:
                        self.sampleformat = SAMPLEFORMAT(value)
                    except ValueError:
                        self.sampleformat = int(value)
                else:
                    try:
                        self.sampleformat = SAMPLEFORMAT(value[0])
                    except ValueError:
                        self.sampleformat = int(value[0])
        elif self.bitspersample == 32 and (
            self.is_indica or (self.index != 0 and self.parent.is_indica)
        ):
            # IndicaLabsImageWriter does not write SampleFormat tag
            self.sampleformat = SAMPLEFORMAT.IEEEFP

        if 322 in tags:  # TileWidth
            self.rowsperstrip = 0
        elif 257 in tags:  # ImageLength
            if 278 not in tags or tags[278].count > 1:  # RowsPerStrip
                self.rowsperstrip = self.imagelength
            self.rowsperstrip = min(self.rowsperstrip, self.imagelength)
            # self.stripsperimage = int(math.floor(
            #    float(self.imagelength + self.rowsperstrip - 1) /
            #    self.rowsperstrip))

        # determine dtype
        dtypestr = TIFF.SAMPLE_DTYPES.get(
            (self.sampleformat, self.bitspersample), None
        )
        if dtypestr is not None:
            dtype = numpy.dtype(dtypestr)
        else:
            dtype = None
        self.dtype = self._dtype = dtype

        # determine shape of data
        imagelength = self.imagelength
        imagewidth = self.imagewidth
        imagedepth = self.imagedepth
        samplesperpixel = self.samplesperpixel

        if self.photometric == 2 or samplesperpixel > 1:  # PHOTOMETRIC.RGB
            if self.planarconfig == 1:
                self.shaped = (
                    1,
                    imagedepth,
                    imagelength,
                    imagewidth,
                    samplesperpixel,
                )
                if imagedepth == 1:
                    self.shape = (imagelength, imagewidth, samplesperpixel)
                    self.axes = 'YXS'
                else:
                    self.shape = (
                        imagedepth,
                        imagelength,
                        imagewidth,
                        samplesperpixel,
                    )
                    self.axes = 'ZYXS'
            else:
                self.shaped = (
                    samplesperpixel,
                    imagedepth,
                    imagelength,
                    imagewidth,
                    1,
                )
                if imagedepth == 1:
                    self.shape = (samplesperpixel, imagelength, imagewidth)
                    self.axes = 'SYX'
                else:
                    self.shape = (
                        samplesperpixel,
                        imagedepth,
                        imagelength,
                        imagewidth,
                    )
                    self.axes = 'SZYX'
        else:
            self.shaped = (1, imagedepth, imagelength, imagewidth, 1)
            if imagedepth == 1:
                self.shape = (imagelength, imagewidth)
                self.axes = 'YX'
            else:
                self.shape = (imagedepth, imagelength, imagewidth)
                self.axes = 'ZYX'

        if not self.databytecounts:
            self.databytecounts = (
                product(self.shape) * (self.bitspersample // 8),
            )
            if self.compression != 1:
                logger().error(f'{self!r} missing ByteCounts tag')

        if imagelength and self.rowsperstrip and not self.is_lsm:
            # fix incorrect number of strip bytecounts and offsets
            maxstrips = (
                int(
                    math.floor(imagelength + self.rowsperstrip - 1)
                    / self.rowsperstrip
                )
                * self.imagedepth
            )
            if self.planarconfig == 2:
                maxstrips *= self.samplesperpixel
            if maxstrips != len(self.databytecounts):
                logger().error(
                    f'{self!r} incorrect StripByteCounts count '
                    f'({len(self.databytecounts)} != {maxstrips})'
                )
                self.databytecounts = self.databytecounts[:maxstrips]
            if maxstrips != len(self.dataoffsets):
                logger().error(
                    f'{self!r} incorrect StripOffsets count '
                    f'({len(self.dataoffsets)} != {maxstrips})'
                )
                self.dataoffsets = self.dataoffsets[:maxstrips]

        value = tags.valueof(42113)  # GDAL_NODATA
        if value is not None and dtype is not None:
            try:
                pytype = type(dtype.type(0).item())
                value = value.replace(',', '.')  # comma decimal separator
                self.nodata = pytype(value)
                if not numpy.can_cast(
                    numpy.min_scalar_type(self.nodata), self.dtype
                ):
                    raise ValueError(
                        f'{self.nodata} is not castable to {self.dtype}'
                    )
            except Exception as exc:
                logger().warning(
                    f'{self!r} parsing GDAL_NODATA tag raised {exc!r:.128}'
                )
                self.nodata = 0

        mcustarts = tags.valueof(65426)
        if mcustarts is not None and self.is_ndpi:
            # use NDPI JPEG McuStarts as tile offsets
            mcustarts = mcustarts.astype(numpy.int64)
            high = tags.valueof(65432)
            if high is not None:
                # McuStartsHighBytes
                high = high.astype(numpy.uint64)
                high <<= 32
                mcustarts += high.astype(numpy.int64)
            fh.seek(self.dataoffsets[0])
            jpegheader = fh.read(mcustarts[0])
            try:
                (
                    self.tilelength,
                    self.tilewidth,
                    self.jpegheader,
                ) = ndpi_jpeg_tile(jpegheader)
            except ValueError as exc:
                logger().warning(
                    f'{self!r} <tifffile.ndpi_jpeg_tile> raised {exc!r:.128}'
                )
            else:
                # TODO: optimize tuple(ndarray.tolist())
                databytecounts = numpy.diff(
                    mcustarts, append=self.databytecounts[0]
                )
                self.databytecounts = tuple(databytecounts.tolist())
                mcustarts += self.dataoffsets[0]
                self.dataoffsets = tuple(mcustarts.tolist())

    @cached_property
    def decode(
        self,
    ) -> Callable[
        ...,
        tuple[
            NDArray[Any] | None,
            tuple[int, int, int, int, int],
            tuple[int, int, int, int],
        ],
    ]:
        """Return decoded segment, its shape, and indices in image.

        The decode function is implemented as a closure and has the following
        signature:

        Parameters:
            data (Union[bytes, None]):
                Encoded bytes of segment (strip or tile) or None for empty
                segments.
            index (int):
                Index of segment in Offsets and Bytecount tag values.
            jpegtables (Optional[bytes]):
                For JPEG compressed segments only, value of JPEGTables tag
                if any.

        Returns:
            - Decoded segment or None for empty segments.
            - Position of segment in image array of normalized shape
              (separate sample, depth, length, width, contig sample).
            - Shape of segment (depth, length, width, contig samples).
              The shape of strips depends on their linear index.

        Raises:
            ValueError or NotImplementedError:
                Decoding is not supported.
            TiffFileError:
                Invalid TIFF structure.

        """
        if self.hash in self.parent._parent._decoders:
            return self.parent._parent._decoders[self.hash]

        def cache(decode, /):
            self.parent._parent._decoders[self.hash] = decode
            return decode

        if self.dtype is None or self._dtype is None:

            def decode_raise_dtype(*args, **kwargs):
                raise ValueError(
                    'data type not supported '
                    f'(SampleFormat {self.sampleformat}, '
                    f'{self.bitspersample}-bit)'
                )

            return cache(decode_raise_dtype)

        if 0 in self.shaped:

            def decode_raise_empty(*args, **kwargs):
                raise ValueError('empty image')

            return cache(decode_raise_empty)

        try:
            if self.compression == 1:
                decompress = None
            else:
                decompress = TIFF.DECOMPRESSORS[self.compression]
            if (
                self.compression in {65000, 65001, 65002}
                and not self.parent.is_eer
            ):
                raise KeyError(self.compression)
        except KeyError as exc:

            def decode_raise_compression(*args, exc=str(exc)[1:-1], **kwargs):
                raise ValueError(f'{exc}')

            return cache(decode_raise_compression)

        try:
            if self.predictor == 1:
                unpredict = None
            else:
                unpredict = TIFF.UNPREDICTORS[self.predictor]
        except KeyError as exc:
            if self.compression in TIFF.IMAGE_COMPRESSIONS:
                logger().warning(
                    f'{self!r} ignoring predictor {self.predictor}'
                )
                unpredict = None

            else:

                def decode_raise_predictor(
                    *args, exc=str(exc)[1:-1], **kwargs
                ):
                    raise ValueError(f'{exc}')

                return cache(decode_raise_predictor)

        if self.tags.get(339) is not None:
            tag = self.tags[339]  # SampleFormat
            if tag.count != 1 and any(i - tag.value[0] for i in tag.value):

                def decode_raise_sampleformat(*args, **kwargs):
                    raise ValueError(
                        f'sample formats do not match {tag.value}'
                    )

                return cache(decode_raise_sampleformat)

        if self.is_subsampled and (
            self.compression not in {6, 7, 34892, 33007}
            or self.planarconfig == 2
        ):

            def decode_raise_subsampling(*args, **kwargs):
                raise NotImplementedError(
                    'chroma subsampling not supported without JPEG compression'
                )

            return cache(decode_raise_subsampling)

        if self.compression == 50001 and self.samplesperpixel == 4:
            # WebP segments may be missing all-opaque alpha channel
            def decompress_webp_rgba(data, out=None):
                return imagecodecs.webp_decode(data, hasalpha=True, out=out)

            decompress = decompress_webp_rgba

        # normalize segments shape to [depth, length, width, contig]
        if self.is_tiled:
            stshape = (
                self.tiledepth,
                self.tilelength,
                self.tilewidth,
                self.samplesperpixel if self.planarconfig == 1 else 1,
            )
        else:
            stshape = (
                1,
                self.rowsperstrip,
                self.imagewidth,
                self.samplesperpixel if self.planarconfig == 1 else 1,
            )

        stdepth, stlength, stwidth, samples = stshape
        _, imdepth, imlength, imwidth, samples = self.shaped

        if self.is_tiled:
            width = (imwidth + stwidth - 1) // stwidth
            length = (imlength + stlength - 1) // stlength
            depth = (imdepth + stdepth - 1) // stdepth

            def indices(
                segmentindex: int, /
            ) -> tuple[
                tuple[int, int, int, int, int], tuple[int, int, int, int]
            ]:
                # return indices and shape of tile in image array
                return (
                    (
                        segmentindex // (width * length * depth),
                        (segmentindex // (width * length)) % depth * stdepth,
                        (segmentindex // width) % length * stlength,
                        segmentindex % width * stwidth,
                        0,
                    ),
                    stshape,
                )

            def reshape(
                data: NDArray[Any],
                indices: tuple[int, int, int, int, int],
                shape: tuple[int, int, int, int],
                /,
            ) -> NDArray[Any]:
                # return reshaped tile or raise TiffFileError
                size = shape[0] * shape[1] * shape[2] * shape[3]
                if data.ndim == 1 and data.size > size:
                    # decompression / unpacking might return too many bytes
                    data = data[:size]
                if data.size == size:
                    # complete tile
                    # data might be non-contiguous; cannot reshape inplace
                    return data.reshape(shape)
                try:
                    # data fills remaining space
                    # found in JPEG/PNG compressed tiles
                    return data.reshape(
                        (
                            min(imdepth - indices[1], shape[0]),
                            min(imlength - indices[2], shape[1]),
                            min(imwidth - indices[3], shape[2]),
                            samples,
                        )
                    )
                except ValueError:
                    pass
                try:
                    # data fills remaining horizontal space
                    # found in tiled GeoTIFF
                    return data.reshape(
                        (
                            min(imdepth - indices[1], shape[0]),
                            min(imlength - indices[2], shape[1]),
                            shape[2],
                            samples,
                        )
                    )
                except ValueError:
                    pass
                raise TiffFileError(
                    f'corrupted tile @ {indices} cannot be reshaped from '
                    f'{data.shape} to {shape}'
                )

            def pad(
                data: NDArray[Any], shape: tuple[int, int, int, int], /
            ) -> tuple[NDArray[Any], tuple[int, int, int, int]]:
                # pad tile to shape
                if data.shape == shape:
                    return data, shape
                padwidth = [(0, i - j) for i, j in zip(shape, data.shape)]
                data = numpy.pad(data, padwidth, constant_values=self.nodata)
                return data, shape

            def pad_none(
                shape: tuple[int, int, int, int], /
            ) -> tuple[int, int, int, int]:
                # return shape of tile
                return shape

        else:
            # strips
            length = (imlength + stlength - 1) // stlength

            def indices(
                segmentindex: int, /
            ) -> tuple[
                tuple[int, int, int, int, int], tuple[int, int, int, int]
            ]:
                # return indices and shape of strip in image array
                indices = (
                    segmentindex // (length * imdepth),
                    (segmentindex // length) % imdepth * stdepth,
                    segmentindex % length * stlength,
                    0,
                    0,
                )
                shape = (
                    stdepth,
                    min(stlength, imlength - indices[2]),
                    stwidth,
                    samples,
                )
                return indices, shape

            def reshape(
                data: NDArray[Any],
                indices: tuple[int, int, int, int, int],
                shape: tuple[int, int, int, int],
                /,
            ) -> NDArray[Any]:
                # return reshaped strip or raise TiffFileError
                size = shape[0] * shape[1] * shape[2] * shape[3]
                if data.ndim == 1 and data.size > size:
                    # decompression / unpacking might return too many bytes
                    data = data[:size]
                if data.size == size:
                    # expected size
                    try:
                        data.shape = shape
                    except AttributeError:
                        # incompatible shape for in-place modification
                        # decoder returned non-contiguous array
                        data = data.reshape(shape)
                    return data
                datashape = data.shape
                try:
                    # too many rows?
                    data.shape = shape[0], -1, shape[2], shape[3]
                    data = data[:, : shape[1]]
                    data.shape = shape
                    return data
                except ValueError:
                    pass
                raise TiffFileError(
                    'corrupted strip cannot be reshaped from '
                    f'{datashape} to {shape}'
                )

            def pad(
                data: NDArray[Any], shape: tuple[int, int, int, int], /
            ) -> tuple[NDArray[Any], tuple[int, int, int, int]]:
                # pad strip length to rowsperstrip
                shape = (shape[0], stlength, shape[2], shape[3])
                if data.shape == shape:
                    return data, shape
                padwidth = [
                    (0, 0),
                    (0, stlength - data.shape[1]),
                    (0, 0),
                    (0, 0),
                ]
                data = numpy.pad(data, padwidth, constant_values=self.nodata)
                return data, shape

            def pad_none(
                shape: tuple[int, int, int, int], /
            ) -> tuple[int, int, int, int]:
                # return shape of strip
                return (shape[0], stlength, shape[2], shape[3])

        if self.compression in {6, 7, 34892, 33007}:
            # JPEG needs special handling
            if self.fillorder == 2:
                logger().debug(f'{self!r} disabling LSB2MSB for JPEG')
            if unpredict:
                logger().debug(f'{self!r} disabling predictor for JPEG')
            if 28672 in self.tags:  # SonyRawFileType
                logger().warning(
                    f'{self!r} SonyRawFileType might need additional '
                    'unpacking (see issue #95)'
                )

            colorspace, outcolorspace = jpeg_decode_colorspace(
                self.photometric,
                self.planarconfig,
                self.extrasamples,
                self.is_jfif,
            )

            def decode_jpeg(
                data: bytes | None,
                index: int,
                /,
                *,
                jpegtables: bytes | None = None,
                jpegheader: bytes | None = None,
                _fullsize: bool = False,
            ) -> tuple[
                NDArray[Any] | None,
                tuple[int, int, int, int, int],
                tuple[int, int, int, int],
            ]:
                # return decoded segment, its shape, and indices in image
                segmentindex, shape = indices(index)
                if data is None:
                    if _fullsize:
                        shape = pad_none(shape)
                    return data, segmentindex, shape
                data_array: NDArray[Any] = imagecodecs.jpeg_decode(
                    data,
                    bitspersample=self.bitspersample,
                    tables=jpegtables,
                    header=jpegheader,
                    colorspace=colorspace,
                    outcolorspace=outcolorspace,
                    shape=shape[1:3],
                )
                data_array = reshape(data_array, segmentindex, shape)
                if _fullsize:
                    data_array, shape = pad(data_array, shape)
                return data_array, segmentindex, shape

            return cache(decode_jpeg)

        if self.compression in {65000, 65001, 65002}:
            # EER decoder requires shape and extra args

            if self.compression == 65002:
                rlebits = int(self.tags.valueof(65007, 7))
                horzbits = int(self.tags.valueof(65008, 2))
                vertbits = int(self.tags.valueof(65009, 2))
            elif self.compression == 65001:
                rlebits = 7
                horzbits = 2
                vertbits = 2
            else:
                rlebits = 8
                horzbits = 2
                vertbits = 2

            def decode_eer(
                data: bytes | None,
                index: int,
                /,
                *,
                jpegtables: bytes | None = None,
                jpegheader: bytes | None = None,
                _fullsize: bool = False,
            ) -> tuple[
                NDArray[Any] | None,
                tuple[int, int, int, int, int],
                tuple[int, int, int, int],
            ]:
                # return decoded eer segment, its shape, and indices in image
                segmentindex, shape = indices(index)
                if data is None:
                    if _fullsize:
                        shape = pad_none(shape)
                    return data, segmentindex, shape
                data_array = decompress(
                    data,
                    shape=shape[1:3],
                    rlebits=rlebits,
                    horzbits=horzbits,
                    vertbits=vertbits,
                    superres=False,
                )  # type: ignore[call-arg, misc]
                return data_array.reshape(shape), segmentindex, shape

            return cache(decode_eer)

        if self.compression == 48124:
            # Jetraw requires pre-allocated output buffer

            def decode_jetraw(
                data: bytes | None,
                index: int,
                /,
                *,
                jpegtables: bytes | None = None,
                jpegheader: bytes | None = None,
                _fullsize: bool = False,
            ) -> tuple[
                NDArray[Any] | None,
                tuple[int, int, int, int, int],
                tuple[int, int, int, int],
            ]:
                # return decoded segment, its shape, and indices in image
                segmentindex, shape = indices(index)
                if data is None:
                    if _fullsize:
                        shape = pad_none(shape)
                    return data, segmentindex, shape
                data_array = numpy.zeros(shape, numpy.uint16)
                decompress(data, out=data_array)  # type: ignore[misc]
                return data_array.reshape(shape), segmentindex, shape

            return cache(decode_jetraw)

        if self.compression in TIFF.IMAGE_COMPRESSIONS:
            # presume codecs always return correct dtype, native byte order...
            if self.fillorder == 2:
                logger().debug(
                    f'{self!r} '
                    f'disabling LSB2MSB for compression {self.compression}'
                )
            if unpredict:
                logger().debug(
                    f'{self!r} '
                    f'disabling predictor for compression {self.compression}'
                )

            def decode_image(
                data: bytes | None,
                index: int,
                /,
                *,
                jpegtables: bytes | None = None,
                jpegheader: bytes | None = None,
                _fullsize: bool = False,
            ) -> tuple[
                NDArray[Any] | None,
                tuple[int, int, int, int, int],
                tuple[int, int, int, int],
            ]:
                # return decoded segment, its shape, and indices in image
                segmentindex, shape = indices(index)
                if data is None:
                    if _fullsize:
                        shape = pad_none(shape)
                    return data, segmentindex, shape
                data_array: NDArray[Any]
                data_array = decompress(data)  # type: ignore[misc]
                # del data
                data_array = reshape(data_array, segmentindex, shape)
                if _fullsize:
                    data_array, shape = pad(data_array, shape)
                return data_array, segmentindex, shape

            return cache(decode_image)

        dtype = numpy.dtype(self.parent.byteorder + self._dtype.char)

        if self.sampleformat == 5:
            # complex integer
            if unpredict is not None:
                raise NotImplementedError(
                    'unpredicting complex integers not supported'
                )

            itype = numpy.dtype(
                f'{self.parent.byteorder}i{self.bitspersample // 16}'
            )
            ftype = numpy.dtype(
                f'{self.parent.byteorder}f{dtype.itemsize // 2}'
            )

            def unpack(data: bytes, /) -> NDArray[Any]:
                # return complex integer as numpy.complex
                return numpy.frombuffer(data, itype).astype(ftype).view(dtype)

        elif self.bitspersample in {8, 16, 32, 64, 128}:
            # regular data types

            if (self.bitspersample * stwidth * samples) % 8:
                raise ValueError('data and sample size mismatch')
            if self.predictor in {3, 34894, 34895}:  # PREDICTOR.FLOATINGPOINT
                # floating-point horizontal differencing decoder needs
                # raw byte order
                dtype = numpy.dtype(self._dtype.char)

            def unpack(data: bytes, /) -> NDArray[Any]:
                # return numpy array from buffer
                try:
                    # read only numpy array
                    return numpy.frombuffer(data, dtype)
                except ValueError:
                    # for example, LZW strips may be missing EOI
                    bps = self.bitspersample // 8
                    size = (len(data) // bps) * bps
                    return numpy.frombuffer(data[:size], dtype)

        elif isinstance(self.bitspersample, tuple):
            # for example, RGB 565
            def unpack(data: bytes, /) -> NDArray[Any]:
                # return numpy array from packed integers
                return unpack_rgb(data, dtype, self.bitspersample)

        elif self.bitspersample == 24 and dtype.char == 'f':
            # float24
            if unpredict is not None:
                # floatpred_decode requires numpy.float24, which does not exist
                raise NotImplementedError('unpredicting float24 not supported')

            def unpack(data: bytes, /) -> NDArray[Any]:
                # return numpy.float32 array from float24
                return imagecodecs.float24_decode(
                    data, byteorder=self.parent.byteorder
                )

        else:
            # bilevel and packed integers
            def unpack(data: bytes, /) -> NDArray[Any]:
                # return NumPy array from packed integers
                return imagecodecs.packints_decode(
                    data, dtype, self.bitspersample, runlen=stwidth * samples
                )

        def decode_other(
            data: bytes | None,
            index: int,
            /,
            *,
            jpegtables: bytes | None = None,
            jpegheader: bytes | None = None,
            _fullsize: bool = False,
        ) -> tuple[
            NDArray[Any] | None,
            tuple[int, int, int, int, int],
            tuple[int, int, int, int],
        ]:
            # return decoded segment, its shape, and indices in image
            segmentindex, shape = indices(index)
            if data is None:
                if _fullsize:
                    shape = pad_none(shape)
                return data, segmentindex, shape
            if self.fillorder == 2:
                data = imagecodecs.bitorder_decode(data)
            if decompress is not None:
                # TODO: calculate correct size for packed integers
                size = shape[0] * shape[1] * shape[2] * shape[3]
                data = decompress(data, out=size * dtype.itemsize)
            data_array = unpack(data)  # type: ignore[arg-type]
            # del data
            data_array = reshape(data_array, segmentindex, shape)
            data_array = data_array.astype('=' + dtype.char, copy=False)
            if unpredict is not None:
                # unpredict is faster with native byte order
                data_array = unpredict(data_array, axis=-2, out=data_array)
            if _fullsize:
                data_array, shape = pad(data_array, shape)
            return data_array, segmentindex, shape

        return cache(decode_other)

    def segments(
        self,
        *,
        lock: threading.RLock | NullContext | None = None,
        maxworkers: int | None = None,
        func: Callable[..., Any] | None = None,  # TODO: type this
        sort: bool = False,
        buffersize: int | None = None,
        _fullsize: bool | None = None,
    ) -> Iterator[
        tuple[
            NDArray[Any] | None,
            tuple[int, int, int, int, int],
            tuple[int, int, int, int],
        ]
    ]:
        """Return iterator over decoded tiles or strips.

        Parameters:
            lock:
                Reentrant lock to synchronize file seeks and reads.
            maxworkers:
                Maximum number of threads to concurrently decode segments.
            func:
                Function to process decoded segment.
            sort:
                Read segments from file in order of their offsets.
            buffersize:
                Approximate number of bytes to read from file in one pass.
                The default is :py:attr:`_TIFF.BUFFERSIZE`.
            _fullsize:
                Internal use.

        Yields:
            - Decoded segment or None for empty segments.
            - Position of segment in image array of normalized shape
              (separate sample, depth, length, width, contig sample).
            - Shape of segment (depth, length, width, contig samples).
              The shape of strips depends on their linear index.

        """
        keyframe = self.keyframe  # self or keyframe
        fh = self.parent.filehandle
        if lock is None:
            lock = fh.lock
        if _fullsize is None:
            _fullsize = keyframe.is_tiled

        decodeargs: dict[str, Any] = {'_fullsize': bool(_fullsize)}
        if keyframe.compression in {6, 7, 34892, 33007}:  # JPEG
            decodeargs['jpegtables'] = self.jpegtables
            decodeargs['jpegheader'] = keyframe.jpegheader

        if func is None:

            def decode(args, decodeargs=decodeargs, decode=keyframe.decode):
                return decode(*args, **decodeargs)

        else:

            def decode(args, decodeargs=decodeargs, decode=keyframe.decode):
                return func(decode(*args, **decodeargs))

        if maxworkers is None or maxworkers < 1:
            maxworkers = keyframe.maxworkers
        if maxworkers < 2:
            for segment in fh.read_segments(
                self.dataoffsets,
                self.databytecounts,
                lock=lock,
                sort=sort,
                buffersize=buffersize,
                flat=True,
            ):
                yield decode(segment)
        else:
            # reduce memory overhead by processing chunks of up to
            # buffersize of segments because ThreadPoolExecutor.map is not
            # collecting iterables lazily
            with ThreadPoolExecutor(maxworkers) as executor:
                for segments in fh.read_segments(
                    self.dataoffsets,
                    self.databytecounts,
                    lock=lock,
                    sort=sort,
                    buffersize=buffersize,
                    flat=False,
                ):
                    yield from executor.map(decode, segments)

    def asarray(
        self,
        *,
        out: OutputType = None,
        squeeze: bool = True,
        lock: threading.RLock | NullContext | None = None,
        maxworkers: int | None = None,
        buffersize: int | None = None,
    ) -> NDArray[Any]:
        """Return image from page as NumPy array.

        Parameters:
            out:
                Specifies how image array is returned.
                By default, a new NumPy array is created.
                If a *numpy.ndarray*, a writable array to which the image
                is copied.
                If *'memmap'*, directly memory-map the image data in the
                file if possible; else create a memory-mapped array in a
                temporary file.
                If a *string* or *open file*, the file used to create a
                memory-mapped array.
            squeeze:
                Remove all length-1 dimensions (except X and Y) from
                image array.
                If *False*, return the image array with normalized
                5-dimensional shape :py:attr:`TiffPage.shaped`.
            lock:
                Reentrant lock to synchronize seeks and reads from file.
                The default is the lock of the parent's file handle.
            maxworkers:
                Maximum number of threads to concurrently decode segments.
                If *None* or *0*, use up to :py:attr:`_TIFF.MAXWORKERS`
                threads. See remarks in :py:meth:`TiffFile.asarray`.
            buffersize:
                Approximate number of bytes to read from file in one pass.
                The default is :py:attr:`_TIFF.BUFFERSIZE`.

        Returns:
            NumPy array of decompressed, unpredicted, and unpacked image data
            read from Strip/Tile Offsets/ByteCounts, formatted according to
            shape and dtype metadata found in tags and arguments.
            Photometric conversion, premultiplied alpha, orientation, and
            colorimetry corrections are not applied.
            Specifically, CMYK images are not converted to RGB, MinIsWhite
            images are not inverted, color palettes are not applied,
            gamma is not corrected, and CFA images are not demosaciced.
            Exception are YCbCr JPEG compressed images, which are converted to
            RGB.

        Raises:
            ValueError:
                Format of image in file is not supported and cannot be decoded.

        """
        keyframe = self.keyframe  # self or keyframe

        if 0 in keyframe.shaped or keyframe._dtype is None:
            return numpy.empty((0,), keyframe.dtype)

        if len(self.dataoffsets) == 0:
            raise TiffFileError('missing data offset')

        fh = self.parent.filehandle
        if lock is None:
            lock = fh.lock

        if (
            isinstance(out, str)
            and out == 'memmap'
            and keyframe.is_memmappable
        ):
            # direct memory map array in file
            with lock:
                closed = fh.closed
                if closed:
                    warnings.warn(
                        f'{self!r} reading array from closed file', UserWarning
                    )
                    fh.open()
                result = fh.memmap_array(
                    keyframe.parent.byteorder + keyframe._dtype.char,
                    keyframe.shaped,
                    offset=self.dataoffsets[0],
                )

        elif keyframe.is_contiguous:
            # read contiguous bytes to array
            if keyframe.is_subsampled:
                raise NotImplementedError('chroma subsampling not supported')
            if out is not None:
                out = create_output(out, keyframe.shaped, keyframe._dtype)
            with lock:
                closed = fh.closed
                if closed:
                    warnings.warn(
                        f'{self!r} reading array from closed file', UserWarning
                    )
                    fh.open()
                fh.seek(self.dataoffsets[0])
                result = fh.read_array(
                    keyframe.parent.byteorder + keyframe._dtype.char,
                    product(keyframe.shaped),
                    out=out,
                )
            if keyframe.fillorder == 2:
                result = imagecodecs.bitorder_decode(result, out=result)
            if keyframe.predictor != 1:
                # predictors without compression
                unpredict = TIFF.UNPREDICTORS[keyframe.predictor]
                if keyframe.predictor == 1:
                    result = unpredict(result, axis=-2, out=result)
                else:
                    # floatpred cannot decode in-place
                    out = unpredict(result, axis=-2, out=result)
                    result[:] = out

        elif (
            keyframe.jpegheader is not None
            and keyframe is self
            and 273 in self.tags  # striped ...
            and self.is_tiled  # but reported as tiled
            # TODO: imagecodecs can decode larger JPEG
            and self.imagewidth <= 65500
            and self.imagelength <= 65500
        ):
            # decode the whole NDPI JPEG strip
            with lock:
                closed = fh.closed
                if closed:
                    warnings.warn(
                        f'{self!r} reading array from closed file', UserWarning
                    )
                    fh.open()
                fh.seek(self.tags[273].value[0])  # StripOffsets
                data = fh.read(self.tags[279].value[0])  # StripByteCounts
            decompress = TIFF.DECOMPRESSORS[self.compression]
            result = decompress(
                data,
                bitspersample=self.bitspersample,
                out=out,
                # shape=(self.imagelength, self.imagewidth)
            )
            del data

        else:
            # decode individual strips or tiles
            with lock:
                closed = fh.closed
                if closed:
                    warnings.warn(
                        f'{self!r} reading array from closed file', UserWarning
                    )
                    fh.open()
                keyframe.decode  # init TiffPage.decode function under lock

            result = create_output(out, keyframe.shaped, keyframe._dtype)

            def func(
                decoderesult: tuple[
                    NDArray[Any] | None,
                    tuple[int, int, int, int, int],
                    tuple[int, int, int, int],
                ],
                keyframe: TiffPage = keyframe,
                out: NDArray[Any] = result,
            ) -> None:
                # copy decoded segments to output array
                segment, (s, d, h, w, _), shape = decoderesult
                if segment is None:
                    out[
                        s, d : d + shape[0], h : h + shape[1], w : w + shape[2]
                    ] = keyframe.nodata
                else:
                    out[
                        s, d : d + shape[0], h : h + shape[1], w : w + shape[2]
                    ] = segment[
                        : keyframe.imagedepth - d,
                        : keyframe.imagelength - h,
                        : keyframe.imagewidth - w,
                    ]
                # except IndexError:
                #     pass  # corrupted file, for example, with too many strips

            for _ in self.segments(
                func=func,
                lock=lock,
                maxworkers=maxworkers,
                buffersize=buffersize,
                sort=True,
                _fullsize=False,
            ):
                pass

        result.shape = keyframe.shaped
        if squeeze:
            try:
                result.shape = keyframe.shape
            except ValueError as exc:
                logger().warning(
                    f'{self!r} <asarray> failed to reshape '
                    f'{result.shape} to {keyframe.shape}, raised {exc!r:.128}'
                )

        if closed:
            # TODO: close file if an exception occurred above
            fh.close()
        return result

    def aszarr(self, **kwargs: Any) -> ZarrTiffStore:
        """Return image from page as Zarr 2 store.

        Parameters:
            **kwarg: Passed to :py:class:`ZarrTiffStore`.

        """
        return ZarrTiffStore(self, **kwargs)

    def asrgb(
        self,
        *,
        uint8: bool = False,
        alpha: Container[int] | None = None,
        **kwargs: Any,
    ) -> NDArray[Any]:
        """Return image as RGB(A). Work in progress. Do not use.

        :meta private:

        """
        data = self.asarray(**kwargs)
        keyframe = self.keyframe  # self or keyframe

        if keyframe.photometric == PHOTOMETRIC.PALETTE:
            colormap = keyframe.colormap
            if colormap is None:
                raise ValueError('no colormap')
            if (
                colormap.shape[1] < 2**keyframe.bitspersample
                or keyframe.dtype is None
                or keyframe.dtype.char not in 'BH'
            ):
                raise ValueError('cannot apply colormap')
            if uint8:
                if colormap.max() > 255:
                    colormap >>= 8
                colormap = colormap.astype(numpy.uint8)
            if 'S' in keyframe.axes:
                data = data[..., 0] if keyframe.planarconfig == 1 else data[0]
            data = apply_colormap(data, colormap)

        elif keyframe.photometric == PHOTOMETRIC.RGB:
            if keyframe.extrasamples:
                if alpha is None:
                    alpha = EXTRASAMPLE
                for i, exs in enumerate(keyframe.extrasamples):
                    if exs in EXTRASAMPLE:
                        if keyframe.planarconfig == 1:
                            data = data[..., [0, 1, 2, 3 + i]]
                        else:
                            data = data[:, [0, 1, 2, 3 + i]]
                        break
            else:
                if keyframe.planarconfig == 1:
                    data = data[..., :3]
                else:
                    data = data[:, :3]
            # TODO: convert to uint8?

        elif keyframe.photometric == PHOTOMETRIC.MINISBLACK:
            raise NotImplementedError
        elif keyframe.photometric == PHOTOMETRIC.MINISWHITE:
            raise NotImplementedError
        elif keyframe.photometric == PHOTOMETRIC.SEPARATED:
            raise NotImplementedError
        else:
            raise NotImplementedError
        return data

    def _gettags(
        self,
        codes: Container[int] | None = None,
        /,
        lock: threading.RLock | None = None,
    ) -> list[tuple[int, TiffTag]]:
        """Return list of (code, TiffTag)."""
        return [
            (tag.code, tag)
            for tag in self.tags
            if codes is None or tag.code in codes
        ]

    def _nextifd(self) -> int:
        """Return offset to next IFD from file."""
        fh = self.parent.filehandle
        tiff = self.parent.tiff
        fh.seek(self.offset)
        tagno = struct.unpack(tiff.tagnoformat, fh.read(tiff.tagnosize))[0]
        fh.seek(self.offset + tiff.tagnosize + tagno * tiff.tagsize)
        return int(
            struct.unpack(tiff.offsetformat, fh.read(tiff.offsetsize))[0]
        )

    def aspage(self) -> TiffPage:
        """Return TiffPage instance."""
        return self

    @property
    def index(self) -> int:
        """Index of page in IFD chain."""
        return self._index[-1]

    @property
    def treeindex(self) -> tuple[int, ...]:
        """Index of page in IFD tree."""
        return self._index

    @property
    def keyframe(self) -> TiffPage:
        """Self."""
        return self

    @keyframe.setter
    def keyframe(self, index: TiffPage) -> None:
        return

    @property
    def name(self) -> str:
        """Name of image array."""
        index = self._index if len(self._index) > 1 else self._index[0]
        return f'TiffPage {index}'

    @property
    def ndim(self) -> int:
        """Number of dimensions in image array."""
        return len(self.shape)

    @cached_property
    def dims(self) -> tuple[str, ...]:
        """Names of dimensions in image array."""
        names = TIFF.AXES_NAMES
        return tuple(names[ax] for ax in self.axes)

    @cached_property
    def sizes(self) -> dict[str, int]:
        """Ordered map of dimension names to lengths."""
        shape = self.shape
        names = TIFF.AXES_NAMES
        return {names[ax]: shape[i] for i, ax in enumerate(self.axes)}

    @cached_property
    def coords(self) -> dict[str, NDArray[Any]]:
        """Ordered map of dimension names to coordinate arrays."""
        resolution = self.get_resolution()
        coords: dict[str, NDArray[Any]] = {}

        for ax, size in zip(self.axes, self.shape):
            name = TIFF.AXES_NAMES[ax]
            value = None
            step: int | float = 1

            if ax == 'X':
                step = resolution[0]
            elif ax == 'Y':
                step = resolution[1]
            elif ax == 'S':
                value = self._sample_names()
            elif ax == 'Z':
                # a ZResolution tag doesn't exist.
                # use XResolution if it agrees with YResolution
                if resolution[0] == resolution[1]:
                    step = resolution[0]

            if value is not None:
                coords[name] = numpy.asarray(value)
            elif step == 0 or step == 1 or size == 0:
                coords[name] = numpy.arange(size)
            else:
                coords[name] = numpy.linspace(
                    0, size / step, size, endpoint=False, dtype=numpy.float32
                )
            assert len(coords[name]) == size
        return coords

    @cached_property
    def attr(self) -> dict[str, Any]:
        """Arbitrary metadata associated with image array."""
        # TODO: what to return?
        return {}

    @cached_property
    def size(self) -> int:
        """Number of elements in image array."""
        return product(self.shape)

    @cached_property
    def nbytes(self) -> int:
        """Number of bytes in image array."""
        if self.dtype is None:
            return 0
        return self.size * self.dtype.itemsize

    @property
    def colormap(self) -> NDArray[numpy.uint16] | None:
        """Value of Colormap tag."""
        return self.tags.valueof(320)

    @property
    def iccprofile(self) -> bytes | None:
        """Value of InterColorProfile tag."""
        return self.tags.valueof(34675)

    @property
    def transferfunction(self) -> NDArray[numpy.uint16] | None:
        """Value of TransferFunction tag."""
        return self.tags.valueof(301)

    def get_resolution(
        self,
        unit: RESUNIT | int | str | None = None,
        scale: float | int | None = None,
    ) -> tuple[int | float, int | float]:
        """Return number of pixels per unit in X and Y dimensions.

        By default, the XResolution and YResolution tag values are returned.
        Missing tag values are set to 1.

        Parameters:
            unit:
                Unit of measurement of returned values.
                The default is the value of the ResolutionUnit tag.
            scale:
                Factor to convert resolution values to meter unit.
                The default is determined from the ResolutionUnit tag.

        """
        scales = {
            1: 1,  # meter, no unit
            2: 100 / 2.54,  # INCH
            3: 100,  # CENTIMETER
            4: 1000,  # MILLIMETER
            5: 1000000,  # MICROMETER
        }
        if unit is not None:
            unit = enumarg(RESUNIT, unit)
            try:
                if scale is None:
                    resolutionunit = self.tags.valueof(296, default=2)
                    scale = scales[resolutionunit]
            except Exception as exc:
                logger().warning(
                    f'{self!r} <get_resolution> raised {exc!r:.128}'
                )
                scale = 1
            else:
                scale2 = scales[unit]
                if scale % scale2 == 0:
                    scale //= scale2
                else:
                    scale /= scale2
        elif scale is None:
            scale = 1

        resolution: list[int | float] = []
        n: int
        d: int
        for code in 282, 283:
            try:
                n, d = self.tags.valueof(code, default=(1, 1))
                if d == 0:
                    value = n * scale
                elif n % d == 0:
                    value = n // d * scale
                else:
                    value = n / d * scale
            except Exception:
                value = 1
            resolution.append(value)
        return resolution[0], resolution[1]

    @cached_property
    def resolution(self) -> tuple[float, float]:
        """Number of pixels per resolutionunit in X and Y directions."""
        # values are returned in (somewhat unexpected) XY order to
        # keep symmetry with the TiffWriter.write resolution argument
        resolution = self.get_resolution()
        return float(resolution[0]), float(resolution[1])

    @property
    def resolutionunit(self) -> int:
        """Unit of measurement for X and Y resolutions."""
        return self.tags.valueof(296, default=2)

    @property
    def datetime(self) -> DateTime | None:
        """Date and time of image creation."""
        value = self.tags.valueof(306)
        if value is None:
            return None
        try:
            return strptime(value)
        except Exception:
            pass
        return None

    @property
    def tile(self) -> tuple[int, ...] | None:
        """Tile depth, length, and width."""
        if not self.is_tiled:
            return None
        if self.tiledepth > 1:
            return (self.tiledepth, self.tilelength, self.tilewidth)
        return (self.tilelength, self.tilewidth)

    @cached_property
    def chunks(self) -> tuple[int, ...]:
        """Shape of images in tiles or strips."""
        shape: list[int] = []
        if self.tiledepth > 1:
            shape.append(self.tiledepth)
        if self.is_tiled:
            shape.extend((self.tilelength, self.tilewidth))
        else:
            shape.extend((self.rowsperstrip, self.imagewidth))
        if self.planarconfig == 1 and self.samplesperpixel > 1:
            shape.append(self.samplesperpixel)
        return tuple(shape)

    @cached_property
    def chunked(self) -> tuple[int, ...]:
        """Shape of chunked image."""
        shape: list[int] = []
        if self.planarconfig == 2 and self.samplesperpixel > 1:
            shape.append(self.samplesperpixel)
        if self.is_tiled:
            if self.imagedepth > 1:
                shape.append(
                    (self.imagedepth + self.tiledepth - 1) // self.tiledepth
                )
            shape.append(
                (self.imagelength + self.tilelength - 1) // self.tilelength
            )
            shape.append(
                (self.imagewidth + self.tilewidth - 1) // self.tilewidth
            )
        else:
            if self.imagedepth > 1:
                shape.append(self.imagedepth)
            shape.append(
                (self.imagelength + self.rowsperstrip - 1) // self.rowsperstrip
            )
            shape.append(1)
        if self.planarconfig == 1 and self.samplesperpixel > 1:
            shape.append(1)
        return tuple(shape)

    @cached_property
    def hash(self) -> int:
        """Checksum to identify pages in same series.

        Pages with the same hash can use the same decode function.
        The hash is calculated from the following properties:
        :py:attr:`TiffFile.tiff`,
        :py:attr:`TiffPage.shaped`,
        :py:attr:`TiffPage.rowsperstrip`,
        :py:attr:`TiffPage.tilewidth`,
        :py:attr:`TiffPage.tilelength`,
        :py:attr:`TiffPage.tiledepth`,
        :py:attr:`TiffPage.sampleformat`,
        :py:attr:`TiffPage.bitspersample`,
        :py:attr:`TiffPage.fillorder`,
        :py:attr:`TiffPage.predictor`,
        :py:attr:`TiffPage.compression`,
        :py:attr:`TiffPage.extrasamples`, and
        :py:attr:`TiffPage.photometric`.

        """
        return hash(
            self.shaped
            + (
                self.parent.tiff,
                self.rowsperstrip,
                self.tilewidth,
                self.tilelength,
                self.tiledepth,
                self.sampleformat,
                self.bitspersample,
                self.fillorder,
                self.predictor,
                self.compression,
                self.extrasamples,
                self.photometric,
            )
        )

    @cached_property
    def pages(self) -> TiffPages | None:
        """Sequence of sub-pages, SubIFDs."""
        if 330 not in self.tags:
            return None
        return TiffPages(self, index=self.index)

    @cached_property
    def maxworkers(self) -> int:
        """Maximum number of threads for decoding segments.

        A value of 0 disables multi-threading also when stacking pages.

        """
        if self.is_contiguous or self.dtype is None:
            return 0
        if self.compression in TIFF.IMAGE_COMPRESSIONS:
            return min(TIFF.MAXWORKERS, len(self.dataoffsets))
        bytecount = product(self.chunks) * self.dtype.itemsize
        if bytecount < 2048:
            # disable multi-threading for small segments
            return 0
        if self.compression == 5 and bytecount < 14336:
            # disable multi-threading for small LZW compressed segments
            return 0
        if len(self.dataoffsets) < 4:
            return 1
        if self.compression != 1 or self.fillorder != 1 or self.predictor != 1:
            if imagecodecs is not None:
                return min(TIFF.MAXWORKERS, len(self.dataoffsets))
        return 2  # optimum for large number of uncompressed tiles

    @cached_property
    def is_contiguous(self) -> bool:
        """Image data is stored contiguously.

        Contiguous image data can be read from
        ``offset=TiffPage.dataoffsets[0]`` with ``size=TiffPage.nbytes``.
        Excludes prediction and fillorder.

        """
        if (
            self.sampleformat == 5
            or self.compression != 1
            or self.bitspersample not in {8, 16, 32, 64}
        ):
            return False
        if 322 in self.tags:  # TileWidth
            if (
                self.imagewidth != self.tilewidth
                or self.imagelength % self.tilelength
                or self.tilewidth % 16
                or self.tilelength % 16
            ):
                return False
            if (
                32997 in self.tags  # ImageDepth
                and 32998 in self.tags  # TileDepth
                and (
                    self.imagelength != self.tilelength
                    or self.imagedepth % self.tiledepth
                )
            ):
                return False
        offsets = self.dataoffsets
        bytecounts = self.databytecounts
        if len(offsets) == 0:
            return False
        if len(offsets) == 1:
            return True
        if self.is_stk or self.is_lsm:
            return True
        if sum(bytecounts) != self.nbytes:
            return False
        if all(
            bytecounts[i] != 0 and offsets[i] + bytecounts[i] == offsets[i + 1]
            for i in range(len(offsets) - 1)
        ):
            return True
        return False

    @cached_property
    def is_final(self) -> bool:
        """Image data are stored in final form. Excludes byte-swapping."""
        return (
            self.is_contiguous
            and self.fillorder == 1
            and self.predictor == 1
            and not self.is_subsampled
        )

    @cached_property
    def is_memmappable(self) -> bool:
        """Image data in file can be memory-mapped to NumPy array."""
        return (
            self.parent.filehandle.is_file
            and self.is_final
            # and (self.bitspersample == 8 or self.parent.isnative)
            # aligned?
            and self.dtype is not None
            and self.dataoffsets[0] % self.dtype.itemsize == 0
        )

    def __repr__(self) -> str:
        index = self._index if len(self._index) > 1 else self._index[0]
        return f'<tifffile.TiffPage {index} @{self.offset}>'

    def __str__(self) -> str:
        return self._str()

    def _str(self, detail: int = 0, width: int = 79) -> str:
        """Return string containing information about TiffPage."""
        if self.keyframe != self:
            return TiffFrame._str(
                self, detail, width  # type: ignore[arg-type]
            )
        attr = ''
        for name in ('memmappable', 'final', 'contiguous'):
            attr = getattr(self, 'is_' + name)
            if attr:
                attr = name.upper()
                break

        def tostr(name: str, /, skip: int = 1) -> str:
            obj = getattr(self, name)
            if obj == skip:
                return ''
            try:
                value = getattr(obj, 'name')
            except AttributeError:
                return ''
            return str(value)

        info = '  '.join(
            s.lower()
            for s in (
                'x'.join(str(i) for i in self.shape),
                f'{SAMPLEFORMAT(self.sampleformat).name}{self.bitspersample}',
                ' '.join(
                    i
                    for i in (
                        PHOTOMETRIC(self.photometric).name,
                        'REDUCED' if self.is_reduced else '',
                        'MASK' if self.is_mask else '',
                        'TILED' if self.is_tiled else '',
                        tostr('compression'),
                        tostr('planarconfig'),
                        tostr('predictor'),
                        tostr('fillorder'),
                    )
                    + (attr,)
                    if i
                ),
                '|'.join(f.upper() for f in sorted(self.flags)),
            )
            if s
        )
        index = self._index if len(self._index) > 1 else self._index[0]
        info = f'TiffPage {index} @{self.offset}  {info}'
        if detail <= 0:
            return info
        info_list = [info, self.tags._str(detail + 1, width=width)]
        if detail > 1:
            for name in ('ndpi',):
                name = name + '_tags'
                attr = getattr(self, name, '')
                if attr:
                    info_list.append(
                        f'{name.upper()}\n'
                        f'{pformat(attr, width=width, height=detail * 8)}'
                    )
        if detail > 3:
            try:
                data = self.asarray()
                info_list.append(
                    f'DATA\n'
                    f'{pformat(data, width=width, height=detail * 8)}'
                )
            except Exception:
                pass
        return '\n\n'.join(info_list)

    def _sample_names(self) -> list[str] | None:
        """Return names of samples."""
        if 'S' not in self.axes:
            return None
        samples = self.shape[self.axes.find('S')]
        extrasamples = len(self.extrasamples)
        if samples < 1 or extrasamples > 2:
            return None
        if self.photometric == 0:
            names = ['WhiteIsZero']
        elif self.photometric == 1:
            names = ['BlackIsZero']
        elif self.photometric == 2:
            names = ['Red', 'Green', 'Blue']
        elif self.photometric == 5:
            names = ['Cyan', 'Magenta', 'Yellow', 'Black']
        elif self.photometric == 6:
            if self.compression in {6, 7, 34892, 33007}:
                # YCBCR -> RGB for JPEG
                names = ['Red', 'Green', 'Blue']
            else:
                names = ['Luma', 'Cb', 'Cr']
        else:
            return None
        if extrasamples > 0:
            names += [enumarg(EXTRASAMPLE, self.extrasamples[0]).name.title()]
        if extrasamples > 1:
            names += [enumarg(EXTRASAMPLE, self.extrasamples[1]).name.title()]
        if len(names) != samples:
            return None
        return names

    @cached_property
    def flags(self) -> set[str]:
        r"""Set of ``is\_\*`` properties that are True."""
        return {
            name.lower()
            for name in TIFF.PAGE_FLAGS
            if getattr(self, 'is_' + name)
        }

    @cached_property
    def andor_tags(self) -> dict[str, Any] | None:
        """Consolidated metadata from Andor tags."""
        if not self.is_andor:
            return None
        result = {'Id': self.tags[4864].value}  # AndorId
        for tag in self.tags:  # list(self.tags.values()):
            code = tag.code
            if not 4864 < code < 5031:
                continue
            name = tag.name
            name = name[5:] if len(name) > 5 else name
            result[name] = tag.value
            # del self.tags[code]
        return result

    @cached_property
    def epics_tags(self) -> dict[str, Any] | None:
        """Consolidated metadata from EPICS areaDetector tags.

        Use the :py:func:`epics_datetime` function to get a datetime object
        from the epicsTSSec and epicsTSNsec tags.

        """
        if not self.is_epics:
            return None
        result = {}
        for tag in self.tags:  # list(self.tags.values()):
            code = tag.code
            if not 65000 <= code < 65500:
                continue
            value = tag.value
            if code == 65000:
                # not a POSIX timestamp
                # https://github.com/bluesky/area-detector-handlers/issues/20
                result['timeStamp'] = float(value)
            elif code == 65001:
                result['uniqueID'] = int(value)
            elif code == 65002:
                result['epicsTSSec'] = int(value)
            elif code == 65003:
                result['epicsTSNsec'] = int(value)
            else:
                key, value = value.split(':', 1)
                result[key] = astype(value)
            # del self.tags[code]
        return result

    @cached_property
    def ndpi_tags(self) -> dict[str, Any] | None:
        """Consolidated metadata from Hamamatsu NDPI tags."""
        # TODO: parse 65449 ini style comments
        if not self.is_ndpi:
            return None
        tags = self.tags
        result = {}
        for name in ('Make', 'Model', 'Software'):
            result[name] = tags[name].value
        for code, name in TIFF.NDPI_TAGS.items():
            if code in tags:
                result[name] = tags[code].value
                # del tags[code]
        if 'McuStarts' in result:
            mcustarts = result['McuStarts']
            if 'McuStartsHighBytes' in result:
                high = result['McuStartsHighBytes'].astype(numpy.uint64)
                high <<= 32
                mcustarts = mcustarts.astype(numpy.uint64)
                mcustarts += high
                del result['McuStartsHighBytes']
            result['McuStarts'] = mcustarts
        return result

    @cached_property
    def geotiff_tags(self) -> dict[str, Any] | None:
        """Consolidated metadata from GeoTIFF tags."""
        if not self.is_geotiff:
            return None
        tags = self.tags

        gkd = tags.valueof(34735)  # GeoKeyDirectoryTag
        if gkd is None or len(gkd) < 2 or gkd[0] != 1:
            logger().warning(f'{self!r} invalid GeoKeyDirectoryTag')
            return {}

        result = {
            'KeyDirectoryVersion': gkd[0],
            'KeyRevision': gkd[1],
            'KeyRevisionMinor': gkd[2],
            # 'NumberOfKeys': gkd[3],
        }
        # deltags = ['GeoKeyDirectoryTag']
        geokeys = TIFF.GEO_KEYS
        geocodes = TIFF.GEO_CODES
        for index in range(gkd[3]):
            try:
                keyid, tagid, count, offset = gkd[
                    4 + index * 4 : index * 4 + 8
                ]
            except Exception as exc:
                logger().warning(
                    f'{self!r} corrupted GeoKeyDirectoryTag '
                    f'raised {exc!r:.128}'
                )
                continue
            if tagid == 0:
                value = offset
            else:
                try:
                    value = tags[tagid].value[offset : offset + count]
                except TiffFileError as exc:
                    logger().warning(
                        f'{self!r} corrupted GeoKeyDirectoryTag {tagid} '
                        f'raised {exc!r:.128}'
                    )
                    continue
                except KeyError as exc:
                    logger().warning(
                        f'{self!r} GeoKeyDirectoryTag {tagid} not found, '
                        f'raised {exc!r:.128}'
                    )
                    continue
                if tagid == 34737 and count > 1 and value[-1] == '|':
                    value = value[:-1]
                value = value if count > 1 else value[0]
            if keyid in geocodes:
                try:
                    value = geocodes[keyid](value)
                except Exception:
                    pass
            try:
                key = geokeys(keyid).name
            except ValueError:
                key = keyid
            result[key] = value

        value = tags.valueof(33920)  # IntergraphMatrixTag
        if value is not None:
            value = numpy.array(value)
            if value.size == 16:
                value = value.reshape((4, 4)).tolist()
            result['IntergraphMatrix'] = value

        value = tags.valueof(33550)  # ModelPixelScaleTag
        if value is not None:
            result['ModelPixelScale'] = numpy.array(value).tolist()

        value = tags.valueof(33922)  # ModelTiepointTag
        if value is not None:
            value = numpy.array(value).reshape((-1, 6)).squeeze().tolist()
            result['ModelTiepoint'] = value

        value = tags.valueof(34264)  # ModelTransformationTag
        if value is not None:
            value = numpy.array(value).reshape((4, 4)).tolist()
            result['ModelTransformation'] = value

        # if 33550 in tags and 33922 in tags:
        #     sx, sy, sz = tags[33550].value  # ModelPixelScaleTag
        #     tiepoints = tags[33922].value  # ModelTiepointTag
        #     transforms = []
        #     for tp in range(0, len(tiepoints), 6):
        #         i, j, k, x, y, z = tiepoints[tp : tp + 6]
        #         transforms.append(
        #             [
        #                 [sx, 0.0, 0.0, x - i * sx],
        #                 [0.0, -sy, 0.0, y + j * sy],
        #                 [0.0, 0.0, sz, z - k * sz],
        #                 [0.0, 0.0, 0.0, 1.0],
        #             ]
        #         )
        #     if len(tiepoints) == 6:
        #         transforms = transforms[0]
        #     result['ModelTransformation'] = transforms

        rpcc = tags.valueof(50844)  # RPCCoefficientTag
        if rpcc is not None:
            result['RPCCoefficient'] = {
                'ERR_BIAS': rpcc[0],
                'ERR_RAND': rpcc[1],
                'LINE_OFF': rpcc[2],
                'SAMP_OFF': rpcc[3],
                'LAT_OFF': rpcc[4],
                'LONG_OFF': rpcc[5],
                'HEIGHT_OFF': rpcc[6],
                'LINE_SCALE': rpcc[7],
                'SAMP_SCALE': rpcc[8],
                'LAT_SCALE': rpcc[9],
                'LONG_SCALE': rpcc[10],
                'HEIGHT_SCALE': rpcc[11],
                'LINE_NUM_COEFF': rpcc[12:33],
                'LINE_DEN_COEFF ': rpcc[33:53],
                'SAMP_NUM_COEFF': rpcc[53:73],
                'SAMP_DEN_COEFF': rpcc[73:],
            }
        return result

    @cached_property
    def shaped_description(self) -> str | None:
        """Description containing array shape if exists, else None."""
        for description in (self.description, self.description1):
            if not description or '"mibi.' in description:
                return None
            if description[:1] == '{' and '"shape":' in description:
                return description
            if description[:6] == 'shape=':
                return description
        return None

    @cached_property
    def imagej_description(self) -> str | None:
        """ImageJ description if exists, else None."""
        for description in (self.description, self.description1):
            if not description:
                return None
            if description[:7] == 'ImageJ=' or description[:7] == 'SCIFIO=':
                return description
        return None

    @cached_property
    def is_jfif(self) -> bool:
        """JPEG compressed segments contain JFIF metadata."""
        if (
            self.compression not in {6, 7, 34892, 33007}
            or len(self.dataoffsets) < 1
            or self.dataoffsets[0] == 0
            or len(self.databytecounts) < 1
            or self.databytecounts[0] < 11
        ):
            return False
        fh = self.parent.filehandle
        fh.seek(self.dataoffsets[0] + 6)
        data = fh.read(4)
        return data == b'JFIF'  # or data == b'Exif'

    @property
    def is_frame(self) -> bool:
        """Object is :py:class:`TiffFrame` instance."""
        return False

    @property
    def is_virtual(self) -> bool:
        """Page does not have IFD structure in file."""
        return False

    @property
    def is_subifd(self) -> bool:
        """Page is SubIFD of another page."""
        return len(self._index) > 1

    @property
    def is_reduced(self) -> bool:
        """Page is reduced image of another image."""
        return bool(self.subfiletype & 0b1)

    @property
    def is_multipage(self) -> bool:
        """Page is part of multi-page image."""
        return bool(self.subfiletype & 0b10)

    @property
    def is_mask(self) -> bool:
        """Page is transparency mask for another image."""
        return bool(self.subfiletype & 0b100)

    @property
    def is_mrc(self) -> bool:
        """Page is part of Mixed Raster Content."""
        return bool(self.subfiletype & 0b1000)

    @property
    def is_tiled(self) -> bool:
        """Page contains tiled image."""
        return self.tilewidth > 0  # return 322 in self.tags  # TileWidth

    @property
    def is_subsampled(self) -> bool:
        """Page contains chroma subsampled image."""
        if self.subsampling is not None:
            return self.subsampling != (1, 1)
        return self.photometric == 6  # YCbCr
        # RGB JPEG usually stored as subsampled YCbCr
        # self.compression == 7
        # and self.photometric == 2
        # and self.planarconfig == 1

    @property
    def is_imagej(self) -> bool:
        """Page contains ImageJ description metadata."""
        return self.imagej_description is not None

    @property
    def is_shaped(self) -> bool:
        """Page contains Tifffile JSON metadata."""
        return self.shaped_description is not None

    @property
    def is_mdgel(self) -> bool:
        """Page contains MDFileTag tag."""
        return (
            37701 not in self.tags  # AgilentBinary
            and 33445 in self.tags  # MDFileTag
        )

    @property
    def is_agilent(self) -> bool:
        """Page contains Agilent Technologies tags."""
        # tag 270 and 285 contain color names
        return 285 in self.tags and 37701 in self.tags  # AgilentBinary

    @property
    def is_mediacy(self) -> bool:
        """Page contains Media Cybernetics Id tag."""
        tag = self.tags.get(50288)  # MC_Id
        try:
            return tag is not None and tag.value[:7] == b'MC TIFF'
        except Exception:
            return False

    @property
    def is_stk(self) -> bool:
        """Page contains UIC1Tag tag."""
        return 33628 in self.tags

    @property
    def is_lsm(self) -> bool:
        """Page contains CZ_LSMINFO tag."""
        return 34412 in self.tags

    @property
    def is_fluoview(self) -> bool:
        """Page contains FluoView MM_STAMP tag."""
        return 34362 in self.tags

    @property
    def is_nih(self) -> bool:
        """Page contains NIHImageHeader tag."""
        return 43314 in self.tags

    @property
    def is_volumetric(self) -> bool:
        """Page contains SGI ImageDepth tag with value > 1."""
        return self.imagedepth > 1

    @property
    def is_vista(self) -> bool:
        """Software tag is 'ISS Vista'."""
        return self.software == 'ISS Vista'

    @property
    def is_metaseries(self) -> bool:
        """Page contains MDS MetaSeries metadata in ImageDescription tag."""
        if self.index != 0 or self.software != 'MetaSeries':
            return False
        d = self.description
        return d.startswith('<MetaData>') and d.endswith('</MetaData>')

    @property
    def is_ome(self) -> bool:
        """Page contains OME-XML in ImageDescription tag."""
        if self.index != 0 or not self.description:
            return False
        return self.description[-10:].strip().endswith('OME>')

    @property
    def is_scn(self) -> bool:
        """Page contains Leica SCN XML in ImageDescription tag."""
        if self.index != 0 or not self.description:
            return False
        return self.description[-10:].strip().endswith('</scn>')

    @property
    def is_micromanager(self) -> bool:
        """Page contains MicroManagerMetadata tag."""
        return 51123 in self.tags

    @property
    def is_andor(self) -> bool:
        """Page contains Andor Technology tags 4864-5030."""
        return 4864 in self.tags

    @property
    def is_pilatus(self) -> bool:
        """Page contains Pilatus tags."""
        return self.software[:8] == 'TVX TIFF' and self.description[:2] == '# '

    @property
    def is_epics(self) -> bool:
        """Page contains EPICS areaDetector tags."""
        return (
            self.description == 'EPICS areaDetector'
            or self.software == 'EPICS areaDetector'
        )

    @property
    def is_tvips(self) -> bool:
        """Page contains TVIPS metadata."""
        return 37706 in self.tags

    @property
    def is_fei(self) -> bool:
        """Page contains FEI_SFEG or FEI_HELIOS tags."""
        return 34680 in self.tags or 34682 in self.tags

    @property
    def is_sem(self) -> bool:
        """Page contains CZ_SEM tag."""
        return 34118 in self.tags

    @property
    def is_svs(self) -> bool:
        """Page contains Aperio metadata."""
        return self.description[:7] == 'Aperio '

    @property
    def is_bif(self) -> bool:
        """Page contains Ventana metadata."""
        try:
            return 700 in self.tags and (
                # avoid reading XMP tag from file at this point
                # b'<iScan' in self.tags[700].value[:4096]
                'Ventana' in self.software
                or self.software[:17] == 'ScanOutputManager'
                or self.description == 'Label Image'
                or self.description == 'Label_Image'
                or self.description == 'Probability_Image'
            )
        except Exception:
            return False

    @property
    def is_scanimage(self) -> bool:
        """Page contains ScanImage metadata."""
        return (
            self.software[:3] == 'SI.'
            or self.description[:6] == 'state.'
            or 'scanimage.SI' in self.description[-256:]
        )

    @property
    def is_indica(self) -> bool:
        """Page contains IndicaLabs metadata."""
        return self.software[:21] == 'IndicaLabsImageWriter'

    @property
    def is_avs(self) -> bool:
        """Page contains Argos AVS XML metadata."""
        try:
            return (
                65000 in self.tags and self.tags.valueof(65000)[:6] == '<Argos'
            )
        except Exception:
            return False

    @property
    def is_qpi(self) -> bool:
        """Page contains PerkinElmer tissue images metadata."""
        # The ImageDescription tag contains XML with a top-level
        # <PerkinElmer-QPI-ImageDescription> element
        return self.software[:15] == 'PerkinElmer-QPI'

    @property
    def is_geotiff(self) -> bool:
        """Page contains GeoTIFF metadata."""
        return 34735 in self.tags  # GeoKeyDirectoryTag

    @property
    def is_gdal(self) -> bool:
        """Page contains GDAL metadata."""
        # startswith '<GDALMetadata>'
        return 42112 in self.tags  # GDAL_METADATA

    @property
    def is_astrotiff(self) -> bool:
        """Page contains AstroTIFF FITS metadata."""
        return (
            self.description[:7] == 'SIMPLE '
            and self.description[-3:] == 'END'
        )

    @property
    def is_streak(self) -> bool:
        """Page contains Hamamatsu streak metadata."""
        return (
            self.description[:1] == '['
            and '],' in self.description[1:32]
            # and self.tags.get(315, '').value[:19] == 'Copyright Hamamatsu'
        )

    @property
    def is_dng(self) -> bool:
        """Page contains DNG metadata."""
        return 50706 in self.tags  # DNGVersion

    @property
    def is_tiffep(self) -> bool:
        """Page contains TIFF/EP metadata."""
        return 37398 in self.tags  # TIFF/EPStandardID

    @property
    def is_sis(self) -> bool:
        """Page contains Olympus SIS metadata."""
        return 33560 in self.tags or 33471 in self.tags

    @property
    def is_ndpi(self) -> bool:
        """Page contains NDPI metadata."""
        return 65420 in self.tags and 271 in self.tags

    @property
    def is_philips(self) -> bool:
        """Page contains Philips DP metadata."""
        return self.software[:10] == 'Philips DP' and self.description[
            -16:
        ].strip().endswith('</DataObject>')

    @property
    def is_eer(self) -> bool:
        """Page contains EER acquisition metadata."""
        return (
            self.parent.is_bigtiff
            and self.compression in {1, 65000, 65001, 65002}
            and 65001 in self.tags
            and self.tags[65001].dtype == 7
        )


@final
class TiffFrame:
    """Lightweight TIFF image file directory (IFD).

    The purpose of TiffFrame is to reduce resource usage and speed up reading
    image data from file compared to TiffPage.
    Properties other than `offset`, `index`, `dataoffsets`, `databytecounts`,
    `subifds`, and `jpegtables` are assumed to be identical with a specified
    TiffPage instance, the keyframe.
    TiffFrame instances have no `tags` property.
    Virtual frames just reference the image data in the file. They may not
    have an IFD structure in the file.

    TiffFrame instances are not thread-safe. All attributes are read-only.

    Parameters:
        parent:
            TiffFile instance to read frame from.
            The file handle position must be at an offset to an IFD structure.
            Only a limited number of tag values are read from file.
        index:
            Index of frame in IFD tree.
        offset:
            Position of frame in file.
        keyframe:
            TiffPage instance with same hash as frame.
        dataoffsets:
            Data offsets of "virtual frame".
        databytecounts:
            Data bytecounts of "virtual frame".

    """

    __slots__ = (
        'parent',
        'offset',
        'dataoffsets',
        'databytecounts',
        'subifds',
        'jpegtables',
        '_keyframe',
        '_index',
    )

    is_mdgel: bool = False
    pages: TiffPages | None = None
    # tags = {}

    parent: TiffFile
    """TiffFile instance frame belongs to."""

    offset: int
    """Position of frame in file."""

    dataoffsets: tuple[int, ...]
    """Positions of strips or tiles in file."""

    databytecounts: tuple[int, ...]
    """Size of strips or tiles in file."""

    subifds: tuple[int, ...] | None
    """Positions of SubIFDs in file."""

    jpegtables: bytes | None
    """JPEG quantization and/or Huffman tables."""

    _keyframe: TiffPage | None
    _index: tuple[int, ...]  # index of frame in IFD tree.

    def __init__(
        self,
        parent: TiffFile,
        /,
        index: int | Sequence[int],
        *,
        offset: int | None = None,
        keyframe: TiffPage | None = None,
        dataoffsets: tuple[int, ...] | None = None,
        databytecounts: tuple[int, ...] | None = None,
    ):
        self._keyframe = None
        self.parent = parent

        self.offset = int(offset) if offset else 0
        self.subifds = None
        self.jpegtables = None
        self.dataoffsets = ()
        self.databytecounts = ()
        if isinstance(index, int):
            self._index = (index,)
        else:
            self._index = tuple(index)

        if dataoffsets is not None and databytecounts is not None:
            # initialize "virtual frame" from offsets and bytecounts
            self.offset = 0 if offset is None else offset
            self.dataoffsets = dataoffsets
            self.databytecounts = databytecounts
            self._keyframe = keyframe
            return

        if offset is None:
            self.offset = parent.filehandle.tell()
        else:
            parent.filehandle.seek(offset)

        if keyframe is None:
            tags = {273, 279, 324, 325, 330, 347}
        elif keyframe.is_contiguous:
            # use databytecounts from keyframe
            tags = {256, 273, 324, 330}
            self.databytecounts = keyframe.databytecounts
        else:
            tags = {256, 273, 279, 324, 325, 330, 347}

        for code, tag in self._gettags(tags):
            if code in {273, 324}:
                self.dataoffsets = tag.value
            elif code in {279, 325}:
                self.databytecounts = tag.value
            elif code == 330:
                self.subifds = tag.value
            elif code == 347:
                self.jpegtables = tag.value
            elif keyframe is None or (
                code == 256 and keyframe.imagewidth != tag.value
            ):
                raise RuntimeError('incompatible keyframe')

        if not self.dataoffsets:
            logger().warning(f'{self!r} is missing required tags')
        elif keyframe is not None and len(self.dataoffsets) != len(
            keyframe.dataoffsets
        ):
            raise RuntimeError('incompatible keyframe')

        if keyframe is not None:
            self.keyframe = keyframe

    def _gettags(
        self,
        codes: Container[int] | None = None,
        /,
        lock: threading.RLock | None = None,
    ) -> list[tuple[int, TiffTag]]:
        """Return list of (code, TiffTag) from file."""
        fh = self.parent.filehandle
        tiff = self.parent.tiff
        unpack = struct.unpack
        rlock: Any = NullContext() if lock is None else lock
        tags = []

        with rlock:
            fh.seek(self.offset)
            try:
                tagno = unpack(tiff.tagnoformat, fh.read(tiff.tagnosize))[0]
                if tagno > 4096:
                    raise ValueError(f'suspicious number of tags {tagno}')
            except Exception as exc:
                raise TiffFileError(
                    f'corrupted tag list @{self.offset}'
                ) from exc

            tagoffset = self.offset + tiff.tagnosize  # fh.tell()
            tagsize = tiff.tagsize
            tagindex = -tagsize
            codeformat = tiff.tagformat1[:2]
            tagbytes = fh.read(tagsize * tagno)

            for _ in range(tagno):
                tagindex += tagsize
                code = unpack(codeformat, tagbytes[tagindex : tagindex + 2])[0]
                if codes and code not in codes:
                    continue
                try:
                    tag = TiffTag.fromfile(
                        self.parent,
                        offset=tagoffset + tagindex,
                        header=tagbytes[tagindex : tagindex + tagsize],
                    )
                except TiffFileError as exc:
                    logger().error(
                        f'{self!r} <TiffTag.fromfile> raised {exc!r:.128}'
                    )
                    continue
                tags.append((code, tag))

        return tags

    def _nextifd(self) -> int:
        """Return offset to next IFD from file."""
        return TiffPage._nextifd(self)  # type: ignore[arg-type]

    def aspage(self) -> TiffPage:
        """Return TiffPage from file.

        Raise ValueError if frame is virtual.

        """
        if self.is_virtual:
            raise ValueError('cannot return virtual frame as page')
        fh = self.parent.filehandle
        closed = fh.closed
        if closed:
            # this is an inefficient resort in case a user calls aspage
            # of a TiffFrame with a closed FileHandle.
            warnings.warn(
                f'{self!r} reading TiffPage from closed file', UserWarning
            )
            fh.open()
        try:
            fh.seek(self.offset)
            page = TiffPage(self.parent, index=self.index)
        finally:
            if closed:
                fh.close()
        return page

    def asarray(self, *args: Any, **kwargs: Any) -> NDArray[Any]:
        """Return image from frame as NumPy array.

        Parameters:
            **kwargs: Arguments passed to :py:meth:`TiffPage.asarray`.

        """
        return TiffPage.asarray(
            self, *args, **kwargs  # type: ignore[arg-type]
        )

    def aszarr(self, **kwargs: Any) -> ZarrTiffStore:
        """Return image from frame as Zarr 2 store.

        Parameters:
            **kwarg: Arguments passed to :py:class:`ZarrTiffStore`.

        """
        return ZarrTiffStore(self, **kwargs)

    def asrgb(self, *args: Any, **kwargs: Any) -> NDArray[Any]:
        """Return image from frame as RGB(A). Work in progress. Do not use.

        :meta private:

        """
        return TiffPage.asrgb(self, *args, **kwargs)  # type: ignore[arg-type]

    def segments(self, *args: Any, **kwargs: Any) -> Iterator[
        tuple[
            NDArray[Any] | None,
            tuple[int, int, int, int, int],
            tuple[int, int, int, int],
        ]
    ]:
        """Return iterator over decoded tiles or strips.

        Parameters:
            **kwargs: Arguments passed to :py:meth:`TiffPage.segments`.

        :meta private:

        """
        return TiffPage.segments(
            self, *args, **kwargs  # type: ignore[arg-type]
        )

    @property
    def index(self) -> int:
        """Index of frame in IFD chain."""
        return self._index[-1]

    @property
    def treeindex(self) -> tuple[int, ...]:
        """Index of frame in IFD tree."""
        return self._index

    @property
    def keyframe(self) -> TiffPage | None:
        """TiffPage with same properties as this frame."""
        return self._keyframe

    @keyframe.setter
    def keyframe(self, keyframe: TiffPage, /) -> None:
        if self._keyframe == keyframe:
            return
        if self._keyframe is not None:
            raise RuntimeError('cannot reset keyframe')
        if len(self.dataoffsets) != len(keyframe.dataoffsets):
            raise RuntimeError('incompatible keyframe')
        if keyframe.is_contiguous:
            self.databytecounts = keyframe.databytecounts
        self._keyframe = keyframe

    @property
    def is_frame(self) -> bool:
        """Object is :py:class:`TiffFrame` instance."""
        return True

    @property
    def is_virtual(self) -> bool:
        """Frame does not have IFD structure in file."""
        return self.offset <= 0

    @property
    def is_subifd(self) -> bool:
        """Frame is SubIFD of another page."""
        return len(self._index) > 1

    @property
    def is_final(self) -> bool:
        assert self._keyframe is not None
        return self._keyframe.is_final

    @property
    def is_contiguous(self) -> bool:
        assert self._keyframe is not None
        return self._keyframe.is_contiguous

    @property
    def is_memmappable(self) -> bool:
        assert self._keyframe is not None
        return self._keyframe.is_memmappable

    @property
    def hash(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.hash

    @property
    def shape(self) -> tuple[int, ...]:
        assert self._keyframe is not None
        return self._keyframe.shape

    @property
    def shaped(self) -> tuple[int, int, int, int, int]:
        assert self._keyframe is not None
        return self._keyframe.shaped

    @property
    def chunks(self) -> tuple[int, ...]:
        assert self._keyframe is not None
        return self._keyframe.chunks

    @property
    def chunked(self) -> tuple[int, ...]:
        assert self._keyframe is not None
        return self._keyframe.chunked

    @property
    def tile(self) -> tuple[int, ...] | None:
        assert self._keyframe is not None
        return self._keyframe.tile

    @property
    def name(self) -> str:
        index = self._index if len(self._index) > 1 else self._index[0]
        return f'TiffFrame {index}'

    @property
    def ndim(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.ndim

    @property
    def dims(self) -> tuple[str, ...]:
        assert self._keyframe is not None
        return self._keyframe.dims

    @property
    def sizes(self) -> dict[str, int]:
        assert self._keyframe is not None
        return self._keyframe.sizes

    @property
    def coords(self) -> dict[str, NDArray[Any]]:
        assert self._keyframe is not None
        return self._keyframe.coords

    @property
    def size(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.size

    @property
    def nbytes(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.nbytes

    @property
    def dtype(self) -> numpy.dtype[Any] | None:
        assert self._keyframe is not None
        return self._keyframe.dtype

    @property
    def axes(self) -> str:
        assert self._keyframe is not None
        return self._keyframe.axes

    def get_resolution(
        self,
        unit: RESUNIT | int | None = None,
        scale: float | int | None = None,
    ) -> tuple[int | float, int | float]:
        assert self._keyframe is not None
        return self._keyframe.get_resolution(unit, scale)

    @property
    def resolution(self) -> tuple[float, float]:
        assert self._keyframe is not None
        return self._keyframe.resolution

    @property
    def resolutionunit(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.resolutionunit

    @property
    def datetime(self) -> DateTime | None:
        # TODO: TiffFrame.datetime can differ from TiffPage.datetime?
        assert self._keyframe is not None
        return self._keyframe.datetime

    @property
    def compression(self) -> int:
        assert self._keyframe is not None
        return self._keyframe.compression

    @property
    def decode(
        self,
    ) -> Callable[
        ...,
        tuple[
            NDArray[Any] | None,
            tuple[int, int, int, int, int],
            tuple[int, int, int, int],
        ],
    ]:
        assert self._keyframe is not None
        return self._keyframe.decode

    def __repr__(self) -> str:
        index = self._index if len(self._index) > 1 else self._index[0]
        return f'<tifffile.TiffFrame {index} @{self.offset}>'

    def __str__(self) -> str:
        return self._str()

    def _str(self, detail: int = 0, width: int = 79) -> str:
        """Return string containing information about TiffFrame."""
        if self._keyframe is None:
            info = ''
            kf = None
        else:
            info = '  '.join(
                s
                for s in (
                    'x'.join(str(i) for i in self.shape),
                    str(self.dtype),
                )
            )
            kf = self._keyframe._str(width=width - 11)
        if detail > 3:
            of = pformat(self.dataoffsets, width=width - 9, height=detail - 3)
            bc = pformat(
                self.databytecounts, width=width - 13, height=detail - 3
            )
            info = f'\n Keyframe {kf}\n Offsets {of}\n Bytecounts {bc}'
        index = self._index if len(self._index) > 1 else self._index[0]
        return f'TiffFrame {index} @{self.offset}  {info}'


@final
class TiffPages(Sequence[TiffPage | TiffFrame]):
    """Sequence of TIFF image file directories (IFD chain).

    TiffPages instances have a state, such as a cache and keyframe, and are not
    thread-safe. All attributes are read-only.

    Parameters:
        arg:
            If a *TiffFile*, the file position must be at offset to offset to
            TiffPage.
            If a *TiffPage* or *TiffFrame*, page offsets are read from the
            SubIFDs tag.
            Only the first page is initially read from the file.
        index:
            Position of IFD chain in IFD tree.

    """

    parent: TiffFile | None = None
    """TiffFile instance pages belongs to."""

    _pages: list[TiffPage | TiffFrame | int]  # list of pages
    _keyframe: TiffPage | None
    _tiffpage: type[TiffPage] | type[TiffFrame]  # class used for reading pages
    _indexed: bool
    _cached: bool
    _cache: bool
    _offset: int
    _nextpageoffset: int | None
    _index: tuple[int, ...] | None

    def __init__(
        self,
        arg: TiffFile | TiffPage | TiffFrame,
        /,
        *,
        index: Sequence[int] | int | None = None,
    ) -> None:
        offset: int
        self.parent = None
        self._pages = []  # cache of TiffPages, TiffFrames, or their offsets
        self._indexed = False  # True if offsets to all pages were read
        self._cached = False  # True if all pages were read into cache
        self._tiffpage = TiffPage  # class used for reading pages
        self._keyframe = None  # page that is currently used as keyframe
        self._cache = False  # do not cache frames or pages (if not keyframe)
        self._offset = 0
        self._nextpageoffset = None

        if index is None:
            self._index = None
        elif isinstance(index, (int, numpy.integer)):
            self._index = (int(index),)
        else:
            self._index = tuple(index)

        if isinstance(arg, TiffFile):
            # read offset to first page from current file position
            self.parent = arg
            fh = self.parent.filehandle
            self._nextpageoffset = fh.tell()
            offset = struct.unpack(
                self.parent.tiff.offsetformat,
                fh.read(self.parent.tiff.offsetsize),
            )[0]
            if offset == 0:
                logger().warning(f'{arg!r} contains no pages')
                self._indexed = True
                return
        elif arg.subifds is not None:
            # use offsets from SubIFDs tag
            offsets = arg.subifds
            self.parent = arg.parent
            fh = self.parent.filehandle
            if len(offsets) == 0 or offsets[0] == 0:
                logger().warning(f'{arg!r} contains invalid SubIFDs')
                self._indexed = True
                return
            offset = offsets[0]
        else:
            self._indexed = True
            return

        self._offset = offset
        if offset >= fh.size:
            logger().warning(
                f'{self!r} invalid offset to first page {offset!r}'
            )
            self._indexed = True
            return

        pageindex: int | tuple[int, ...] = (
            0 if self._index is None else self._index + (0,)
        )

        # read and cache first page
        fh.seek(offset)
        page = TiffPage(self.parent, index=pageindex)
        self._pages.append(page)
        self._keyframe = page
        if self._nextpageoffset is None:
            # offsets from SubIFDs tag
            self._pages.extend(offsets[1:])
            self._indexed = True
            self._cached = True

    @property
    def pages(self) -> list[TiffPage | TiffFrame | int]:
        """Deprecated. Use the TiffPages sequence interface.

        :meta private:

        """
        warnings.warn(
            '<tifffile.TiffPages.pages> is deprecated since 2024.5.22. '
            'Use the TiffPages sequence interface.',
            DeprecationWarning,
            stacklevel=2,
        )
        return self._pages

    @property
    def first(self) -> TiffPage:
        """First page as TiffPage if exists, else raise IndexError."""
        return cast(TiffPage, self._pages[0])

    @property
    def is_multipage(self) -> bool:
        """IFD chain contains more than one page."""
        try:
            self._seek(1)
            return True
        except IndexError:
            return False

    @property
    def cache(self) -> bool:
        """Pages and frames are being cached.

        When set to *False*, the cache is cleared.

        """
        return self._cache

    @cache.setter
    def cache(self, value: bool, /) -> None:
        value = bool(value)
        if self._cache and not value:
            self._clear()
        self._cache = value

    @property
    def useframes(self) -> bool:
        """Use TiffFrame (True) or TiffPage (False)."""
        return self._tiffpage == TiffFrame

    @useframes.setter
    def useframes(self, value: bool, /) -> None:
        self._tiffpage = TiffFrame if value else TiffPage

    @property
    def keyframe(self) -> TiffPage | None:
        """TiffPage used as keyframe for new TiffFrames."""
        return self._keyframe

    def set_keyframe(self, index: int, /) -> None:
        """Set keyframe to TiffPage specified by `index`.

        If not found in the cache, the TiffPage at `index` is loaded from file
        and added to the cache.

        """
        if not isinstance(index, (int, numpy.integer)):
            raise TypeError(f'indices must be integers, not {type(index)}')
        index = int(index)
        if index < 0:
            index %= len(self)
        if self._keyframe is not None and self._keyframe.index == index:
            return
        if index == 0:
            self._keyframe = cast(TiffPage, self._pages[0])
            return
        if self._indexed or index < len(self._pages):
            page = self._pages[index]
            if isinstance(page, TiffPage):
                self._keyframe = page
                return
            if isinstance(page, TiffFrame):
                # remove existing TiffFrame
                self._pages[index] = page.offset
        # load TiffPage from file
        tiffpage = self._tiffpage
        self._tiffpage = TiffPage
        try:
            self._keyframe = cast(TiffPage, self._getitem(index))
        finally:
            self._tiffpage = tiffpage
        # always cache keyframes
        self._pages[index] = self._keyframe

    @property
    def next_page_offset(self) -> int | None:
        """Offset where offset to new page can be stored."""
        if not self._indexed:
            self._seek(-1)
        return self._nextpageoffset

    def get(
        self,
        key: int,
        /,
        default: TiffPage | TiffFrame | None = None,
        *,
        validate: int = 0,
        cache: bool = False,
        aspage: bool = True,
    ) -> TiffPage | TiffFrame:
        """Return specified page from cache or file.

        The specified TiffPage or TiffFrame is read from file if it is not
        found in the cache.

        Parameters:
            key:
                Index of requested page in IFD chain.
            default:
                Page or frame to return if key is out of bounds.
                By default, an IndexError is raised if key is out of bounds.
            validate:
                If non-zero, raise RuntimeError if value does not match hash
                of TiffPage or TiffFrame.
            cache:
                Store returned page in cache for future use.
            aspage:
                Return TiffPage instance.

        """
        try:
            return self._getitem(
                key, validate=validate, cache=cache, aspage=aspage
            )
        except IndexError:
            if default is None:
                raise
        return default

    def _load(self, keyframe: TiffPage | bool | None = True, /) -> None:
        """Read all remaining pages from file."""
        assert self.parent is not None
        if self._cached:
            return
        pages = self._pages
        if not pages:
            return
        if not self._indexed:
            self._seek(-1)
        if not self._cache:
            return
        fh = self.parent.filehandle
        if keyframe is not None:
            keyframe = self._keyframe
        for i, page in enumerate(pages):
            if isinstance(page, (int, numpy.integer)):
                pageindex: int | tuple[int, ...] = (
                    i if self._index is None else self._index + (i,)
                )
                fh.seek(page)
                page = self._tiffpage(
                    self.parent, index=pageindex, keyframe=keyframe
                )
                pages[i] = page
        self._cached = True

    def _load_virtual_frames(self) -> None:
        """Calculate virtual TiffFrames."""
        assert self.parent is not None
        pages = self._pages
        try:
            if len(pages) > 1:
                raise ValueError('pages already loaded')
            page = cast(TiffPage, pages[0])
            if not page.is_contiguous:
                raise ValueError('data not contiguous')
            self._seek(4)
            # following pages are int
            delta = cast(int, pages[2]) - cast(int, pages[1])
            if (
                cast(int, pages[3]) - cast(int, pages[2]) != delta
                or cast(int, pages[4]) - cast(int, pages[3]) != delta
            ):
                raise ValueError('page offsets not equidistant')
            page1 = self._getitem(1, validate=page.hash)
            offsetoffset = page1.dataoffsets[0] - page1.offset
            if offsetoffset < 0 or offsetoffset > delta:
                raise ValueError('page offsets not equidistant')
            pages = [page, page1]
            filesize = self.parent.filehandle.size - delta

            for index, offset in enumerate(
                range(page1.offset + delta, filesize, delta)
            ):
                index += 2
                d = index * delta
                dataoffsets = tuple(i + d for i in page.dataoffsets)
                offset_or_none = offset if offset < 2**31 - 1 else None
                pages.append(
                    TiffFrame(
                        page.parent,
                        index=(
                            index
                            if self._index is None
                            else self._index + (index,)
                        ),
                        offset=offset_or_none,
                        dataoffsets=dataoffsets,
                        databytecounts=page.databytecounts,
                        keyframe=page,
                    )
                )
            self._pages = pages
            self._cache = True
            self._cached = True
            self._indexed = True
        except Exception as exc:
            if self.parent.filehandle.size >= 2147483648:
                logger().warning(
                    f'{self!r} <_load_virtual_frames> raised {exc!r:.128}'
                )

    def _clear(self, fully: bool = True, /) -> None:
        """Delete all but first page from cache. Set keyframe to first page."""
        pages = self._pages
        if not pages:
            return
        self._keyframe = cast(TiffPage, pages[0])
        if fully:
            # delete all but first TiffPage/TiffFrame
            for i, page in enumerate(pages[1:]):
                if not isinstance(page, int) and page.offset is not None:
                    pages[i + 1] = page.offset
        else:
            # delete only TiffFrames
            for i, page in enumerate(pages):
                if isinstance(page, TiffFrame) and page.offset is not None:
                    pages[i] = page.offset
        self._cached = False

    def _seek(self, index: int, /) -> int:
        """Seek file to offset of page specified by index and return offset."""
        assert self.parent is not None

        pages = self._pages
        lenpages = len(pages)
        if lenpages == 0:
            raise IndexError('index out of range')

        fh = self.parent.filehandle
        if fh.closed:
            raise ValueError('seek of closed file')

        if self._indexed or 0 <= index < lenpages:
            page = pages[index]
            offset = page if isinstance(page, int) else page.offset
            return fh.seek(offset)

        tiff = self.parent.tiff
        offsetformat = tiff.offsetformat
        offsetsize = tiff.offsetsize
        tagnoformat = tiff.tagnoformat
        tagnosize = tiff.tagnosize
        tagsize = tiff.tagsize
        unpack = struct.unpack

        page = pages[-1]
        offset = page if isinstance(page, int) else page.offset

        while lenpages < 2**32:
            # read offsets to pages from file until index is reached
            fh.seek(offset)
            # skip tags
            try:
                tagno = int(unpack(tagnoformat, fh.read(tagnosize))[0])
                if tagno > 4096:
                    raise TiffFileError(f'suspicious number of tags {tagno}')
            except Exception as exc:
                logger().error(
                    f'{self!r} corrupted tag list of page '
                    f'{lenpages} @{offset} raised {exc!r:.128}',
                )
                del pages[-1]
                lenpages -= 1
                self._indexed = True
                break
            self._nextpageoffset = offset + tagnosize + tagno * tagsize
            fh.seek(self._nextpageoffset)

            # read offset to next page
            try:
                offset = int(unpack(offsetformat, fh.read(offsetsize))[0])
            except Exception as exc:
                logger().error(
                    f'{self!r} invalid offset to page '
                    f'{lenpages + 1} @{self._nextpageoffset} '
                    f'raised {exc!r:.128}'
                )
                self._indexed = True
                break
            if offset == 0:
                self._indexed = True
                break
            if offset >= fh.size:
                logger().error(f'{self!r} invalid page offset {offset!r}')
                self._indexed = True
                break

            pages.append(offset)
            lenpages += 1
            if 0 <= index < lenpages:
                break

            # detect some circular references
            if lenpages == 100:
                for i, p in enumerate(pages[:-1]):
                    if offset == (p if isinstance(p, int) else p.offset):
                        index = i
                        self._pages = pages[: i + 1]
                        self._indexed = True
                        logger().error(
                            f'{self!r} invalid circular reference to IFD '
                            f'{i} at {offset=}'
                        )
                        break

        if index >= lenpages:
            raise IndexError('index out of range')

        page = pages[index]
        return fh.seek(page if isinstance(page, int) else page.offset)

    def _getlist(
        self,
        key: int | slice | Iterable[int] | None = None,
        /,
        useframes: bool = True,
        validate: bool = True,
    ) -> list[TiffPage | TiffFrame]:
        """Return specified pages as list of TiffPages or TiffFrames.

        The first item is a TiffPage, and is used as a keyframe for
        following TiffFrames.

        """
        getitem = self._getitem
        _useframes = self.useframes

        if key is None:
            key = iter(range(len(self)))
        elif isinstance(key, (int, numpy.integer)):
            # return single TiffPage
            key = int(key)
            self.useframes = False
            if key == 0:
                return [self.first]
            try:
                return [getitem(key)]
            finally:
                self.useframes = _useframes
        elif isinstance(key, slice):
            start, stop, _ = key.indices(2**31 - 1)
            if not self._indexed and max(stop, start) > len(self._pages):
                self._seek(-1)
            key = iter(range(*key.indices(len(self._pages))))
        elif isinstance(key, Iterable):
            key = iter(key)
        else:
            raise TypeError(
                f'key must be an integer, slice, or iterable, not {type(key)}'
            )

        # use first page as keyframe
        assert self._keyframe is not None
        keyframe = self._keyframe
        self.set_keyframe(next(key))
        validhash = self._keyframe.hash if validate else 0
        if useframes:
            self.useframes = True
        try:
            pages = [getitem(i, validate=validhash) for i in key]
            pages.insert(0, self._keyframe)
        finally:
            # restore state
            self._keyframe = keyframe
            if useframes:
                self.useframes = _useframes
        return pages

    def _getitem(
        self,
        key: int,
        /,
        *,
        validate: int = 0,  # hash
        cache: bool = False,
        aspage: bool = False,
    ) -> TiffPage | TiffFrame:
        """Return specified page from cache or file."""
        assert self.parent is not None
        key = int(key)
        pages = self._pages

        if key < 0:
            key %= len(self)
        elif self._indexed and key >= len(pages):
            raise IndexError(f'index {key} out of range({len(pages)})')

        tiffpage = TiffPage if aspage else self._tiffpage

        if key < len(pages):
            page = pages[key]
            if self._cache and not aspage:
                if not isinstance(page, (int, numpy.integer)):
                    if validate and validate != page.hash:
                        raise RuntimeError('page hash mismatch')
                    return page
            elif isinstance(page, (TiffPage, tiffpage)):
                # page is not an int
                if (
                    validate
                    and validate != page.hash  # type: ignore[union-attr]
                ):
                    raise RuntimeError('page hash mismatch')
                return page  # type: ignore[return-value]

        pageindex: int | tuple[int, ...] = (
            key if self._index is None else self._index + (key,)
        )
        self._seek(key)
        page = tiffpage(self.parent, index=pageindex, keyframe=self._keyframe)
        assert isinstance(page, (TiffPage, TiffFrame))
        if validate and validate != page.hash:
            raise RuntimeError('page hash mismatch')
        if self._cache or cache:
            pages[key] = page
        return page

    @overload
    def __getitem__(self, key: int, /) -> TiffPage | TiffFrame: ...

    @overload
    def __getitem__(
        self, key: slice | Iterable[int], /
    ) -> list[TiffPage | TiffFrame]: ...

    def __getitem__(
        self, key: int | slice | Iterable[int], /
    ) -> TiffPage | TiffFrame | list[TiffPage | TiffFrame]:
        pages = self._pages
        getitem = self._getitem

        if isinstance(key, (int, numpy.integer)):
            key = int(key)
            if key == 0:
                return cast(TiffPage, pages[key])
            return getitem(key)

        if isinstance(key, slice):
            start, stop, _ = key.indices(2**31 - 1)
            if not self._indexed and max(stop, start) > len(pages):
                self._seek(-1)
            return [getitem(i) for i in range(*key.indices(len(pages)))]

        if isinstance(key, Iterable):
            return [getitem(k) for k in key]

        raise TypeError('key must be an integer, slice, or iterable')

    def __iter__(self) -> Iterator[TiffPage | TiffFrame]:
        i = 0
        while True:
            try:
                yield self._getitem(i)
                i += 1
            except IndexError:
                break
        if self._cache:
            self._cached = True

    def __bool__(self) -> bool:
        """Return True if file contains any pages."""
        return len(self._pages) > 0

    def __len__(self) -> int:
        """Return number of pages in file."""
        if not self._indexed:
            self._seek(-1)
        return len(self._pages)

    def __repr__(self) -> str:
        return f'<tifffile.TiffPages @{self._offset}>'


@final
class TiffTag:
    """TIFF tag structure.

    TiffTag instances are not thread-safe. All attributes are read-only.

    Parameters:
        parent:
            TIFF file tag belongs to.
        offset:
            Position of tag structure in file.
        code:
            Decimal code of tag.
        dtype:
            Data type of tag value item.
        count:
            Number of items in tag value.
        valueoffset:
            Position of tag value in file.

    """

    __slots__ = (
        'parent',
        'offset',
        'code',
        'dtype',
        'count',
        '_value',
        'valueoffset',
    )

    parent: TiffFile | TiffWriter
    """TIFF file tag belongs to."""

    offset: int
    """Position of tag structure in file."""

    code: int
    """Decimal code of tag."""

    dtype: int
    """:py:class:`DATATYPE` of tag value item."""

    count: int
    """Number of items in tag value."""

    valueoffset: int
    """Position of tag value in file."""

    _value: Any

    def __init__(
        self,
        parent: TiffFile | TiffWriter,
        offset: int,
        code: int,
        dtype: DATATYPE | int,
        count: int,
        value: Any,
        valueoffset: int,
        /,
    ) -> None:
        self.parent = parent
        self.offset = int(offset)
        self.code = int(code)
        self.count = int(count)
        self._value = value
        self.valueoffset = valueoffset
        try:
            self.dtype = DATATYPE(dtype)
        except ValueError:
            self.dtype = int(dtype)

    @classmethod
    def fromfile(
        cls,
        parent: TiffFile,
        /,
        *,
        offset: int | None = None,
        header: bytes | None = None,
        validate: bool = True,
    ) -> TiffTag:
        """Return TiffTag instance from file.

        Parameters:
            parent:
                TiffFile instance tag is read from.
            offset:
                Position of tag structure in file.
                The default is the position of the file handle.
            header:
                Tag structure as bytes.
                The default is read from the file.
            validate:
                Raise TiffFileError if data type or value offset are invalid.

        Raises:
            TiffFileError:
                Data type or value offset are invalid and `validate` is *True*.

        """
        tiff = parent.tiff
        fh = parent.filehandle

        if header is None:
            if offset is None:
                offset = fh.tell()
            else:
                fh.seek(offset)
            header = fh.read(tiff.tagsize)
        elif offset is None:
            offset = fh.tell()

        valueoffset = offset + tiff.tagsize - tiff.tagoffsetthreshold
        code, dtype, count, value = struct.unpack(
            tiff.tagformat1 + tiff.tagformat2[1:], header
        )

        try:
            valueformat = TIFF.DATA_FORMATS[dtype]
        except KeyError as exc:
            msg = (
                f'<tifffile.TiffTag {code} @{offset}> '
                f'invalid data type {dtype!r}'
            )
            if validate:
                raise TiffFileError(msg) from exc
            logger().error(msg)
            return cls(parent, offset, code, dtype, count, None, 0)

        valuesize = count * struct.calcsize(valueformat)
        if (
            valuesize > tiff.tagoffsetthreshold
            or code in TIFF.TAG_READERS  # TODO: only works with offsets?
        ):
            valueoffset = struct.unpack(tiff.offsetformat, value)[0]
            if validate and code in TIFF.TAG_LOAD:
                value = TiffTag._read_value(
                    parent, offset, code, dtype, count, valueoffset
                )
            elif valueoffset < 8 or valueoffset + valuesize > fh.size:
                msg = (
                    f'<tifffile.TiffTag {code} @{offset}> '
                    f'invalid value offset {valueoffset}'
                )
                if validate:
                    raise TiffFileError(msg)
                logger().warning(msg)
                value = None
            elif code in TIFF.TAG_LOAD:
                value = TiffTag._read_value(
                    parent, offset, code, dtype, count, valueoffset
                )
            else:
                value = None
        elif dtype in {1, 2, 7}:
            # BYTES, ASCII, UNDEFINED
            value = value[:valuesize]
        elif (
            tiff.is_ndpi
            and count == 1
            and dtype in {4, 9, 13}
            and value[4:] != b'\x00\x00\x00\x00'
        ):
            # NDPI IFD or LONG, for example, in StripOffsets or StripByteCounts
            value = struct.unpack('<Q', value)
        else:
            fmt = (
                f'{tiff.byteorder}'
                f'{count * int(valueformat[0])}'
                f'{valueformat[1]}'
            )
            value = struct.unpack(fmt, value[:valuesize])

        value = TiffTag._process_value(value, code, dtype, offset)

        return cls(parent, offset, code, dtype, count, value, valueoffset)

    @staticmethod
    def _read_value(
        parent: TiffFile | TiffWriter,
        offset: int,
        code: int,
        dtype: int,
        count: int,
        valueoffset: int,
        /,
    ) -> Any:
        """Read tag value from file."""
        try:
            valueformat = TIFF.DATA_FORMATS[dtype]
        except KeyError as exc:
            raise TiffFileError(
                f'<tifffile.TiffTag {code} @{offset}> '
                f'invalid data type {dtype!r}'
            ) from exc

        fh = parent.filehandle
        byteorder = parent.tiff.byteorder
        offsetsize = parent.tiff.offsetsize

        valuesize = count * struct.calcsize(valueformat)
        if valueoffset < 8 or valueoffset + valuesize > fh.size:
            raise TiffFileError(
                f'<tifffile.TiffTag {code} @{offset}> '
                f'invalid value offset {valueoffset}'
            )
        # if valueoffset % 2:
        #     logger().warning(
        #         f'<tifffile.TiffTag {code} @{offset}> '
        #         'value does not begin on word boundary'
        #     )

        fh.seek(valueoffset)
        if code in TIFF.TAG_READERS:
            readfunc = TIFF.TAG_READERS[code]
            try:
                value = readfunc(fh, byteorder, dtype, count, offsetsize)
            except Exception as exc:
                logger().warning(
                    f'<tifffile.TiffTag {code} @{offset}> raised {exc!r:.128}'
                )
            else:
                return value

        if dtype in {1, 2, 7}:
            # BYTES, ASCII, UNDEFINED
            value = fh.read(valuesize)
            if len(value) != valuesize:
                logger().warning(
                    f'<tifffile.TiffTag {code} @{offset}> '
                    'could not read all values'
                )
        elif code not in TIFF.TAG_TUPLE and count > 1024:
            value = read_numpy(fh, byteorder, dtype, count, offsetsize)
        else:
            value = struct.unpack(
                f'{byteorder}{count * int(valueformat[0])}{valueformat[1]}',
                fh.read(valuesize),
            )
        return value

    @staticmethod
    def _process_value(
        value: Any, code: int, dtype: int, offset: int, /
    ) -> Any:
        """Process tag value."""
        if (
            value is None
            or dtype == 1  # BYTE
            or dtype == 7  # UNDEFINED
            or code in TIFF.TAG_READERS
            or not isinstance(value, (bytes, str, tuple))
        ):
            return value

        if dtype == 2:
            # TIFF ASCII fields can contain multiple strings,
            #   each terminated with a NUL
            value = value.rstrip(b'\x00')
            try:
                value = value.decode('utf-8').strip()
            except UnicodeDecodeError:
                try:
                    value = value.decode('cp1252').strip()
                except UnicodeDecodeError as exc:
                    logger().warning(
                        f'<tifffile.TiffTag {code} @{offset}> '
                        f'coercing invalid ASCII to bytes, due to {exc!r:.128}'
                    )
            return value

        if code in TIFF.TAG_ENUM:
            t = TIFF.TAG_ENUM[code]
            try:
                value = tuple(t(v) for v in value)
            except ValueError as exc:
                if code not in {259, 317}:  # ignore compression/predictor
                    logger().warning(
                        f'<tifffile.TiffTag {code} @{offset}> '
                        f'raised {exc!r:.128}'
                    )

        if len(value) == 1 and code not in TIFF.TAG_TUPLE:
            value = value[0]

        return value

    @property
    def value(self) -> Any:
        """Value of tag, delay-loaded from file if necessary."""
        if self._value is None:
            # print(
            #     f'_read_value {self.code} {TIFF.TAGS.get(self.code)} '
            #     f'{self.dtype}[{self.count}] @{self.valueoffset} '
            # )
            fh = self.parent.filehandle
            with fh.lock:
                closed = fh.closed
                if closed:
                    # this is an inefficient resort in case a user delay loads
                    # tag values from a TiffPage with a closed FileHandle.
                    warnings.warn(
                        f'{self!r} reading value from closed file', UserWarning
                    )
                    fh.open()
                try:
                    value = TiffTag._read_value(
                        self.parent,
                        self.offset,
                        self.code,
                        self.dtype,
                        self.count,
                        self.valueoffset,
                    )
                finally:
                    if closed:
                        fh.close()
            self._value = TiffTag._process_value(
                value,
                self.code,
                self.dtype,
                self.offset,
            )
        return self._value

    @value.setter
    def value(self, value: Any, /) -> None:
        self._value = value

    @property
    def dtype_name(self) -> str:
        """Name of data type of tag value."""
        try:
            return self.dtype.name  # type: ignore[attr-defined]
        except AttributeError:
            return f'TYPE{self.dtype}'

    @property
    def name(self) -> str:
        """Name of tag from :py:attr:`_TIFF.TAGS` registry."""
        return TIFF.TAGS.get(self.code, str(self.code))

    @property
    def dataformat(self) -> str:
        """Data type as `struct.pack` format."""
        return TIFF.DATA_FORMATS[self.dtype]

    @property
    def valuebytecount(self) -> int:
        """Number of bytes of tag value in file."""
        return self.count * struct.calcsize(TIFF.DATA_FORMATS[self.dtype])

    def astuple(self) -> TagTuple:
        """Return tag code, dtype, count, and encoded value.

        The encoded value is read from file if necessary.

        """
        if isinstance(self.value, bytes):
            value = self.value
        else:
            tiff = self.parent.tiff
            dataformat = TIFF.DATA_FORMATS[self.dtype]
            count = self.count * int(dataformat[0])
            fmt = f'{tiff.byteorder}{count}{dataformat[1]}'
            try:
                if self.dtype == 2:
                    # ASCII
                    value = struct.pack(fmt, self.value.encode('ascii'))
                    if len(value) != count:
                        raise ValueError
                elif count == 1 and not isinstance(self.value, tuple):
                    value = struct.pack(fmt, self.value)
                else:
                    value = struct.pack(fmt, *self.value)
            except Exception as exc:
                if tiff.is_ndpi and count == 1:
                    raise ValueError(
                        'cannot pack 64-bit NDPI value to 32-bit dtype'
                    ) from exc
                fh = self.parent.filehandle
                pos = fh.tell()
                fh.seek(self.valueoffset)
                value = fh.read(struct.calcsize(fmt))
                fh.seek(pos)
        return self.code, int(self.dtype), self.count, value, True

    def overwrite(
        self,
        value: Any,
        /,
        *,
        dtype: DATATYPE | int | str | None = None,
        erase: bool = True,
    ) -> TiffTag:
        """Write new tag value to file and return new TiffTag instance.

        Warning: changing tag values in TIFF files might result in corrupted
        files or have unexpected side effects.

        The packed value is appended to the file if it is longer than the
        old value. The file position is left where it was.

        Overwriting tag values in NDPI files > 4 GB is only supported if
        single integer values and new offsets do not exceed the 32-bit range.

        Parameters:
            value:
                New tag value to write.
                Must be compatible with the `struct.pack` formats corresponding
                to the tag's data type.
            dtype:
                New tag data type. By default, the data type is not changed.
            erase:
                Overwrite previous tag values in file with zeros.

        Raises:
            struct.error:
                Value is not compatible with dtype or new offset exceeds
                TIFF size limit.
            ValueError:
                Invalid value or dtype, or old integer value in NDPI files
                exceeds 32-bit range.

        """
        if self.offset < 8 or self.valueoffset < 8:
            raise ValueError(f'cannot rewrite tag at offset {self.offset} < 8')

        if hasattr(value, 'filehandle'):
            # passing a TiffFile instance is deprecated and no longer required
            # since 2021.7.30
            raise TypeError(
                'TiffTag.overwrite got an unexpected TiffFile instance '
                'as first argument'
            )

        fh = self.parent.filehandle
        tiff = self.parent.tiff
        if tiff.is_ndpi:
            # only support files < 4GB
            if self.count == 1 and self.dtype in {4, 13}:
                if isinstance(self.value, tuple):
                    v = self.value[0]
                else:
                    v = self.value
                if v > 4294967295:
                    raise ValueError('cannot patch NDPI > 4 GB files')
            tiff = TIFF.CLASSIC_LE

        if value is None:
            value = b''
        if dtype is None:
            dtype = self.dtype
        elif isinstance(dtype, str):
            if len(dtype) > 1 and dtype[0] in '<>|=':
                dtype = dtype[1:]
            try:
                dtype = TIFF.DATA_DTYPES[dtype]
            except KeyError as exc:
                raise ValueError(f'unknown data type {dtype!r}') from exc
        else:
            dtype = enumarg(DATATYPE, dtype)

        packedvalue: bytes | None = None
        dataformat: str
        try:
            dataformat = TIFF.DATA_FORMATS[dtype]
        except KeyError as exc:
            raise ValueError(f'unknown data type {dtype!r}') from exc

        if dtype == 2:
            # strings
            if isinstance(value, str):
                # enforce 7-bit ASCII on Unicode strings
                try:
                    value = value.encode('ascii')
                except UnicodeEncodeError as exc:
                    raise ValueError(
                        'TIFF strings must be 7-bit ASCII'
                    ) from exc
            elif not isinstance(value, bytes):
                raise ValueError('TIFF strings must be 7-bit ASCII')
            if len(value) == 0 or value[-1:] != b'\x00':
                value += b'\x00'
            count = len(value)
            value = (value,)

        elif isinstance(value, bytes):
            # pre-packed binary data
            dtsize = struct.calcsize(dataformat)
            if len(value) % dtsize:
                raise ValueError('invalid packed binary data')
            count = len(value) // dtsize
            packedvalue = value
            value = (value,)

        else:
            try:
                count = len(value)
            except TypeError:
                value = (value,)
                count = 1
            if dtype in {5, 10}:
                if count < 2 or count % 2:
                    raise ValueError('invalid RATIONAL value')
                count //= 2  # rational

        if packedvalue is None:
            packedvalue = struct.pack(
                f'{tiff.byteorder}{count * int(dataformat[0])}{dataformat[1]}',
                *value,
            )
        newsize = len(packedvalue)
        oldsize = self.count * struct.calcsize(TIFF.DATA_FORMATS[self.dtype])
        valueoffset = self.valueoffset

        pos = fh.tell()
        try:
            if dtype != self.dtype:
                # rewrite data type
                fh.seek(self.offset + 2)
                fh.write(struct.pack(tiff.byteorder + 'H', dtype))

            if oldsize <= tiff.tagoffsetthreshold:
                if newsize <= tiff.tagoffsetthreshold:
                    # inline -> inline: overwrite
                    fh.seek(self.offset + 4)
                    fh.write(struct.pack(tiff.tagformat2, count, packedvalue))
                else:
                    # inline -> separate: append to file
                    fh.seek(0, os.SEEK_END)
                    valueoffset = fh.tell()
                    if valueoffset % 2:
                        # value offset must begin on a word boundary
                        fh.write(b'\x00')
                        valueoffset += 1
                    # write new offset
                    fh.seek(self.offset + 4)
                    fh.write(
                        struct.pack(
                            tiff.tagformat2,
                            count,
                            struct.pack(tiff.offsetformat, valueoffset),
                        )
                    )
                    # write new value
                    fh.seek(valueoffset)
                    fh.write(packedvalue)

            elif newsize <= tiff.tagoffsetthreshold:
                # separate -> inline: erase old value
                valueoffset = (
                    self.offset + 4 + struct.calcsize(tiff.tagformat2[:2])
                )
                fh.seek(self.offset + 4)
                fh.write(struct.pack(tiff.tagformat2, count, packedvalue))
                if erase:
                    fh.seek(self.valueoffset)
                    fh.write(b'\x00' * oldsize)
            elif newsize <= oldsize or self.valueoffset + oldsize == fh.size:
                # separate -> separate smaller: overwrite, erase remaining
                fh.seek(self.offset + 4)
                fh.write(struct.pack(tiff.tagformat2[:2], count))
                fh.seek(self.valueoffset)
                fh.write(packedvalue)
                if erase and oldsize - newsize > 0:
                    fh.write(b'\x00' * (oldsize - newsize))
            else:
                # separate -> separate larger: erase old value, append to file
                fh.seek(0, os.SEEK_END)
                valueoffset = fh.tell()
                if valueoffset % 2:
                    # value offset must begin on a word boundary
                    fh.write(b'\x00')
                    valueoffset += 1
                # write offset
                fh.seek(self.offset + 4)
                fh.write(
                    struct.pack(
                        tiff.tagformat2,
                        count,
                        struct.pack(tiff.offsetformat, valueoffset),
                    )
                )
                # write value
                fh.seek(valueoffset)
                fh.write(packedvalue)
                if erase:
                    fh.seek(self.valueoffset)
                    fh.write(b'\x00' * oldsize)

        finally:
            fh.seek(pos)  # must restore file position

        return TiffTag(
            self.parent,
            self.offset,
            self.code,
            dtype,
            count,
            value,
            valueoffset,
        )

    def _fix_lsm_bitspersample(self) -> None:
        """Correct LSM bitspersample tag.

        Old LSM writers may use a separate region for two 16-bit values,
        although they fit into the tag value element of the tag.

        """
        if self.code != 258 or self.count != 2:
            return
        # TODO: test this case; need example file
        logger().warning(f'{self!r} correcting LSM bitspersample tag')
        value = struct.pack('<HH', *self.value)
        self.valueoffset = struct.unpack('<I', value)[0]
        self.parent.filehandle.seek(self.valueoffset)
        self.value = struct.unpack('<HH', self.parent.filehandle.read(4))

    def __repr__(self) -> str:
        name = '|'.join(TIFF.TAGS.getall(self.code, []))
        if name:
            name = ' ' + name
        return f'<tifffile.TiffTag {self.code}{name} @{self.offset}>'

    def __str__(self) -> str:
        return self._str()

    def _str(self, detail: int = 0, width: int = 79) -> str:
        """Return string containing information about TiffTag."""
        height = 1 if detail <= 0 else 8 * detail
        dtype = self.dtype_name
        if self.count > 1:
            dtype += f'[{self.count}]'
        name = '|'.join(TIFF.TAGS.getall(self.code, []))
        if name:
            name = f'{self.code} {name} @{self.offset}'
        else:
            name = f'{self.code} @{self.offset}'
        line = f'TiffTag {name} {dtype} @{self.valueoffset} '
        line = line[:width]
        try:
            value = self.value
        except TiffFileError:
            value = 'CORRUPTED'
        else:
            try:
                if self.count == 1:
                    value = enumstr(value)
                else:
                    value = pformat(tuple(enumstr(v) for v in value))
            except Exception:
                if not isinstance(value, (tuple, list)):
                    pass
                elif height == 1:
                    value = value[:256]
                elif len(value) > 2048:
                    value = (
                        value[:1024] + value[-1024:]  # type: ignore[operator]
                    )
                value = pformat(value, width=width, height=height)
        if detail <= 0:
            line += '= '
            line += value[:width]
            line = line[:width]
        else:
            line += '\n' + value
        return line


@final
class TiffTags:
    """Multidict-like interface to TiffTag instances in TiffPage.

    Differences to a regular dict:

    - values are instances of :py:class:`TiffTag`.
    - keys are :py:attr:`TiffTag.code` (int).
    - multiple values can be stored per key.
    - can be indexed by :py:attr:`TiffTag.name` (`str`), slower than by key.
    - `iter()` returns values instead of keys.
    - `values()` and `items()` contain all values sorted by offset.
    - `len()` returns number of all values.
    - `get()` takes optional index argument.
    - some functions are not implemented, such as, `update` and `pop`.

    """

    __slots__ = ('_dict', '_list')

    _dict: dict[int, TiffTag]
    _list: list[dict[int, TiffTag]]

    def __init__(self) -> None:
        self._dict = {}
        self._list = [self._dict]

    def add(self, tag: TiffTag, /) -> None:
        """Add tag."""
        code = tag.code
        for d in self._list:
            if code not in d:
                d[code] = tag
                break
        else:
            self._list.append({code: tag})

    def keys(self) -> list[int]:
        """Return codes of all tags."""
        return list(self._dict.keys())

    def values(self) -> list[TiffTag]:
        """Return all tags in order they are stored in file."""
        tags = (t for d in self._list for t in d.values())
        return sorted(tags, key=lambda t: t.offset)

    def items(self) -> list[tuple[int, TiffTag]]:
        """Return all (code, tag) pairs in order tags are stored in file."""
        items = (i for d in self._list for i in d.items())
        return sorted(items, key=lambda i: i[1].offset)

    def valueof(
        self,
        key: int | str,
        /,
        default: Any = None,
        index: int | None = None,
    ) -> Any:
        """Return value of tag by code or name if exists, else default.

        Parameters:
            key:
                Code or name of tag to return.
            default:
                Another value to return if specified tag is corrupted or
                not found.
            index:
                Specifies tag in case of multiple tags with identical code.
                The default is the first tag.

        """
        tag = self.get(key, default=None, index=index)
        if tag is None:
            return default
        try:
            return tag.value
        except TiffFileError:
            return default  # corrupted tag

    def get(
        self,
        key: int | str,
        /,
        default: TiffTag | None = None,
        index: int | None = None,
    ) -> TiffTag | None:
        """Return tag by code or name if exists, else default.

        Parameters:
            key:
                Code or name of tag to return.
            default:
                Another tag to return if specified tag is corrupted or
                not found.
            index:
                Specifies tag in case of multiple tags with identical code.
                The default is the first tag.

        """
        if index is None:
            if key in self._dict:
                return self._dict[cast(int, key)]
            if not isinstance(key, str):
                return default
            index = 0
        try:
            tags = self._list[index]
        except IndexError:
            return default
        if key in tags:
            return tags[cast(int, key)]
        if not isinstance(key, str):
            return default
        for tag in tags.values():
            if tag.name == key:
                return tag
        return default

    def getall(
        self,
        key: int | str,
        /,
        default: Any = None,
    ) -> list[TiffTag] | None:
        """Return list of all tags by code or name if exists, else default.

        Parameters:
            key:
                Code or name of tags to return.
            default:
                Value to return if no tags are found.

        """
        result: list[TiffTag] = []
        for tags in self._list:
            if key in tags:
                result.append(tags[cast(int, key)])
            else:
                break
        if result:
            return result
        if not isinstance(key, str):
            return default
        for tags in self._list:
            for tag in tags.values():
                if tag.name == key:
                    result.append(tag)
                    break
            if not result:
                break
        return result if result else default

    def __getitem__(self, key: int | str, /) -> TiffTag:
        """Return first tag by code or name. Raise KeyError if not found."""
        if key in self._dict:
            return self._dict[cast(int, key)]
        if not isinstance(key, str):
            raise KeyError(key)
        for tag in self._dict.values():
            if tag.name == key:
                return tag
        raise KeyError(key)

    def __setitem__(self, code: int, tag: TiffTag, /) -> None:
        """Add tag."""
        assert tag.code == code
        self.add(tag)

    def __delitem__(self, key: int | str, /) -> None:
        """Delete all tags by code or name."""
        found = False
        for tags in self._list:
            if key in tags:
                found = True
                del tags[cast(int, key)]
            else:
                break
        if found:
            return
        if not isinstance(key, str):
            raise KeyError(key)
        for tags in self._list:
            for tag in tags.values():
                if tag.name == key:
                    del tags[tag.code]
                    found = True
                    break
            else:
                break
        if not found:
            raise KeyError(key)
        return

    def __contains__(self, item: object, /) -> bool:
        """Return if tag is in map."""
        if item in self._dict:
            return True
        if not isinstance(item, str):
            return False
        for tag in self._dict.values():
            if tag.name == item:
                return True
        return False

    def __iter__(self) -> Iterator[TiffTag]:
        """Return iterator over all tags."""
        return iter(self.values())

    def __len__(self) -> int:
        """Return number of tags."""
        size = 0
        for d in self._list:
            size += len(d)
        return size

    def __repr__(self) -> str:
        return f'<tifffile.TiffTags @0x{id(self):016X}>'

    def __str__(self) -> str:
        return self._str()

    def _str(self, detail: int = 0, width: int = 79) -> str:
        """Return string with information about TiffTags."""
        info = []
        tlines = []
        vlines = []
        for tag in self:
            value = tag._str(width=width + 1)
            tlines.append(value[:width].strip())
            if detail > 0 and len(value) > width:
                try:
                    value = tag.value
                except Exception:
                    # delay load failed or closed file
                    continue
                if tag.code in {273, 279, 324, 325}:
                    if detail < 1:
                        value = value[:256]
                    elif len(value) > 1024:
                        value = value[:512] + value[-512:]
                    value = pformat(value, width=width, height=detail * 3)
                else:
                    value = pformat(value, width=width, height=detail * 8)
                if tag.count > 1:
                    vlines.append(
                        f'{tag.name} {tag.dtype_name}[{tag.count}]\n{value}'
                    )
                else:
                    vlines.append(f'{tag.name}\n{value}')
        info.append('\n'.join(tlines))
        if detail > 0 and vlines:
            info.append('\n')
            info.append('\n\n'.join(vlines))
        return '\n'.join(info)


@final
class TiffTagRegistry:
    """Registry of TIFF tag codes and names.

    Map tag codes and names to names and codes respectively.
    One tag code may be registered with several names, for example, 34853 is
    used for GPSTag or OlympusSIS2.
    Different tag codes may be registered with the same name, for example,
    37387 and 41483 are both named FlashEnergy.

    Parameters:
        arg: Mapping of codes to names.

    Examples:
        >>> tags = TiffTagRegistry([(34853, 'GPSTag'), (34853, 'OlympusSIS2')])
        >>> tags.add(37387, 'FlashEnergy')
        >>> tags.add(41483, 'FlashEnergy')
        >>> tags['GPSTag']
        34853
        >>> tags[34853]
        'GPSTag'
        >>> tags.getall(34853)
        ['GPSTag', 'OlympusSIS2']
        >>> tags.getall('FlashEnergy')
        [37387, 41483]
        >>> len(tags)
        4

    """

    __slots__ = ('_dict', '_list')

    _dict: dict[int | str, str | int]
    _list: list[dict[int | str, str | int]]

    def __init__(
        self,
        arg: TiffTagRegistry | dict[int, str] | Sequence[tuple[int, str]],
        /,
    ) -> None:
        self._dict = {}
        self._list = [self._dict]
        self.update(arg)

    def update(
        self,
        arg: TiffTagRegistry | dict[int, str] | Sequence[tuple[int, str]],
        /,
    ) -> None:
        """Add mapping of codes to names to registry.

        Parameters:
            arg: Mapping of codes to names.

        """
        if isinstance(arg, TiffTagRegistry):
            self._list.extend(arg._list)
            return
        if isinstance(arg, dict):
            arg = list(arg.items())
        for code, name in arg:
            self.add(code, name)

    def add(self, code: int, name: str, /) -> None:
        """Add code and name to registry."""
        for d in self._list:
            if code in d and d[code] == name:
                break
            if code not in d and name not in d:
                d[code] = name
                d[name] = code
                break
        else:
            self._list.append({code: name, name: code})

    def items(self) -> list[tuple[int, str]]:
        """Return all registry items as (code, name)."""
        items = (
            i for d in self._list for i in d.items() if isinstance(i[0], int)
        )
        return sorted(items, key=lambda i: i[0])  # type: ignore[arg-type]

    @overload
    def get(self, key: int, /, default: None) -> str | None: ...

    @overload
    def get(self, key: str, /, default: None) -> int | None: ...

    @overload
    def get(self, key: int, /, default: str) -> str: ...

    def get(
        self, key: int | str, /, default: str | None = None
    ) -> str | int | None:
        """Return first code or name if exists, else default.

        Parameters:
            key: tag code or name to lookup.
            default: value to return if key is not found.

        """
        for d in self._list:
            if key in d:
                return d[key]
        return default

    @overload
    def getall(self, key: int, /, default: None) -> list[str] | None: ...

    @overload
    def getall(self, key: str, /, default: None) -> list[int] | None: ...

    @overload
    def getall(self, key: int, /, default: list[str]) -> list[str]: ...

    def getall(
        self, key: int | str, /, default: list[str] | None = None
    ) -> list[str] | list[int] | None:
        """Return list of all codes or names if exists, else default.

        Parameters:
            key: tag code or name to lookup.
            default: value to return if key is not found.

        """
        result = [d[key] for d in self._list if key in d]
        return result if result else default  # type: ignore[return-value]

    @overload
    def __getitem__(self, key: int, /) -> str: ...

    @overload
    def __getitem__(self, key: str, /) -> int: ...

    def __getitem__(self, key: int | str, /) -> int | str:
        """Return first code or name. Raise KeyError if not found."""
        for d in self._list:
            if key in d:
                return d[key]
        raise KeyError(key)

    def __delitem__(self, key: int | str, /) -> None:
        """Delete all tags of code or name."""
        found = False
        for d in self._list:
            if key in d:
                found = True
                value = d[key]
                del d[key]
                del d[value]
        if not found:
            raise KeyError(key)

    def __contains__(self, item: int | str, /) -> bool:
        """Return if code or name is in registry."""
        for d in self._list:
            if item in d:
                return True
        return False

    def __iter__(self) -> Iterator[tuple[int, str]]:
        """Return iterator over all items in registry."""
        return iter(self.items())

    def __len__(self) -> int:
        """Return number of registered tags."""
        size = 0
        for d in self._list:
            size += len(d)
        return size // 2

    def __repr__(self) -> str:
        return f'<tifffile.TiffTagRegistry @0x{id(self):016X}>'

    def __str__(self) -> str:
        return 'TiffTagRegistry(((\n  {}\n))'.format(
            ',\n  '.join(f'({code}, {name!r})' for code, name in self.items())
        )


@final
class TiffPageSeries(Sequence[TiffPage | TiffFrame | None]):
    """Sequence of TIFF pages making up multi-dimensional image.

    Many TIFF based formats, such as OME-TIFF, use series of TIFF pages to
    store chunks of larger, multi-dimensional images.
    The image shape and position of chunks in the multi-dimensional image is
    defined in format-specific metadata.
    All pages in a series must have the same :py:meth:`TiffPage.hash`,
    that is, the same shape, data type, and storage properties.
    Items of a series may be None (missing) or instances of
    :py:class:`TiffPage` or :py:class:`TiffFrame`, possibly belonging to
    different files.

    Parameters:
        pages:
            List of TiffPage, TiffFrame, or None.
            The file handles of TiffPages or TiffFrames may not be open.
        shape:
            Shape of image array in series.
        dtype:
            Data type of image array in series.
        axes:
            Character codes for dimensions in shape.
            Length must match shape.
        attr:
            Arbitrary metadata associated with series.
        index:
            Index of series in multi-series files.
        parent:
            TiffFile instance series belongs to.
        name:
            Name of series.
        kind:
            Nature of series, such as, 'ome' or 'imagej'.
        truncated:
            Series is truncated, for example, ImageJ hyperstack > 4 GB.
        multifile:
            Series contains pages from multiple files.
        squeeze:
            Remove length-1 dimensions (except X and Y) from shape and axes
            by default.
        transform:
            Function to transform image data after decoding.

    """

    levels: list[TiffPageSeries]
    """Multi-resolution, pyramidal levels. ``levels[0] is self``."""

    parent: TiffFile | None
    """TiffFile instance series belongs to."""

    keyframe: TiffPage
    """TiffPage of series."""

    dtype: numpy.dtype[Any]
    """Data type (native byte order) of image array in series."""

    kind: str
    """Nature of series."""

    name: str
    """Name of image series from metadata."""

    transform: Callable[[NDArray[Any]], NDArray[Any]] | None
    """Function to transform image data after decoding."""

    is_multifile: bool
    """Series contains pages from multiple files."""

    is_truncated: bool
    """Series contains single page describing multi-dimensional image."""

    _pages: list[TiffPage | TiffFrame | None]
    # List of pages in series.
    # Might contain only first page of contiguous series

    _index: int  # index of series in multi-series files
    _squeeze: bool
    _axes: str
    _axes_squeezed: str
    _shape: tuple[int, ...]
    _shape_squeezed: tuple[int, ...]
    _len: int
    _attr: dict[str, Any]

    def __init__(
        self,
        pages: Sequence[TiffPage | TiffFrame | None],
        /,
        shape: Sequence[int] | None = None,
        dtype: DTypeLike | None = None,
        axes: str | None = None,
        *,
        attr: dict[str, Any] | None = None,
        coords: Mapping[str, NDArray[Any] | None] | None = None,
        index: int | None = None,
        parent: TiffFile | None = None,
        name: str | None = None,
        kind: str | None = None,
        truncated: bool = False,
        multifile: bool = False,
        squeeze: bool = True,
        transform: Callable[[NDArray[Any]], NDArray[Any]] | None = None,
    ) -> None:
        self._shape = ()
        self._shape_squeezed = ()
        self._axes = ''
        self._axes_squeezed = ''
        self._attr = {} if attr is None else dict(attr)

        self._index = int(index) if index else 0
        self._pages = list(pages)
        self.levels = [self]
        npages = len(self._pages)
        try:
            # find open TiffPage
            keyframe = next(
                p.keyframe
                for p in self._pages
                if p is not None
                and p.keyframe is not None
                and not p.keyframe.parent.filehandle.closed
            )
        except StopIteration:
            keyframe = next(
                p.keyframe
                for p in self._pages
                if p is not None and p.keyframe is not None
            )

        if shape is None:
            shape = keyframe.shape
        if axes is None:
            axes = keyframe.axes
        if dtype is None:
            dtype = keyframe.dtype

        self.dtype = numpy.dtype(dtype)
        self.kind = kind if kind else ''
        self.name = name if name else ''
        self.transform = transform
        self.keyframe = keyframe
        self.is_multifile = bool(multifile)
        self.is_truncated = bool(truncated)

        if parent is not None:
            self.parent = parent
        elif self._pages:
            self.parent = self.keyframe.parent
        else:
            self.parent = None

        self._set_dimensions(shape, axes, coords, squeeze)

        if not truncated and npages == 1:
            s = product(keyframe.shape)
            if s > 0:
                self._len = int(product(self.shape) // s)
            else:
                self._len = npages
        else:
            self._len = npages

    def _set_dimensions(
        self,
        shape: Sequence[int],
        axes: str,
        coords: Mapping[str, NDArray[Any] | None] | None = None,
        squeeze: bool = True,
        /,
    ) -> None:
        """Set shape, axes, and coords."""
        self._squeeze = bool(squeeze)
        self._shape = tuple(shape)
        self._axes = axes
        self._shape_squeezed, self._axes_squeezed, _ = squeeze_axes(
            shape, axes
        )

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of image array in series."""
        return self._shape_squeezed if self._squeeze else self._shape

    @property
    def axes(self) -> str:
        """Character codes for dimensions in image array."""
        return self._axes_squeezed if self._squeeze else self._axes

    @property
    def coords(self) -> dict[str, NDArray[Any]]:
        """Ordered map of dimension names to coordinate arrays."""
        raise NotImplementedError
        # return {
        #     name: numpy.arange(size)
        #     for name, size in zip(self.dims, self.shape)
        # }

    def get_shape(self, squeeze: bool | None = None) -> tuple[int, ...]:
        """Return default, squeezed, or expanded shape of series.

        Parameters:
            squeeze: Remove length-1 dimensions from shape.

        """
        if squeeze is None:
            squeeze = self._squeeze
        return self._shape_squeezed if squeeze else self._shape

    def get_axes(self, squeeze: bool | None = None) -> str:
        """Return default, squeezed, or expanded axes of series.

        Parameters:
            squeeze: Remove length-1 dimensions from axes.

        """
        if squeeze is None:
            squeeze = self._squeeze
        return self._axes_squeezed if squeeze else self._axes

    def get_coords(
        self, squeeze: bool | None = None
    ) -> dict[str, NDArray[Any]]:
        """Return default, squeezed, or expanded coords of series.

        Parameters:
            squeeze: Remove length-1 dimensions from coords.

        """
        raise NotImplementedError

    def asarray(
        self, *, level: int | None = None, **kwargs: Any
    ) -> NDArray[Any]:
        """Return images from series of pages as NumPy array.

        Parameters:
            level:
                Pyramid level to return.
                By default, the base layer is returned.
            **kwargs:
                Additional arguments passed to :py:meth:`TiffFile.asarray`.

        """
        if self.parent is None:
            raise ValueError('no parent')
        if level is not None:
            return self.levels[level].asarray(**kwargs)
        result = self.parent.asarray(series=self, **kwargs)
        if self.transform is not None:
            result = self.transform(result)
        return result

    def aszarr(
        self, *, level: int | None = None, **kwargs: Any
    ) -> ZarrTiffStore:
        """Return image array from series of pages as Zarr 2 store.

        Parameters:
            level:
                Pyramid level to return.
                By default, a multi-resolution store is returned.
            **kwargs:
                Additional arguments passed to :py:class:`ZarrTiffStore`.

        """
        if self.parent is None:
            raise ValueError('no parent')
        return ZarrTiffStore(self, level=level, **kwargs)

    @cached_property
    def dataoffset(self) -> int | None:
        """Offset to contiguous image data in file."""
        if not self._pages:
            return None
        pos = 0
        for page in self._pages:
            if page is None or len(page.dataoffsets) == 0:
                return None
            if not page.is_final:
                return None
            if not pos:
                pos = page.dataoffsets[0] + page.nbytes
                continue
            if pos != page.dataoffsets[0]:
                return None
            pos += page.nbytes

        page = self._pages[0]
        if page is None or len(page.dataoffsets) == 0:
            return None
        offset = page.dataoffsets[0]
        if (
            len(self._pages) == 1
            and isinstance(page, TiffPage)
            and (page.is_imagej or page.is_shaped or page.is_stk)
        ):
            # truncated files
            return offset
        if pos == offset + product(self.shape) * self.dtype.itemsize:
            return offset
        return None

    @property
    def is_pyramidal(self) -> bool:
        """Series contains multiple resolutions."""
        return len(self.levels) > 1

    @cached_property
    def attr(self) -> dict[str, Any]:
        """Arbitrary metadata associated with series."""
        return self._attr

    @property
    def ndim(self) -> int:
        """Number of array dimensions."""
        return len(self.shape)

    @property
    def dims(self) -> tuple[str, ...]:
        """Names of dimensions in image array."""
        # return tuple(self.coords.keys())
        return tuple(
            unique_strings(TIFF.AXES_NAMES.get(ax, ax) for ax in self.axes)
        )

    @property
    def sizes(self) -> dict[str, int]:
        """Ordered map of dimension names to lengths."""
        # return dict(zip(self.coords.keys(), self.shape))
        return dict(zip(self.dims, self.shape))

    @cached_property
    def size(self) -> int:
        """Number of elements in array."""
        return product(self.shape)

    @cached_property
    def nbytes(self) -> int:
        """Number of bytes in array."""
        return self.size * self.dtype.itemsize

    @property
    def pages(self) -> TiffPageSeries:
        # sequence of TiffPages or TiffFrame in series
        # a workaround to keep the old interface working
        return self

    def _getitem(self, key: int, /) -> TiffPage | TiffFrame | None:
        """Return specified page of series from cache or file."""
        key = int(key)
        if key < 0:
            key %= self._len
        if len(self._pages) == 1 and 0 < key < self._len:
            page = self._pages[0]
            assert page is not None
            assert self.parent is not None
            return self.parent.pages._getitem(page.index + key)
        return self._pages[key]

    @overload
    def __getitem__(
        self, key: int | numpy.integer[Any], /
    ) -> TiffPage | TiffFrame | None: ...

    @overload
    def __getitem__(
        self, key: slice | Iterable[int], /
    ) -> list[TiffPage | TiffFrame | None]: ...

    def __getitem__(
        self, key: int | numpy.integer[Any] | slice | Iterable[int], /
    ) -> TiffPage | TiffFrame | list[TiffPage | TiffFrame | None] | None:
        """Return specified page(s)."""
        if isinstance(key, (int, numpy.integer)):
            return self._getitem(int(key))
        if isinstance(key, slice):
            return [self._getitem(i) for i in range(*key.indices(self._len))]
        if isinstance(key, Iterable) and not isinstance(key, str):
            return [self._getitem(k) for k in key]
        raise TypeError('key must be an integer, slice, or iterable')

    def __iter__(self) -> Iterator[TiffPage | TiffFrame | None]:
        """Return iterator over pages in series."""
        if len(self._pages) == self._len:
            yield from self._pages
        else:
            assert self.parent is not None and self._pages[0] is not None
            pages = self.parent.pages
            index = self._pages[0].index
            for i in range(self._len):
                yield pages[index + i]

    def __len__(self) -> int:
        """Return number of pages in series."""
        return self._len

    def __repr__(self) -> str:
        return f'<tifffile.TiffPageSeries {self._index} {self.kind}>'

    def __str__(self) -> str:
        s = '  '.join(
            s
            for s in (
                snipstr(f'{self.name!r}', 20) if self.name else '',
                'x'.join(str(i) for i in self.shape),
                str(self.dtype),
                self.axes,
                self.kind,
                (f'{len(self.levels)} Levels') if self.is_pyramidal else '',
                f'{len(self)} Pages',
                (f'@{self.dataoffset}') if self.dataoffset else '',
            )
            if s
        )
        return f'TiffPageSeries {self._index}  {s}'


# TODO: this interface does not expose index keys except in __getitem__
class ZarrStore(MutableMapping[str, bytes]):
    """Zarr 2 store base class.

    ZarrStore instances must be closed with :py:meth:`ZarrStore.close`,
    which is automatically called when using the 'with' context manager.

    Parameters:
        fillvalue:
            Value to use for missing chunks of Zarr 2 store.
            The default is 0.
        chunkmode:
            Specifies how to chunk data.

    References:
        1. https://zarr.readthedocs.io/en/stable/spec/v2.html
        2. https://forum.image.sc/t/multiscale-arrays-v0-1/37930

    """

    _store: dict[str, Any]
    _fillvalue: int | float
    _chunkmode: int

    def __init__(
        self,
        /,
        *,
        fillvalue: int | float | None = None,
        chunkmode: CHUNKMODE | int | str | None = None,
    ) -> None:
        self._store = {}
        self._fillvalue = 0 if fillvalue is None else fillvalue
        if chunkmode is None:
            self._chunkmode = CHUNKMODE(0)
        else:
            self._chunkmode = enumarg(CHUNKMODE, chunkmode)

    def __enter__(self) -> ZarrStore:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def close(self) -> None:
        """Close ZarrStore."""

    def flush(self) -> None:
        """Flush ZarrStore."""
        raise PermissionError('ZarrStore is read-only')

    def clear(self) -> None:
        """Clear ZarrStore."""
        raise PermissionError('ZarrStore is read-only')

    def keys(self) -> KeysView[str]:
        """Return keys in ZarrStore."""
        return self._store.keys()

    def items(self) -> ItemsView[str, Any]:
        """Return items in ZarrStore."""
        return self._store.items()

    def values(self) -> ValuesView[Any]:
        """Return values in ZarrStore."""
        return self._store.values()

    def __iter__(self) -> Iterator[str]:
        return iter(self._store)

    def __len__(self) -> int:
        return len(self._store)

    def __contains__(self, key: object, /) -> bool:
        if key in self._store:
            return True
        assert isinstance(key, str)
        return self._contains(key)

    def _contains(self, key: str, /) -> bool:
        """Return if key is in store."""
        raise NotImplementedError

    def __delitem__(self, key: object, /) -> None:
        raise PermissionError('ZarrStore is read-only')

    def __getitem__(self, key: str, /) -> Any:
        if key in self._store:
            return self._store[key]
        if key[-7:] == '.zarray' or key[-7:] == '.zgroup':
            # catch '.zarray' and 'attribute/.zarray'
            raise KeyError(key)
        return self._getitem(key)

    def _getitem(self, key: str, /) -> NDArray[Any]:
        """Return chunk from file."""
        raise NotImplementedError

    def __setitem__(self, key: str, value: bytes, /) -> None:
        if key in self._store:
            raise KeyError(key)
        if key[-7:] == '.zarray' or key[-7:] == '.zgroup':
            # catch '.zarray' and 'attribute/.zarray'
            raise KeyError(key)
        return self._setitem(key, value)

    def _setitem(self, key: str, value: bytes, /) -> None:
        """Write chunk to file."""
        raise NotImplementedError

    @property
    def is_multiscales(self) -> bool:
        """Return if ZarrStore is multi-scales."""
        return b'multiscales' in self._store['.zattrs']

    @staticmethod
    def _empty_chunk(
        shape: tuple[int, ...],
        dtype: DTypeLike,
        fillvalue: int | float | None,
        /,
    ) -> NDArray[Any]:
        """Return empty chunk."""
        if fillvalue is None or fillvalue == 0:
            # return bytes(product(shape) * dtype.itemsize)
            return numpy.zeros(shape, dtype)
        chunk = numpy.empty(shape, dtype)
        chunk[:] = fillvalue
        return chunk  # .tobytes()

    @staticmethod
    def _dtype_str(dtype: numpy.dtype[Any], /) -> str:
        """Return dtype as string with native byte order."""
        if dtype.itemsize == 1:
            byteorder = '|'
        else:
            byteorder = {'big': '>', 'little': '<'}[sys.byteorder]
        return byteorder + dtype.str[1:]

    @staticmethod
    def _json(obj: Any, /) -> bytes:
        """Serialize object to JSON formatted string."""
        return json.dumps(
            obj,
            indent=1,
            sort_keys=True,
            ensure_ascii=True,
            separators=(',', ': '),
        ).encode('ascii')

    @staticmethod
    def _value(value: Any, dtype: numpy.dtype[Any], /) -> Any:
        """Return value which is serializable to JSON."""
        if value is None:
            return value
        if dtype.kind == 'b':
            return bool(value)
        if dtype.kind in 'ui':
            return int(value)
        if dtype.kind == 'f':
            if numpy.isnan(value):
                return 'NaN'
            if numpy.isposinf(value):
                return 'Infinity'
            if numpy.isneginf(value):
                return '-Infinity'
            return float(value)
        if dtype.kind == 'c':
            value = numpy.array(value, dtype)
            return (
                ZarrStore._value(value.real, dtype.type().real.dtype),
                ZarrStore._value(value.imag, dtype.type().imag.dtype),
            )
        return value

    @staticmethod
    def _ndindex(
        shape: tuple[int, ...], chunks: tuple[int, ...], /
    ) -> Iterator[str]:
        """Return iterator over all chunk index strings."""
        assert len(shape) == len(chunks)
        chunked = tuple(
            i // j + (1 if i % j else 0) for i, j in zip(shape, chunks)
        )
        for indices in numpy.ndindex(chunked):
            yield '.'.join(str(index) for index in indices)


@final
class ZarrTiffStore(ZarrStore):
    """Zarr 2 store interface to image array in TiffPage or TiffPageSeries.

    ZarrTiffStore is using a TiffFile instance for reading and decoding chunks.
    Therefore, ZarrTiffStore instances cannot be pickled.

    For writing, image data must be stored in uncompressed, unpredicted,
    and unpacked form. Sparse strips and tiles are not written.

    Parameters:
        arg:
            TIFF page or series to wrap as Zarr 2 store.
        level:
            Pyramidal level to wrap. The default is 0.
        chunkmode:
            Use strips or tiles (0) or whole page data (2) as chunks.
            The default is 0.
        fillvalue:
            Value to use for missing chunks. The default is 0.
        zattrs:
            Additional attributes to store in `.zattrs`.
        multiscales:
            Create a multiscales compatible Zarr 2 group store.
            By default, create a Zarr 2 array store for pages and non-pyramidal
            series.
        lock:
            Reentrant lock to synchronize seeks and reads from file.
            By default, the lock of the parent's file handle is used.
        squeeze:
            Remove length-1 dimensions from shape of TiffPageSeries.
        maxworkers:
            Maximum number of threads to concurrently decode strips or tiles
            if `chunkmode=2`.
            If *None* or *0*, use up to :py:attr:`_TIFF.MAXWORKERS` threads.
        buffersize:
            Approximate number of bytes to read from file in one pass
            if `chunkmode=2`. The default is :py:attr:`_TIFF.BUFFERSIZE`.
        _openfiles:
            Internal API.

    """

    _data: list[TiffPageSeries]
    _filecache: FileCache
    _transform: Callable[[NDArray[Any]], NDArray[Any]] | None
    _maxworkers: int | None
    _buffersize: int | None
    _squeeze: bool | None
    _writable: bool
    _multiscales: bool

    def __init__(
        self,
        arg: TiffPage | TiffFrame | TiffPageSeries,
        /,
        *,
        level: int | None = None,
        chunkmode: CHUNKMODE | int | str | None = None,
        fillvalue: int | float | None = None,
        zattrs: dict[str, Any] | None = None,
        multiscales: bool | None = None,
        lock: threading.RLock | NullContext | None = None,
        squeeze: bool | None = None,
        maxworkers: int | None = None,
        buffersize: int | None = None,
        _openfiles: int | None = None,
    ) -> None:
        super().__init__(fillvalue=fillvalue, chunkmode=chunkmode)

        if self._chunkmode not in {0, 2}:
            raise NotImplementedError(f'{self._chunkmode!r} not implemented')

        self._squeeze = None if squeeze is None else bool(squeeze)
        self._maxworkers = maxworkers
        self._buffersize = buffersize

        if isinstance(arg, TiffPageSeries):
            self._data = arg.levels
            self._transform = arg.transform
            if multiscales is not None and not multiscales:
                level = 0
            if level is not None:
                self._data = [self._data[level]]
            name = arg.name
        else:
            self._data = [TiffPageSeries([arg])]
            self._transform = None
            name = 'Unnamed'

        fh = self._data[0].keyframe.parent._parent.filehandle
        self._writable = fh.writable() and self._chunkmode == 0
        if lock is None:
            fh.set_lock(True)
            lock = fh.lock
        self._filecache = FileCache(size=_openfiles, lock=lock)

        zattrs = {} if zattrs is None else dict(zattrs)
        # TODO: Zarr Encoding Specification
        # https://xarray.pydata.org/en/stable/internals/zarr-encoding-spec.html

        if multiscales or len(self._data) > 1:
            # multiscales
            self._multiscales = True
            if '_ARRAY_DIMENSIONS' in zattrs:
                array_dimensions = zattrs.pop('_ARRAY_DIMENSIONS')
            else:
                array_dimensions = list(self._data[0].get_axes(squeeze))
            self._store['.zgroup'] = ZarrStore._json({'zarr_format': 2})
            self._store['.zattrs'] = ZarrStore._json(
                {
                    # TODO: use https://ngff.openmicroscopy.org/latest/
                    'multiscales': [
                        {
                            'version': '0.1',
                            'name': name,
                            'datasets': [
                                {'path': str(i)}
                                for i in range(len(self._data))
                            ],
                            # 'axes': [...]
                            # 'type': 'unknown',
                            'metadata': {},
                        }
                    ],
                    **zattrs,
                }
            )
            shape0 = self._data[0].get_shape(squeeze)
            for level, series in enumerate(self._data):
                keyframe = series.keyframe
                keyframe.decode  # cache decode function
                shape = series.get_shape(squeeze)
                dtype = series.dtype
                if fillvalue is None:
                    self._fillvalue = fillvalue = keyframe.nodata
                if self._chunkmode:
                    chunks = keyframe.shape
                else:
                    chunks = keyframe.chunks
                self._store[f'{level}/.zattrs'] = ZarrStore._json(
                    {
                        '_ARRAY_DIMENSIONS': [
                            (f'{ax}{level}' if i != j else ax)
                            for ax, i, j in zip(
                                array_dimensions, shape, shape0
                            )
                        ]
                    }
                )
                self._store[f'{level}/.zarray'] = ZarrStore._json(
                    {
                        'zarr_format': 2,
                        'shape': shape,
                        'chunks': ZarrTiffStore._chunks(chunks, shape),
                        'dtype': ZarrStore._dtype_str(dtype),
                        'compressor': None,
                        'fill_value': ZarrStore._value(fillvalue, dtype),
                        'order': 'C',
                        'filters': None,
                    }
                )
                if self._writable:
                    self._writable = ZarrTiffStore._is_writable(keyframe)
        else:
            self._multiscales = False
            series = self._data[0]
            keyframe = series.keyframe
            keyframe.decode  # cache decode function
            shape = series.get_shape(squeeze)
            dtype = series.dtype
            if fillvalue is None:
                self._fillvalue = fillvalue = keyframe.nodata
            if self._chunkmode:
                chunks = keyframe.shape
            else:
                chunks = keyframe.chunks
            if '_ARRAY_DIMENSIONS' not in zattrs:
                zattrs['_ARRAY_DIMENSIONS'] = list(series.get_axes(squeeze))
            self._store['.zattrs'] = ZarrStore._json(zattrs)
            self._store['.zarray'] = ZarrStore._json(
                {
                    'zarr_format': 2,
                    'shape': shape,
                    'chunks': ZarrTiffStore._chunks(chunks, shape),
                    'dtype': ZarrStore._dtype_str(dtype),
                    'compressor': None,
                    'fill_value': ZarrStore._value(fillvalue, dtype),
                    'order': 'C',
                    'filters': None,
                }
            )
            if self._writable:
                self._writable = ZarrTiffStore._is_writable(keyframe)

    def close(self) -> None:
        """Close open file handles."""
        if hasattr(self, '_filecache'):
            self._filecache.clear()

    def write_fsspec(
        self,
        jsonfile: str | os.PathLike[Any] | TextIO,
        /,
        url: str,
        *,
        groupname: str | None = None,
        templatename: str | None = None,
        compressors: dict[COMPRESSION | int, str | None] | None = None,
        version: int | None = None,
        _shape: Sequence[int] | None = None,
        _axes: Sequence[str] | None = None,
        _index: Sequence[int] | None = None,
        _append: bool = False,
        _close: bool = True,
    ) -> None:
        """Write fsspec ReferenceFileSystem as JSON to file.

        Parameters:
            jsonfile:
                Name or open file handle of output JSON file.
            url:
                Remote location of TIFF file(s) without file name(s).
            groupname:
                Zarr 2 group name.
            templatename:
                Version 1 URL template name. The default is 'u'.
            compressors:
                Mapping of :py:class:`COMPRESSION` codes to Numcodecs codec
                names.
            version:
                Version of fsspec file to write. The default is 0.
            _shape:
                Shape of file sequence (experimental).
            _axes:
                Axes of file sequence (experimental).
            _index
                Index of file in sequence (experimental).
            _append:
                If *True*, only write index keys and values (experimental).
            _close:
                If *True*, no more appends (experimental).

        Raises:
            ValueError:
                ZarrTiffStore cannot be represented as ReferenceFileSystem
                due to features that are not supported by Zarr 2, Numcodecs,
                or Imagecodecs:

                - compressors, such as CCITT
                - filters, such as bitorder reversal, packed integers
                - dtypes, such as float24, complex integers
                - JPEGTables in multi-page series
                - incomplete chunks, such as `imagelength % rowsperstrip != 0`

                Files containing incomplete tiles may fail at runtime.

        Notes:
            Parameters `_shape`,  `_axes`, `_index`, `_append`, and `_close`
            are an experimental API for joining the ReferenceFileSystems of
            multiple files of a TiffSequence.

        References:
            - `fsspec ReferenceFileSystem format
              <https://github.com/fsspec/kerchunk>`_

        """
        compressors = {
            1: None,
            8: 'zlib',
            32946: 'zlib',
            34925: 'lzma',
            50013: 'zlib',  # pixtiff
            5: 'imagecodecs_lzw',
            7: 'imagecodecs_jpeg',
            22610: 'imagecodecs_jpegxr',
            32773: 'imagecodecs_packbits',
            33003: 'imagecodecs_jpeg2k',
            33004: 'imagecodecs_jpeg2k',
            33005: 'imagecodecs_jpeg2k',
            33007: 'imagecodecs_jpeg',
            34712: 'imagecodecs_jpeg2k',
            34887: 'imagecodecs_lerc',
            34892: 'imagecodecs_jpeg',
            34933: 'imagecodecs_png',
            34934: 'imagecodecs_jpegxr',
            48124: 'imagecodecs_jetraw',
            50000: 'imagecodecs_zstd',  # numcodecs.zstd fails w/ unknown sizes
            50001: 'imagecodecs_webp',
            50002: 'imagecodecs_jpegxl',
            52546: 'imagecodecs_jpegxl',
            **({} if compressors is None else compressors),
        }

        for series in self._data:
            errormsg = ' not supported by the fsspec ReferenceFileSystem'
            keyframe = series.keyframe
            if (
                keyframe.compression in {65000, 65001, 65002}
                and keyframe.parent.is_eer
            ):
                compressors[keyframe.compression] = 'imagecodecs_eer'
            if keyframe.compression not in compressors:
                raise ValueError(f'{keyframe.compression!r} is' + errormsg)
            if keyframe.fillorder != 1:
                raise ValueError(f'{keyframe.fillorder!r} is' + errormsg)
            if keyframe.sampleformat not in {1, 2, 3, 6}:
                # TODO: support float24 and cint via filters?
                raise ValueError(f'{keyframe.sampleformat!r} is' + errormsg)
            if (
                keyframe.bitspersample
                not in {
                    8,
                    16,
                    32,
                    64,
                    128,
                }
                and keyframe.compression
                not in {
                    # JPEG
                    7,
                    33007,
                    34892,
                }
                and compressors[keyframe.compression] != 'imagecodecs_eer'
            ):
                raise ValueError(
                    f'BitsPerSample {keyframe.bitspersample} is' + errormsg
                )
            if (
                not self._chunkmode
                and not keyframe.is_tiled
                and keyframe.imagelength % keyframe.rowsperstrip
            ):
                raise ValueError('incomplete chunks are' + errormsg)
            if self._chunkmode and not keyframe.is_final:
                raise ValueError(f'{self._chunkmode!r} is' + errormsg)
            if keyframe.jpegtables is not None and len(series.pages) > 1:
                raise ValueError(
                    'JPEGTables in multi-page files are' + errormsg
                )

        if url is None:
            url = ''
        elif url and url[-1] != '/':
            url += '/'
        url = url.replace('\\', '/')

        if groupname is None:
            groupname = ''
        elif groupname and groupname[-1] != '/':
            groupname += '/'

        byteorder: ByteOrder | None = '<' if sys.byteorder == 'big' else '>'
        if (
            self._data[0].keyframe.parent.byteorder != byteorder
            or self._data[0].keyframe.dtype is None
            or self._data[0].keyframe.dtype.itemsize == 1
        ):
            byteorder = None

        index: str
        _shape = [] if _shape is None else list(_shape)
        _axes = [] if _axes is None else list(_axes)
        if len(_shape) != len(_axes):
            raise ValueError('len(_shape) != len(_axes)')
        if _index is None:
            index = ''
        elif len(_shape) != len(_index):
            raise ValueError('len(_shape) != len(_index)')
        elif _index:
            index = '.'.join(str(i) for i in _index)
            index += '.'

        refs: dict[str, Any] = {}
        refzarr: dict[str, Any]
        if version == 1:
            if _append:
                raise ValueError('cannot append to version 1')
            if templatename is None:
                templatename = 'u'
            refs['version'] = 1
            refs['templates'] = {}
            refs['gen'] = []
            templates = {}
            if self._data[0].is_multifile:
                i = 0
                for page in self._data[0].pages:
                    if page is None or page.keyframe is None:
                        continue
                    fname = page.keyframe.parent.filehandle.name
                    if fname in templates:
                        continue
                    key = f'{templatename}{i}'
                    templates[fname] = f'{{{{{key}}}}}'
                    refs['templates'][key] = url + fname
                    i += 1
            else:
                fname = self._data[0].keyframe.parent.filehandle.name
                key = f'{templatename}'
                templates[fname] = f'{{{{{key}}}}}'
                refs['templates'][key] = url + fname

            refs['refs'] = refzarr = {}
        else:
            refzarr = refs

        if not _append:
            if groupname:
                # TODO: support nested groups
                refzarr['.zgroup'] = ZarrStore._json(
                    {'zarr_format': 2}
                ).decode()

            for key, value in self._store.items():
                if '.zattrs' in key and _axes:
                    value = json.loads(value)
                    if '_ARRAY_DIMENSIONS' in value:
                        value['_ARRAY_DIMENSIONS'] = (
                            _axes + value['_ARRAY_DIMENSIONS']
                        )
                    value = ZarrStore._json(value)
                elif '.zarray' in key:
                    level = int(key.split('/')[0]) if '/' in key else 0
                    keyframe = self._data[level].keyframe
                    value = json.loads(value)
                    if _shape:
                        value['shape'] = _shape + value['shape']
                        value['chunks'] = [1] * len(_shape) + value['chunks']
                    codec_id = compressors[keyframe.compression]
                    if codec_id == 'imagecodecs_jpeg':
                        # TODO: handle JPEG color spaces
                        jpegtables = keyframe.jpegtables
                        if jpegtables is None:
                            tables = None
                        else:
                            import base64

                            tables = base64.b64encode(jpegtables).decode()
                        jpegheader = keyframe.jpegheader
                        if jpegheader is None:
                            header = None
                        else:
                            import base64

                            header = base64.b64encode(jpegheader).decode()
                        (
                            colorspace_jpeg,
                            colorspace_data,
                        ) = jpeg_decode_colorspace(
                            keyframe.photometric,
                            keyframe.planarconfig,
                            keyframe.extrasamples,
                            keyframe.is_jfif,
                        )
                        value['compressor'] = {
                            'id': codec_id,
                            'tables': tables,
                            'header': header,
                            'bitspersample': keyframe.bitspersample,
                            'colorspace_jpeg': colorspace_jpeg,
                            'colorspace_data': colorspace_data,
                        }
                    elif (
                        codec_id == 'imagecodecs_webp'
                        and keyframe.samplesperpixel == 4
                    ):
                        value['compressor'] = {
                            'id': codec_id,
                            'hasalpha': True,
                        }
                    elif codec_id == 'imagecodecs_eer':
                        if keyframe.compression == 65002:
                            rlebits = int(keyframe.tags.valueof(65007, 7))
                            horzbits = int(keyframe.tags.valueof(65008, 2))
                            vertbits = int(keyframe.tags.valueof(65009, 2))
                        elif keyframe.compression == 65001:
                            rlebits = 7
                            horzbits = 2
                            vertbits = 2
                        else:
                            rlebits = 8
                            horzbits = 2
                            vertbits = 2
                        value['compressor'] = {
                            'id': codec_id,
                            'shape': keyframe.chunks,
                            'rlebits': rlebits,
                            'horzbits': horzbits,
                            'vertbits': vertbits,
                        }
                    elif codec_id is not None:
                        value['compressor'] = {'id': codec_id}
                    if byteorder is not None:
                        value['dtype'] = byteorder + value['dtype'][1:]
                    if keyframe.predictor > 1:
                        # predictors need access to chunk shape and dtype
                        # requires imagecodecs > 2021.8.26 to read
                        if keyframe.predictor in {2, 34892, 34893}:
                            filter_id = 'imagecodecs_delta'
                        else:
                            filter_id = 'imagecodecs_floatpred'
                        if keyframe.predictor <= 3:
                            dist = 1
                        elif keyframe.predictor in {34892, 34894}:
                            dist = 2
                        else:
                            dist = 4
                        if (
                            keyframe.planarconfig == 1
                            and keyframe.samplesperpixel > 1
                        ):
                            axis = -2
                        else:
                            axis = -1
                        value['filters'] = [
                            {
                                'id': filter_id,
                                'axis': axis,
                                'dist': dist,
                                'shape': value['chunks'],
                                'dtype': value['dtype'],
                            }
                        ]
                    value = ZarrStore._json(value)

                refzarr[groupname + key] = value.decode()

        fh: TextIO
        if hasattr(jsonfile, 'write'):
            fh = jsonfile  # type: ignore[assignment]
        else:
            fh = open(jsonfile, 'w', encoding='utf-8')

        if version == 1:
            fh.write(json.dumps(refs, indent=1).rsplit('}"', 1)[0] + '}"')
            indent = '  '
        elif _append:
            indent = ' '
        else:
            fh.write(json.dumps(refs, indent=1)[:-2])
            indent = ' '

        for key, value in self._store.items():
            if '.zarray' in key:
                value = json.loads(value)
                shape = value['shape']
                chunks = value['chunks']
                levelstr = (key.split('/')[0] + '/') if '/' in key else ''
                for chunkindex in ZarrStore._ndindex(shape, chunks):
                    key = levelstr + chunkindex
                    keyframe, page, _, offset, bytecount = self._parse_key(key)
                    key = levelstr + index + chunkindex
                    if page and self._chunkmode and offset is None:
                        offset = page.dataoffsets[0]
                        bytecount = keyframe.nbytes
                    if offset and bytecount:
                        fname = keyframe.parent.filehandle.name
                        if version == 1:
                            fname = templates[fname]
                        else:
                            fname = f'{url}{fname}'
                        fh.write(
                            f',\n{indent}"{groupname}{key}": '
                            f'["{fname}", {offset}, {bytecount}]'
                        )

        # TODO: support nested groups
        if version == 1:
            fh.write('\n }\n}')
        elif _close:
            fh.write('\n}')

        if not hasattr(jsonfile, 'write'):
            fh.close()

    def _contains(self, key: str, /) -> bool:
        """Return if key is in store."""
        try:
            _, page, _, offset, bytecount = self._parse_key(key)
        except (KeyError, IndexError):
            return False
        if self._chunkmode and offset is None:
            return True
        return (
            page is not None
            and offset is not None
            and bytecount is not None
            and offset > 0
            and bytecount > 0
        )

    def _getitem(self, key: str, /) -> NDArray[Any]:
        """Return chunk from file."""
        keyframe, page, chunkindex, offset, bytecount = self._parse_key(key)

        if page is None or offset == 0 or bytecount == 0:
            raise KeyError(key)

        fh = page.parent.filehandle

        if self._chunkmode:
            if offset is not None:
                # contiguous image data in page or series
                # create virtual frame instead of loading page from file
                assert bytecount is not None
                page = TiffFrame(
                    page.parent,
                    index=0,
                    keyframe=keyframe,
                    dataoffsets=(offset,),
                    databytecounts=(bytecount,),
                )
            self._filecache.open(fh)
            chunk = page.asarray(
                lock=self._filecache.lock,
                maxworkers=self._maxworkers,
                buffersize=self._buffersize,
            )
            self._filecache.close(fh)
            if self._transform is not None:
                chunk = self._transform(chunk)
            return chunk

        assert offset is not None and bytecount is not None
        chunk_bytes = self._filecache.read(fh, offset, bytecount)

        decodeargs: dict[str, Any] = {'_fullsize': True}
        if page.jpegtables is not None:
            decodeargs['jpegtables'] = page.jpegtables
        if keyframe.jpegheader is not None:
            decodeargs['jpegheader'] = keyframe.jpegheader

        assert chunkindex is not None
        chunk = keyframe.decode(
            chunk_bytes, chunkindex, **decodeargs  # type: ignore[assignment]
        )[0]
        assert chunk is not None
        if self._transform is not None:
            chunk = self._transform(chunk)

        if self._chunkmode:
            chunks = keyframe.shape
        else:
            chunks = keyframe.chunks
        if chunk.size != product(chunks):
            raise RuntimeError(f'{chunk.size} != {product(chunks)}')
        return chunk  # .tobytes()

    def _setitem(self, key: str, value: bytes, /) -> None:
        """Write chunk to file."""
        if not self._writable:
            raise PermissionError('ZarrStore is read-only')
        keyframe, page, chunkindex, offset, bytecount = self._parse_key(key)
        if (
            page is None
            or offset is None
            or offset == 0
            or bytecount is None
            or bytecount == 0
        ):
            return
        if bytecount < len(value):
            value = value[:bytecount]
        self._filecache.write(page.parent.filehandle, offset, value)

    def _parse_key(self, key: str, /) -> tuple[
        TiffPage,
        TiffPage | TiffFrame | None,
        int | None,
        int | None,
        int | None,
    ]:
        """Return keyframe, page, index, offset, and bytecount from key.

        Raise KeyError if key is not valid.

        """
        if self._multiscales:
            try:
                level, key = key.split('/')
                series = self._data[int(level)]
            except (ValueError, IndexError) as exc:
                raise KeyError(key) from exc
        else:
            series = self._data[0]
        keyframe = series.keyframe
        pageindex, chunkindex = self._indices(key, series)
        if series.dataoffset is not None:
            # contiguous or truncated
            page = series[0]
            if page is None or page.dtype is None or page.keyframe is None:
                return keyframe, None, chunkindex, 0, 0
            offset = pageindex * page.size * page.dtype.itemsize
            try:
                offset += page.dataoffsets[chunkindex]
            except IndexError as exc:
                raise KeyError(key) from exc
            if self._chunkmode:
                bytecount = page.size * page.dtype.itemsize
                return page.keyframe, page, chunkindex, offset, bytecount
        elif self._chunkmode:
            with self._filecache.lock:
                page = series[pageindex]
            if page is None or page.keyframe is None:
                return keyframe, None, None, 0, 0
            return page.keyframe, page, None, None, None
        else:
            with self._filecache.lock:
                page = series[pageindex]
            if page is None or page.keyframe is None:
                return keyframe, None, chunkindex, 0, 0
            try:
                offset = page.dataoffsets[chunkindex]
            except IndexError:
                # raise KeyError(key) from exc
                # issue #249: Philips may be missing last row of tiles
                return page.keyframe, page, chunkindex, 0, 0
        try:
            bytecount = page.databytecounts[chunkindex]
        except IndexError as exc:
            raise KeyError(key) from exc
        return page.keyframe, page, chunkindex, offset, bytecount

    def _indices(self, key: str, series: TiffPageSeries, /) -> tuple[int, int]:
        """Return page and strile indices from Zarr chunk index."""
        keyframe = series.keyframe
        shape = series.get_shape(self._squeeze)
        try:
            indices = [int(i) for i in key.split('.')]
        except ValueError as exc:
            raise KeyError(key) from exc
        assert len(indices) == len(shape)
        if self._chunkmode:
            chunked = (1,) * len(keyframe.shape)
        else:
            chunked = keyframe.chunked
        p = 1
        for i, s in enumerate(shape[::-1]):
            p *= s
            if p == keyframe.size:
                i = len(indices) - i - 1
                frames_indices = indices[:i]
                strile_indices = indices[i:]
                frames_chunked = shape[:i]
                strile_chunked = list(shape[i:])  # updated later
                break
        else:
            raise RuntimeError
        if len(strile_chunked) == len(keyframe.shape):
            strile_chunked = list(chunked)
        else:
            # get strile_chunked including singleton dimensions
            i = len(strile_indices) - 1
            j = len(keyframe.shape) - 1
            while True:
                if strile_chunked[i] == keyframe.shape[j]:
                    strile_chunked[i] = chunked[j]
                    i -= 1
                    j -= 1
                elif strile_chunked[i] == 1:
                    i -= 1
                else:
                    raise RuntimeError('shape does not match page shape')
                if i < 0 or j < 0:
                    break
            assert product(strile_chunked) == product(chunked)
        if len(frames_indices) > 0:
            frameindex = int(
                numpy.ravel_multi_index(frames_indices, frames_chunked)
            )
        else:
            frameindex = 0
        if len(strile_indices) > 0:
            strileindex = int(
                numpy.ravel_multi_index(strile_indices, strile_chunked)
            )
        else:
            strileindex = 0
        return frameindex, strileindex

    @staticmethod
    def _chunks(
        chunks: tuple[int, ...], shape: tuple[int, ...], /
    ) -> tuple[int, ...]:
        """Return chunks with same length as shape."""
        ndim = len(shape)
        if ndim == 0:
            return ()  # empty array
        if 0 in shape:
            return (1,) * ndim
        newchunks = []
        i = ndim - 1
        j = len(chunks) - 1
        while True:
            if j < 0:
                newchunks.append(1)
                i -= 1
            elif shape[i] > 1 and chunks[j] > 1:
                newchunks.append(chunks[j])
                i -= 1
                j -= 1
            elif shape[i] == chunks[j]:  # both 1
                newchunks.append(1)
                i -= 1
                j -= 1
            elif shape[i] == 1:
                newchunks.append(1)
                i -= 1
            elif chunks[j] == 1:
                newchunks.append(1)
                j -= 1
            else:
                raise RuntimeError
            if i < 0 or ndim == len(newchunks):
                break
        # assert ndim == len(newchunks)
        return tuple(newchunks[::-1])

    @staticmethod
    def _is_writable(keyframe: TiffPage) -> bool:
        """Return True if chunks are writable."""
        return (
            keyframe.compression == 1
            and keyframe.fillorder == 1
            and keyframe.sampleformat in {1, 2, 3, 6}
            and keyframe.bitspersample in {8, 16, 32, 64, 128}
            # and (
            #     keyframe.rowsperstrip == 0
            #     or keyframe.imagelength % keyframe.rowsperstrip == 0
            # )
        )

    def __enter__(self) -> ZarrTiffStore:
        return self

    def __repr__(self) -> str:
        return f'<tifffile.ZarrTiffStore @0x{id(self):016X}>'


@final
class ZarrFileSequenceStore(ZarrStore):
    """Zarr 2 store interface to image array in FileSequence.

    Parameters:
        filesequence:
            FileSequence instance to wrap as Zarr 2 store.
            Files in containers are not supported.
        fillvalue:
            Value to use for missing chunks. The default is 0.
        chunkmode:
            Currently only one chunk per file is supported.
        chunkshape:
            Shape of chunk in each file.
            Must match ``FileSequence.imread(file, **imreadargs).shape``.
        chunkdtype:
            Data type of chunk in each file.
            Must match ``FileSequence.imread(file, **imreadargs).dtype``.
        axestiled:
            Axes to be tiled. Map stacked sequence axis to chunk axis.
        zattrs:
            Additional attributes to store in `.zattrs`.
        imreadargs:
            Arguments passed to :py:attr:`FileSequence.imread`.
        **kwargs:
            Arguments passed to :py:attr:`FileSequence.imread`in addition
            to `imreadargs`.

    Notes:
        If `chunkshape` or `chunkdtype` are *None* (default), their values
        are determined by reading the first file with
        ``FileSequence.imread(arg.files[0], **imreadargs)``.

    """

    imread: Callable[..., NDArray[Any]]
    """Function to read image array from single file."""

    _lookup: dict[tuple[int, ...], str]
    _chunks: tuple[int, ...]
    _dtype: numpy.dtype[Any]
    _tiled: TiledSequence
    _commonpath: str
    _kwargs: dict[str, Any]

    def __init__(
        self,
        filesequence: FileSequence,
        /,
        *,
        fillvalue: int | float | None = None,
        chunkmode: CHUNKMODE | int | str | None = None,
        chunkshape: Sequence[int] | None = None,
        chunkdtype: DTypeLike | None = None,
        axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
        zattrs: dict[str, Any] | None = None,
        imreadargs: dict[str, Any] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(fillvalue=fillvalue, chunkmode=chunkmode)

        if self._chunkmode not in {0, 3}:
            raise ValueError(f'invalid chunkmode {self._chunkmode!r}')

        if not isinstance(filesequence, FileSequence):
            raise TypeError('not a FileSequence')

        if filesequence._container:
            raise NotImplementedError('cannot open container as Zarr 2 store')

        # TODO: deprecate kwargs?
        if imreadargs is not None:
            kwargs |= imreadargs

        self._kwargs = kwargs
        self._imread = filesequence.imread
        self._commonpath = filesequence.commonpath()

        if chunkshape is None or chunkdtype is None:
            chunk = filesequence.imread(filesequence[0], **kwargs)
            self._chunks = chunk.shape
            self._dtype = chunk.dtype
        else:
            self._chunks = tuple(chunkshape)
            self._dtype = numpy.dtype(chunkdtype)
            chunk = None

        self._tiled = TiledSequence(
            filesequence.shape, self._chunks, axestiled=axestiled
        )
        self._lookup = dict(
            zip(self._tiled.indices(filesequence.indices), filesequence)
        )

        zattrs = {} if zattrs is None else dict(zattrs)
        # TODO: add _ARRAY_DIMENSIONS to ZarrFileSequenceStore
        # if '_ARRAY_DIMENSIONS' not in zattrs:
        #     zattrs['_ARRAY_DIMENSIONS'] = list(...)

        self._store['.zattrs'] = ZarrStore._json(zattrs)
        self._store['.zarray'] = ZarrStore._json(
            {
                'zarr_format': 2,
                'shape': self._tiled.shape,
                'chunks': self._tiled.chunks,
                'dtype': ZarrStore._dtype_str(self._dtype),
                'compressor': None,
                'fill_value': ZarrStore._value(fillvalue, self._dtype),
                'order': 'C',
                'filters': None,
            }
        )

    def _contains(self, key: str, /) -> bool:
        """Return if key is in store."""
        try:
            indices = tuple(int(i) for i in key.split('.'))
        except Exception:
            return False
        return indices in self._lookup

    def _getitem(self, key: str, /) -> NDArray[Any]:
        """Return chunk from file."""
        indices = tuple(int(i) for i in key.split('.'))
        filename = self._lookup.get(indices, None)
        if filename is None:
            raise KeyError(key)
        return self._imread(filename, **self._kwargs)

    def _setitem(self, key: str, value: bytes, /) -> None:
        raise PermissionError('ZarrStore is read-only')

    def write_fsspec(
        self,
        jsonfile: str | os.PathLike[Any] | TextIO,
        /,
        url: str,
        *,
        quote: bool | None = None,
        groupname: str | None = None,
        templatename: str | None = None,
        codec_id: str | None = None,
        version: int | None = None,
        _append: bool = False,
        _close: bool = True,
    ) -> None:
        """Write fsspec ReferenceFileSystem as JSON to file.

        Parameters:
            jsonfile:
                Name or open file handle of output JSON file.
            url:
                Remote location of TIFF file(s) without file name(s).
            quote:
                Quote file names, that is, replace ' ' with '%20'.
                The default is True.
            groupname:
                Zarr 2 group name.
            templatename:
                Version 1 URL template name. The default is 'u'.
            codec_id:
                Name of Numcodecs codec to decode files or chunks.
            version:
                Version of fsspec file to write. The default is 0.
            _append, _close:
                Experimental API.

        References:
            - `fsspec ReferenceFileSystem format
              <https://github.com/fsspec/kerchunk>`_

        """
        from urllib.parse import quote as quote_

        kwargs = self._kwargs.copy()

        if codec_id is not None:
            pass
        elif self._imread is imread:
            codec_id = 'tifffile'
        elif 'imagecodecs' in self._imread.__module__:
            if (
                self._imread.__name__ != 'imread'
                or 'codec' not in self._kwargs
            ):
                raise ValueError('cannot determine codec_id')
            codec = kwargs.pop('codec')
            if isinstance(codec, (list, tuple)):
                codec = codec[0]
            if callable(codec):
                codec = codec.__name__.split('_')[0]
            codec_id = {
                'apng': 'imagecodecs_apng',
                'avif': 'imagecodecs_avif',
                'gif': 'imagecodecs_gif',
                'heif': 'imagecodecs_heif',
                'jpeg': 'imagecodecs_jpeg',
                'jpeg8': 'imagecodecs_jpeg',
                'jpeg12': 'imagecodecs_jpeg',
                'jpeg2k': 'imagecodecs_jpeg2k',
                'jpegls': 'imagecodecs_jpegls',
                'jpegxl': 'imagecodecs_jpegxl',
                'jpegxr': 'imagecodecs_jpegxr',
                'ljpeg': 'imagecodecs_ljpeg',
                'lerc': 'imagecodecs_lerc',
                # 'npy': 'imagecodecs_npy',
                'png': 'imagecodecs_png',
                'qoi': 'imagecodecs_qoi',
                'tiff': 'imagecodecs_tiff',
                'webp': 'imagecodecs_webp',
                'zfp': 'imagecodecs_zfp',
            }[codec]
        else:
            # TODO: choose codec from filename
            raise ValueError('cannot determine codec_id')

        if url is None:
            url = ''
        elif url and url[-1] != '/':
            url += '/'

        if groupname is None:
            groupname = ''
        elif groupname and groupname[-1] != '/':
            groupname += '/'

        refs: dict[str, Any] = {}
        if version == 1:
            if _append:
                raise ValueError('cannot append to version 1 files')
            if templatename is None:
                templatename = 'u'
            refs['version'] = 1
            refs['templates'] = {templatename: url}
            refs['gen'] = []
            refs['refs'] = refzarr = {}
            url = f'{{{{{templatename}}}}}'
        else:
            refzarr = refs

        if groupname and not _append:
            refzarr['.zgroup'] = ZarrStore._json({'zarr_format': 2}).decode()

        for key, value in self._store.items():
            if '.zarray' in key:
                value = json.loads(value)
                # TODO: make kwargs serializable
                value['compressor'] = {'id': codec_id, **kwargs}
                value = ZarrStore._json(value)
            refzarr[groupname + key] = value.decode()

        fh: TextIO
        if hasattr(jsonfile, 'write'):
            fh = jsonfile  # type: ignore[assignment]
        else:
            fh = open(jsonfile, 'w', encoding='utf-8')

        if version == 1:
            fh.write(json.dumps(refs, indent=1).rsplit('}"', 1)[0] + '}"')
            indent = '  '
        elif _append:
            fh.write(',\n')
            fh.write(json.dumps(refs, indent=1)[2:-2])
            indent = ' '
        else:
            fh.write(json.dumps(refs, indent=1)[:-2])
            indent = ' '

        prefix = len(self._commonpath)

        for key, value in self._store.items():
            if '.zarray' in key:
                value = json.loads(value)
                for index, filename in sorted(
                    self._lookup.items(), key=lambda x: x[0]
                ):
                    filename = filename[prefix:].replace('\\', '/')
                    if quote is None or quote:
                        filename = quote_(filename)
                    if filename[0] == '/':
                        filename = filename[1:]
                    indexstr = '.'.join(str(i) for i in index)
                    fh.write(
                        f',\n{indent}"{groupname}{indexstr}": '
                        f'["{url}{filename}"]'
                    )

        if version == 1:
            fh.write('\n }\n}')
        elif _close:
            fh.write('\n}')

        if not hasattr(jsonfile, 'write'):
            fh.close()

    def __enter__(self) -> ZarrFileSequenceStore:
        return self

    def __repr__(self) -> str:
        return f'<tifffile.ZarrFileSequenceStore @0x{id(self):016X}>'

    def __str__(self) -> str:
        return '\n '.join(
            (
                self.__class__.__name__,
                'shape: {}'.format(
                    ', '.join(str(i) for i in self._tiled.shape)
                ),
                'chunks: {}'.format(
                    ', '.join(str(i) for i in self._tiled.chunks)
                ),
                f'dtype: {self._dtype}',
                f'fillvalue: {self._fillvalue}',
            )
        )


class FileSequence(Sequence[str]):
    r"""Sequence of files containing compatible array data.

    Parameters:
        imread:
            Function to read image array from single file.
        files:
            Glob filename pattern or sequence of file names.
            If *None*, use '\*'.
            All files must contain array data of same shape and dtype.
            Binary streams are not supported.
        container:
            Name or open instance of ZIP file in which files are stored.
        sort:
            Function to sort file names if `files` is a pattern.
            The default is :py:func:`natural_sorted`.
            If *False*, disable sorting.
        parse:
            Function to parse sequence of sorted file names to dims, shape,
            chunk indices, and filtered file names.
            The default is :py:func:`parse_filenames` if `kwargs`
            contains `'pattern'`.
        **kwargs:
            Additional arguments passed to `parse` function.

    Examples:
        >>> filenames = ['temp_C001T002.tif', 'temp_C001T001.tif']
        >>> ims = TiffSequence(filenames, pattern=r'_(C)(\d+)(T)(\d+)')
        >>> ims[0]
        'temp_C001T002.tif'
        >>> ims.shape
        (1, 2)
        >>> ims.axes
        'CT'

    """

    imread: Callable[..., NDArray[Any]]
    """Function to read image array from single file."""

    shape: tuple[int, ...]
    """Shape of file series. Excludes shape of chunks in files."""

    axes: str
    """Character codes for dimensions in shape."""

    dims: tuple[str, ...]
    """Names of dimensions in shape."""

    indices: tuple[tuple[int, ...]]
    """Indices of files in shape."""

    _files: list[str]  # list of file names
    _container: Any  # TODO: container type?

    def __init__(
        self,
        imread: Callable[..., NDArray[Any]],
        files: (
            str | os.PathLike[Any] | Sequence[str | os.PathLike[Any]] | None
        ),
        *,
        container: str | os.PathLike[Any] | None = None,
        sort: Callable[..., Any] | bool | None = None,
        parse: Callable[..., Any] | None = None,
        **kwargs: Any,
    ) -> None:
        sort_func: Callable[..., list[str]] | None = None

        if files is None:
            files = '*'
        if sort is None:
            sort_func = natural_sorted
        elif callable(sort):
            sort_func = sort
        elif sort:
            sort_func = natural_sorted
        # elif not sort:
        #     sort_func = None

        self._container = container
        if container is not None:
            import fnmatch

            if isinstance(container, (str, os.PathLike)):
                import zipfile

                self._container = zipfile.ZipFile(container)
            elif not hasattr(self._container, 'open'):
                raise ValueError('invalid container')
            if isinstance(files, str):
                files = fnmatch.filter(self._container.namelist(), files)
                if sort_func is not None:
                    files = sort_func(files)
        elif isinstance(files, os.PathLike):
            files = [os.fspath(files)]
            if sort is not None and sort_func is not None:
                files = sort_func(files)
        elif isinstance(files, str):
            files = glob.glob(files)
            if sort_func is not None:
                files = sort_func(files)

        files = [os.fspath(f) for f in files]  # type: ignore[union-attr]
        if not files:
            raise ValueError('no files found')

        if not callable(imread):
            raise ValueError('invalid imread function')

        if container:
            # redefine imread to read from container
            def imread_(
                fname: str, _imread: Any = imread, **kwargs: Any
            ) -> NDArray[Any]:
                with self._container.open(fname) as handle1:
                    with io.BytesIO(handle1.read()) as handle2:
                        return _imread(handle2, **kwargs)

            imread = imread_

        if parse is None and kwargs.get('pattern', None):
            parse = parse_filenames

        if parse:
            try:
                dims, shape, indices, files = parse(files, **kwargs)
            except ValueError as exc:
                raise ValueError('failed to parse file names') from exc
        else:
            dims = ('sequence',)
            shape = (len(files),)
            indices = tuple((i,) for i in range(len(files)))

        assert isinstance(files, list) and isinstance(files[0], str)
        codes = TIFF.AXES_CODES
        axes = ''.join(codes.get(dim.lower(), dim[0].upper()) for dim in dims)

        self._files = files
        self.imread = imread
        self.axes = axes
        self.dims = tuple(dims)
        self.shape = tuple(shape)
        self.indices = indices

    def asarray(
        self,
        *,
        imreadargs: dict[str, Any] | None = None,
        chunkshape: tuple[int, ...] | None = None,
        chunkdtype: DTypeLike | None = None,
        axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
        out_inplace: bool | None = None,
        ioworkers: int | None = 1,
        out: OutputType = None,
        **kwargs: Any,
    ) -> NDArray[Any]:
        """Return images from files as NumPy array.

        Parameters:
            imreadargs:
                Arguments passed to :py:attr:`FileSequence.imread`.
            chunkshape:
                Shape of chunk in each file.
                Must match ``FileSequence.imread(file, **imreadargs).shape``.
                By default, this is determined by reading the first file.
            chunkdtype:
                Data type of chunk in each file.
                Must match ``FileSequence.imread(file, **imreadargs).dtype``.
                By default, this is determined by reading the first file.
            axestiled:
                Axes to be tiled.
                Map stacked sequence axis to chunk axis.
            ioworkers:
                Maximum number of threads to execute
                :py:attr:`FileSequence.imread` asynchronously.
                If *0*, use up to :py:attr:`_TIFF.MAXIOWORKERS` threads.
                Using threads can significantly improve runtime when reading
                many small files from a network share.
            out_inplace:
                :py:attr:`FileSequence.imread` decodes directly to the output
                instead of returning an array, which is copied to the output.
                Not all imread functions support this, especially in
                non-contiguous cases.
            out:
                Specifies how image array is returned.
                By default, create a new array.
                If a *numpy.ndarray*, a writable array to which the images
                are copied.
                If *'memmap'*, create a memory-mapped array in a temporary
                file.
                If a *string* or *open file*, the file used to create a
                memory-mapped array.
            **kwargs:
                Arguments passed to :py:attr:`FileSequence.imread` in
                addition to `imreadargs`.

        Raises:
            IndexError, ValueError: Array shapes do not match.

        """
        # TODO: deprecate kwargs?
        files = self._files
        if imreadargs is not None:
            kwargs |= imreadargs

        if ioworkers is None or ioworkers < 1:
            ioworkers = TIFF.MAXIOWORKERS
        ioworkers = min(len(files), ioworkers)
        assert isinstance(ioworkers, int)  # mypy bug?

        if out_inplace is None and self.imread == imread:
            out_inplace = True
        else:
            out_inplace = bool(out_inplace)

        if chunkshape is None or chunkdtype is None:
            im = self.imread(files[0], **kwargs)
            chunkshape = im.shape
            chunkdtype = im.dtype
            del im
        chunkdtype = numpy.dtype(chunkdtype)
        assert chunkshape is not None

        if axestiled:
            tiled = TiledSequence(self.shape, chunkshape, axestiled=axestiled)
            result = create_output(out, tiled.shape, chunkdtype)

            def func(index: tuple[int | slice, ...], fname: str) -> None:
                # read single image from file into result
                # if index is None:
                #     return
                if out_inplace:
                    self.imread(fname, out=result[index], **kwargs)
                else:
                    im = self.imread(fname, **kwargs)
                    result[index] = im
                    del im  # delete memory-mapped file

            if ioworkers < 2:
                for index, fname in zip(tiled.slices(self.indices), files):
                    func(index, fname)
            else:
                with ThreadPoolExecutor(ioworkers) as executor:
                    for _ in executor.map(
                        func, tiled.slices(self.indices), files
                    ):
                        pass
        else:
            shape = self.shape + chunkshape
            result = create_output(out, shape, chunkdtype)
            result = result.reshape(-1, *chunkshape)

            def func(index: tuple[int | slice, ...], fname: str) -> None:
                # read single image from file into result
                if index is None:
                    return
                index_ = int(
                    numpy.ravel_multi_index(
                        index,  # type: ignore[arg-type]
                        self.shape,
                    )
                )
                if out_inplace:
                    self.imread(fname, out=result[index_], **kwargs)
                else:
                    im = self.imread(fname, **kwargs)
                    result[index_] = im
                    del im  # delete memory-mapped file

            if ioworkers < 2:
                for index, fname in zip(self.indices, files):
                    func(index, fname)
            else:
                with ThreadPoolExecutor(ioworkers) as executor:
                    for _ in executor.map(func, self.indices, files):
                        pass

            result.shape = shape

        return result

    def aszarr(self, **kwargs: Any) -> ZarrFileSequenceStore:
        """Return images from files as Zarr 2 store.

        Parameters:
            **kwargs: Arguments passed to :py:class:`ZarrFileSequenceStore`.

        """
        return ZarrFileSequenceStore(self, **kwargs)

    def close(self) -> None:
        """Close open files."""
        if self._container is not None:
            self._container.close()
        self._container = None

    def commonpath(self) -> str:
        """Return longest common sub-path of each file in sequence."""
        if len(self._files) == 1:
            commonpath = os.path.dirname(self._files[0])
        else:
            commonpath = os.path.commonpath(self._files)
        return commonpath

    @property
    def files(self) -> list[str]:
        """Deprecated. Use the FileSequence sequence interface.

        :meta private:

        """
        warnings.warn(
            '<tifffile.FileSequence.files> is deprecated since 2024.5.22. '
            'Use the FileSequence sequence interface.',
            DeprecationWarning,
            stacklevel=2,
        )
        return self._files

    @property
    def files_missing(self) -> int:
        """Number of empty chunks."""
        return product(self.shape) - len(self._files)

    def __iter__(self) -> Iterator[str]:
        """Return iterator over all file names."""
        return iter(self._files)

    def __len__(self) -> int:
        return len(self._files)

    @overload
    def __getitem__(self, key: int, /) -> str: ...

    @overload
    def __getitem__(self, key: slice, /) -> list[str]: ...

    def __getitem__(self, key: int | slice, /) -> str | list[str]:
        return self._files[key]

    def __enter__(self) -> FileSequence:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        return f'<tifffile.FileSequence @0x{id(self):016X}>'

    def __str__(self) -> str:
        file = str(self._container) if self._container else self._files[0]
        file = os.path.split(file)[-1]
        return '\n '.join(
            (
                self.__class__.__name__,
                file,
                f'files: {len(self._files)} ({self.files_missing} missing)',
                'shape: {}'.format(', '.join(str(i) for i in self.shape)),
                'dims: {}'.format(', '.join(s for s in self.dims)),
                # f'axes: {self.axes}',
            )
        )


@final
class TiffSequence(FileSequence):
    r"""Sequence of TIFF files containing compatible array data.

    Same as :py:class:`FileSequence` with the :py:func:`imread` function,
    `'\*.tif'` glob pattern, and `out_inplace` enabled by default.

    """

    def __init__(
        self,
        files: (
            str | os.PathLike[Any] | Sequence[str | os.PathLike[Any]] | None
        ) = None,
        *,
        imread: Callable[..., NDArray[Any]] = imread,
        **kwargs: Any,
    ) -> None:
        super().__init__(imread, '*.tif' if files is None else files, **kwargs)

    def __repr__(self) -> str:
        return f'<tifffile.TiffSequence @0x{id(self):016X}>'


@final
class TiledSequence:
    """Tiled sequence of chunks.

    Transform a sequence of stacked chunks to tiled chunks.

    Parameters:
        stackshape:
            Shape of stacked sequence excluding chunks.
        chunkshape:
            Shape of chunks.
        axestiled:
            Axes to be tiled. Map stacked sequence axis
            to chunk axis. By default, the sequence is not tiled.
        axes:
            Character codes for dimensions in stackshape and chunkshape.

    Examples:
        >>> ts = TiledSequence((1, 2), (3, 4), axestiled={1: 0}, axes='ABYX')
        >>> ts.shape
        (1, 6, 4)
        >>> ts.chunks
        (1, 3, 4)
        >>> ts.axes
        'AYX'

    """

    chunks: tuple[int, ...]
    """Shape of chunks in tiled sequence."""
    # with same number of dimensions as shape

    shape: tuple[int, ...]
    """Shape of tiled sequence including chunks."""

    axes: str | tuple[str, ...] | None
    """Dimensions codes of tiled sequence."""

    shape_squeezed: tuple[int, ...]
    """Shape of tiled sequence with length-1 dimensions removed."""

    axes_squeezed: str | tuple[str, ...] | None
    """Dimensions codes of tiled sequence with length-1 dimensions removed."""

    _stackdims: int
    """Number of dimensions in stack excluding chunks."""

    _chunkdims: int
    """Number of dimensions in chunks."""

    _shape_untiled: tuple[int, ...]
    """Shape of untiled sequence (stackshape + chunkshape)."""

    _axestiled: tuple[tuple[int, int], ...]
    """Map axes to tile from stack to chunks."""

    def __init__(
        self,
        stackshape: Sequence[int],
        chunkshape: Sequence[int],
        /,
        *,
        axestiled: dict[int, int] | Sequence[tuple[int, int]] | None = None,
        axes: str | Sequence[str] | None = None,
    ) -> None:
        self._stackdims = len(stackshape)
        self._chunkdims = len(chunkshape)
        self._shape_untiled = tuple(stackshape) + tuple(chunkshape)
        if axes is not None and len(axes) != len(self._shape_untiled):
            raise ValueError(
                'axes length does not match stackshape + chunkshape'
            )

        if axestiled:
            axestiled = dict(axestiled)
            for ax0, ax1 in axestiled.items():
                axestiled[ax0] = ax1 + self._stackdims
            self._axestiled = tuple(reversed(sorted(axestiled.items())))

            axes_list = [] if axes is None else list(axes)
            shape = list(self._shape_untiled)
            chunks = [1] * self._stackdims + list(chunkshape)
            used = set()
            for ax0, ax1 in self._axestiled:
                if ax0 in used or ax1 in used:
                    raise ValueError('duplicate axis')
                used.add(ax0)
                used.add(ax1)
                shape[ax1] *= stackshape[ax0]
            for ax0, ax1 in self._axestiled:
                del shape[ax0]
                del chunks[ax0]
                if axes_list:
                    del axes_list[ax0]
            self.shape = tuple(shape)
            self.chunks = tuple(chunks)
            if axes is None:
                self.axes = None
            elif isinstance(axes, str):
                self.axes = ''.join(axes_list)
            else:
                self.axes = tuple(axes_list)
        else:
            self._axestiled = ()
            self.shape = self._shape_untiled
            self.chunks = (1,) * self._stackdims + tuple(chunkshape)
            if axes is None:
                self.axes = None
            elif isinstance(axes, str):
                self.axes = axes
            else:
                self.axes = tuple(axes)

        assert len(self.shape) == len(self.chunks)
        if self.axes is not None:
            assert len(self.shape) == len(self.axes)

        if self.axes is None:
            self.shape_squeezed = tuple(i for i in self.shape if i > 1)
            self.axes_squeezed = None
        else:
            keep = ('X', 'Y', 'width', 'length', 'height')
            self.shape_squeezed = tuple(
                i
                for i, ax in zip(self.shape, self.axes)
                if i > 1 or ax in keep
            )
            squeezed = tuple(
                ax
                for i, ax in zip(self.shape, self.axes)
                if i > 1 or ax in keep
            )
            self.axes_squeezed = (
                ''.join(squeezed) if isinstance(self.axes, str) else squeezed
            )

    def indices(
        self, indices: Iterable[Sequence[int]], /
    ) -> Iterator[tuple[int, ...]]:
        """Return iterator over chunk indices of tiled sequence.

        Parameters:
            indices: Indices of chunks in stacked sequence.

        """
        chunkindex = [0] * self._chunkdims
        for index in indices:
            if index is None:
                yield None
            else:
                if len(index) != self._stackdims:
                    raise ValueError(f'{len(index)} != {self._stackdims}')
                index = list(index) + chunkindex
                for ax0, ax1 in self._axestiled:
                    index[ax1] = index[ax0]
                for ax0, ax1 in self._axestiled:
                    del index[ax0]
                yield tuple(index)

    def slices(
        self, indices: Iterable[Sequence[int]] | None = None, /
    ) -> Iterator[tuple[int | slice, ...]]:
        """Return iterator over slices of chunks in tiled sequence.

        Parameters:
            indices: Indices of chunks in stacked sequence.

        """
        wholeslice: list[int | slice]
        chunkslice: list[int | slice] = [slice(None)] * self._chunkdims

        if indices is None:
            indices = numpy.ndindex(self._shape_untiled[: self._stackdims])

        for index in indices:
            if index is None:
                yield None
            else:
                assert len(index) == self._stackdims
                wholeslice = [*index, *chunkslice]
                for ax0, ax1 in self._axestiled:
                    j = self._shape_untiled[ax1]
                    i = cast(int, wholeslice[ax0]) * j
                    wholeslice[ax1] = slice(i, i + j)
                for ax0, ax1 in self._axestiled:
                    del wholeslice[ax0]
                yield tuple(wholeslice)

    @property
    def ndim(self) -> int:
        """Number of dimensions of tiled sequence excluding chunks."""
        return len(self.shape)

    @property
    def is_tiled(self) -> bool:
        """Sequence is tiled."""
        return bool(self._axestiled)


@final
class FileHandle:
    """Binary file handle.

    A limited, special purpose binary file handle that can:

    - handle embedded files (for example, LSM within LSM files).
    - re-open closed files (for multi-file formats, such as OME-TIFF).
    - read and write NumPy arrays and records from file-like objects.

    When initialized from another file handle, do not use the other handle
    unless this FileHandle is closed.

    FileHandle instances are not thread-safe.

    Parameters:
        file:
            File name or seekable binary stream, such as open file,
            BytesIO, or fsspec OpenFile.
        mode:
            File open mode if `file` is file name.
            The default is 'rb'. Files are always opened in binary mode.
        name:
            Name of file if `file` is binary stream.
        offset:
            Start position of embedded file.
            The default is the current file position.
        size:
            Size of embedded file.
            The default is the number of bytes from `offset` to
            the end of the file.

    """

    # TODO: make FileHandle a subclass of IO[bytes]

    __slots__ = (
        '_fh',
        '_file',
        '_mode',
        '_name',
        '_dir',
        '_lock',
        '_offset',
        '_size',
        '_close',
    )

    _file: str | os.PathLike[Any] | FileHandle | IO[bytes] | None
    _fh: IO[bytes] | None
    _mode: str
    _name: str
    _dir: str
    _offset: int
    _size: int
    _close: bool
    _lock: threading.RLock | NullContext

    def __init__(
        self,
        file: str | os.PathLike[Any] | FileHandle | IO[bytes],
        /,
        mode: (
            Literal['r', 'r+', 'w', 'x', 'rb', 'r+b', 'wb', 'xb'] | None
        ) = None,
        *,
        name: str | None = None,
        offset: int | None = None,
        size: int | None = None,
    ) -> None:
        self._mode = 'rb' if mode is None else mode
        self._fh = None
        self._file = file  # reference to original argument for re-opening
        self._name = name if name else ''
        self._dir = ''
        self._offset = -1 if offset is None else offset
        self._size = -1 if size is None else size
        self._close = True
        self._lock = NullContext()
        self.open()
        assert self._fh is not None

    def open(self) -> None:
        """Open or re-open file."""
        if self._fh is not None:
            return  # file is open

        if isinstance(self._file, os.PathLike):
            self._file = os.fspath(self._file)

        if isinstance(self._file, str):
            # file name
            if self._mode[-1:] != 'b':
                self._mode += 'b'
            if self._mode not in {'rb', 'r+b', 'wb', 'xb'}:
                raise ValueError(f'invalid mode {self._mode}')
            self._file = os.path.realpath(self._file)
            self._dir, self._name = os.path.split(self._file)
            self._fh = open(self._file, self._mode, encoding=None)
            self._close = True
            self._offset = max(0, self._offset)
        elif isinstance(self._file, FileHandle):
            # FileHandle
            self._fh = self._file._fh
            self._offset = max(0, self._offset)
            self._offset += self._file._offset
            self._close = False
            if not self._name:
                if self._offset:
                    name, ext = os.path.splitext(self._file._name)
                    self._name = f'{name}@{self._offset}{ext}'
                else:
                    self._name = self._file._name
            self._mode = self._file._mode
            self._dir = self._file._dir
        elif hasattr(self._file, 'seek'):
            # binary stream: open file, BytesIO, fsspec LocalFileOpener
            # cast to IO[bytes] even it might not be
            self._fh = cast(IO[bytes], self._file)
            try:
                self._fh.tell()
            except Exception as exc:
                raise ValueError('binary stream is not seekable') from exc

            if self._offset < 0:
                self._offset = self._fh.tell()
            self._close = False
            if not self._name:
                try:
                    self._dir, self._name = os.path.split(self._fh.name)
                except AttributeError:
                    try:
                        self._dir, self._name = os.path.split(
                            self._fh.path  # type: ignore[attr-defined]
                        )
                    except AttributeError:
                        self._name = 'Unnamed binary stream'
            try:
                self._mode = self._fh.mode
            except AttributeError:
                pass
        elif hasattr(self._file, 'open'):
            # fsspec OpenFile
            _file: Any = self._file
            self._fh = cast(IO[bytes], _file.open())
            try:
                self._fh.tell()
            except Exception as exc:
                try:
                    self._fh.close()
                except Exception:
                    pass
                raise ValueError('OpenFile is not seekable') from exc

            if self._offset < 0:
                self._offset = self._fh.tell()
            self._close = True
            if not self._name:
                try:
                    self._dir, self._name = os.path.split(_file.path)
                except AttributeError:
                    self._name = 'Unnamed binary stream'
            try:
                self._mode = _file.mode
            except AttributeError:
                pass

        else:
            raise ValueError(
                'the first parameter must be a file name '
                'or seekable binary file object, '
                f'not {type(self._file)!r}'
            )

        assert self._fh is not None

        if self._offset:
            self._fh.seek(self._offset)

        if self._size < 0:
            pos = self._fh.tell()
            self._fh.seek(self._offset, os.SEEK_END)
            self._size = self._fh.tell()
            self._fh.seek(pos)

    def close(self) -> None:
        """Close file handle."""
        if self._close and self._fh is not None:
            try:
                self._fh.close()
            except Exception:
                # PermissionError on MacOS. See issue #184
                pass
            self._fh = None

    def fileno(self) -> int:
        """Return underlying file descriptor if exists, else raise OSError."""
        assert self._fh is not None
        try:
            return self._fh.fileno()
        except (OSError, AttributeError) as exc:
            raise OSError(
                f'{type(self._fh)} does not have a file descriptor'
            ) from exc

    def writable(self) -> bool:
        """Return True if stream supports writing."""
        assert self._fh is not None
        if hasattr(self._fh, 'writable'):
            return self._fh.writable()
        return False

    def seekable(self) -> bool:
        """Return True if stream supports random access."""
        return True

    def tell(self) -> int:
        """Return file's current position."""
        assert self._fh is not None
        return self._fh.tell() - self._offset

    def seek(self, offset: int, /, whence: int = 0) -> int:
        """Set file's current position.

        Parameters:
            offset:
                Position of file handle relative to position indicated
                by `whence`.
            whence:
                Relative position of `offset`.
                0 (`os.SEEK_SET`) beginning of file (default).
                1 (`os.SEEK_CUR`) current position.
                2 (`os.SEEK_END`) end of file.

        """
        assert self._fh is not None
        if self._offset:
            if whence == 0:
                return (
                    self._fh.seek(self._offset + offset, whence) - self._offset
                )
            if whence == 2 and self._size > 0:
                return (
                    self._fh.seek(self._offset + self._size + offset, 0)
                    - self._offset
                )
        return self._fh.seek(offset, whence)

    def read(self, size: int = -1, /) -> bytes:
        """Return bytes read from file.

        Parameters:
            size:
                Number of bytes to read from file.
                By default, read until the end of the file.

        """
        if size < 0 and self._offset:
            size = self._size
        assert self._fh is not None
        return self._fh.read(size)

    def readinto(self, buffer: bytes, /) -> int:
        """Read bytes from file into buffer.

        Parameters:
            buffer: Buffer to read into.

        Returns:
            Number of bytes read from file.

        """
        assert self._fh is not None
        return self._fh.readinto(buffer)  # type: ignore[attr-defined]

    def write(self, buffer: bytes, /) -> int:
        """Write bytes to file and return number of bytes written.

        Parameters:
            buffer: Bytes to write to file.

        Returns:
            Number of bytes written.

        """
        assert self._fh is not None
        return self._fh.write(buffer)

    def flush(self) -> None:
        """Flush write buffers of stream if applicable."""
        assert self._fh is not None
        if hasattr(self._fh, 'flush'):
            self._fh.flush()

    def memmap_array(
        self,
        dtype: DTypeLike,
        shape: tuple[int, ...],
        offset: int = 0,
        *,
        mode: str = 'r',
        order: str = 'C',
    ) -> NDArray[Any]:
        """Return `numpy.memmap` of array data stored in file.

        Parameters:
            dtype:
                Data type of array in file.
            shape:
                Shape of array in file.
            offset:
                Start position of array-data in file.
            mode:
                File is opened in this mode. The default is read-only.
            order:
                Order of ndarray memory layout. The default is 'C'.

        """
        if not self.is_file:
            raise ValueError('cannot memory-map file without fileno')
        assert self._fh is not None
        return numpy.memmap(
            self._fh,  # type: ignore[call-overload]
            dtype=dtype,
            mode=mode,
            offset=self._offset + offset,
            shape=shape,
            order=order,
        )

    def read_array(
        self,
        dtype: DTypeLike,
        count: int = -1,
        offset: int = 0,
        *,
        out: NDArray[Any] | None = None,
    ) -> NDArray[Any]:
        """Return NumPy array from file in native byte order.

        Parameters:
            dtype:
                Data type of array to read.
            count:
                Number of items to read. By default, all items are read.
            offset:
                Start position of array-data in file.
            out:
                NumPy array to read into. By default, a new array is created.

        """
        dtype = numpy.dtype(dtype)

        if count < 0:
            nbytes = self._size if out is None else out.nbytes
            count = nbytes // dtype.itemsize
        else:
            nbytes = count * dtype.itemsize

        result = numpy.empty(count, dtype) if out is None else out

        if result.nbytes != nbytes:
            raise ValueError('size mismatch')

        assert self._fh is not None

        if offset:
            self._fh.seek(self._offset + offset)

        try:
            n = self._fh.readinto(result)  # type: ignore[attr-defined]
        except AttributeError:
            result[:] = numpy.frombuffer(self._fh.read(nbytes), dtype).reshape(
                result.shape
            )
            n = nbytes

        if n != nbytes:
            raise ValueError(f'failed to read {nbytes} bytes, got {n}')

        if not result.dtype.isnative:
            if not dtype.isnative:
                result.byteswap(True)
            result = result.view(result.dtype.newbyteorder())
        elif result.dtype.isnative != dtype.isnative:
            result.byteswap(True)

        if out is not None:
            if hasattr(out, 'flush'):
                out.flush()

        return result

    def read_record(
        self,
        dtype: DTypeLike,
        shape: tuple[int, ...] | int | None = 1,
        *,
        byteorder: Literal['S', '<', '>', '=', '|'] | None = None,
    ) -> numpy.recarray[Any, Any]:
        """Return NumPy record from file.

        Parameters:
            dtype:
                Data type of record array to read.
            shape:
                Shape of record array to read.
            byteorder:
                Byte order of record array to read.

        """
        assert self._fh is not None

        dtype = numpy.dtype(dtype)
        if byteorder is not None:
            dtype = dtype.newbyteorder(byteorder)

        try:
            record = numpy.rec.fromfile(  # type: ignore[call-overload]
                self._fh, dtype, shape
            )
        except Exception:
            if shape is None:
                shape = self._size // dtype.itemsize
            size = product(sequence(shape)) * dtype.itemsize
            # data = bytearray(size)
            # n = self._fh.readinto(data)
            # data = data[:n]
            # TODO: record is not writable
            data = self._fh.read(size)
            record = numpy.rec.fromstring(
                data,
                dtype,
                shape,
            )
        return record[0] if shape == 1 else record

    def write_empty(self, size: int, /) -> int:
        """Append null-bytes to file.

        The file position must be at the end of the file.

        Parameters:
            size: Number of null-bytes to write to file.

        """
        if size < 1:
            return 0
        assert self._fh is not None
        self._fh.seek(size - 1, os.SEEK_CUR)
        self._fh.write(b'\x00')
        return size

    def write_array(
        self,
        data: NDArray[Any],
        dtype: DTypeLike = None,
        /,
    ) -> int:
        """Write NumPy array to file in C contiguous order.

        Parameters:
            data: Array to write to file.

        """
        assert self._fh is not None
        pos = self._fh.tell()
        # writing non-contiguous arrays is very slow
        data = numpy.ascontiguousarray(data, dtype)
        try:
            data.tofile(self._fh)
        except io.UnsupportedOperation:
            # numpy cannot write to BytesIO
            self._fh.write(data.tobytes())
        return self._fh.tell() - pos

    def read_segments(
        self,
        offsets: Sequence[int],
        bytecounts: Sequence[int],
        /,
        indices: Sequence[int] | None = None,
        *,
        sort: bool = True,
        lock: threading.RLock | NullContext | None = None,
        buffersize: int | None = None,
        flat: bool = True,
    ) -> (
        Iterator[tuple[bytes | None, int]]
        | Iterator[list[tuple[bytes | None, int]]]
    ):
        """Return iterator over segments read from file and their indices.

        The purpose of this function is to

        - reduce small or random reads.
        - reduce acquiring reentrant locks.
        - synchronize seeks and reads.
        - limit size of segments read into memory at once.
          (ThreadPoolExecutor.map is not collecting iterables lazily).

        Parameters:
            offsets:
                Offsets of segments to read from file.
            bytecounts:
                Byte counts of segments to read from file.
            indices:
                Indices of segments in image.
                The default is `range(len(offsets))`.
            sort:
                Read segments from file in order of their offsets.
            lock:
                Reentrant lock to synchronize seeks and reads.
            buffersize:
                Approximate number of bytes to read from file in one pass.
                The default is :py:attr:`_TIFF.BUFFERSIZE`.
            flat:
                If *True*, return iterator over individual (segment, index)
                tuples.
                Else, return an iterator over a list of (segment, index)
                tuples that were acquired in one pass.

        Yields:
            Individual or lists of `(segment, index)` tuples.

        """
        # TODO: Cythonize this?
        assert self._fh is not None
        length = len(offsets)
        if length < 1:
            return
        if length == 1:
            index = 0 if indices is None else indices[0]
            if bytecounts[index] > 0 and offsets[index] > 0:
                if lock is None:
                    lock = self._lock
                with lock:
                    self.seek(offsets[index])
                    data = self._fh.read(bytecounts[index])
            else:
                data = None
            yield (data, index) if flat else [(data, index)]
            return

        if lock is None:
            lock = self._lock
        if buffersize is None:
            buffersize = TIFF.BUFFERSIZE

        if indices is None:
            segments = [(i, offsets[i], bytecounts[i]) for i in range(length)]
        else:
            segments = [
                (indices[i], offsets[i], bytecounts[i]) for i in range(length)
            ]
        if sort:
            segments = sorted(segments, key=lambda x: x[1])

        iscontig = True
        for i in range(length - 1):
            _, offset, bytecount = segments[i]
            nextoffset = segments[i + 1][1]
            if offset == 0 or bytecount == 0 or nextoffset == 0:
                continue
            if offset + bytecount != nextoffset:
                iscontig = False
                break

        seek = self.seek
        read = self._fh.read
        result: list[tuple[bytes | None, int]]

        if iscontig:
            # consolidate reads
            i = 0
            while i < length:
                j = i
                offset = -1
                bytecount = 0
                while bytecount <= buffersize and i < length:
                    _, o, b = segments[i]
                    if o > 0 and b > 0:
                        if offset < 0:
                            offset = o
                        bytecount += b
                    i += 1

                if offset < 0:
                    data = None
                else:
                    with lock:
                        seek(offset)
                        data = read(bytecount)
                start = 0
                stop = 0
                result = []
                while j < i:
                    index, offset, bytecount = segments[j]
                    if offset > 0 and bytecount > 0:
                        stop += bytecount
                        result.append(
                            (data[start:stop], index)  # type: ignore[index]
                        )
                        start = stop
                    else:
                        result.append((None, index))
                    j += 1
                if flat:
                    yield from result
                else:
                    yield result
            return

        i = 0
        while i < length:
            result = []
            size = 0
            with lock:
                while size <= buffersize and i < length:
                    index, offset, bytecount = segments[i]
                    if offset > 0 and bytecount > 0:
                        seek(offset)
                        result.append((read(bytecount), index))
                        # buffer = bytearray(bytecount)
                        # n = fh.readinto(buffer)
                        # data.append(buffer[:n])
                        size += bytecount
                    else:
                        result.append((None, index))
                    i += 1
            if flat:
                yield from result
            else:
                yield result

    def __enter__(self) -> FileHandle:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()
        self._file = None

    # TODO: this may crash the Python interpreter under certain conditions
    # def __getattr__(self, name: str, /) -> Any:
    #     """Return attribute from underlying file object."""
    #     if self._offset:
    #         warnings.warn(
    #             '<tifffile.FileHandle> '
    #             f'{name} not implemented for embedded files',
    #             UserWarning,
    #         )
    #     return getattr(self._fh, name)

    def __repr__(self) -> str:
        return f'<tifffile.FileHandle {snipstr(self._name, 32)!r}>'

    def __str__(self) -> str:
        return '\n '.join(
            (
                'FileHandle',
                self._name,
                self._dir,
                f'{self._size} bytes',
                'closed' if self._fh is None else 'open',
            )
        )

    @property
    def name(self) -> str:
        """Name of file or stream."""
        return self._name

    @property
    def dirname(self) -> str:
        """Directory in which file is stored."""
        return self._dir

    @property
    def path(self) -> str:
        """Absolute path of file."""
        return os.path.join(self._dir, self._name)

    @property
    def extension(self) -> str:
        """File name extension of file or stream."""
        name, ext = os.path.splitext(self._name.lower())
        if ext and name.endswith('.ome'):
            ext = '.ome' + ext
        return ext

    @property
    def size(self) -> int:
        """Size of file in bytes."""
        return self._size

    @property
    def closed(self) -> bool:
        """File is closed."""
        return self._fh is None

    @property
    def lock(self) -> threading.RLock | NullContext:
        """Reentrant lock to synchronize reads and writes."""
        return self._lock

    @lock.setter
    def lock(self, value: bool, /) -> None:
        self.set_lock(value)

    def set_lock(self, value: bool, /) -> None:
        if bool(value) == isinstance(self._lock, NullContext):
            self._lock = threading.RLock() if value else NullContext()

    @property
    def has_lock(self) -> bool:
        """A reentrant lock is currently used to sync reads and writes."""
        return not isinstance(self._lock, NullContext)

    @property
    def is_file(self) -> bool:
        """File has fileno and can be memory-mapped."""
        try:
            self._fh.fileno()  # type: ignore[union-attr]
            return True
        except Exception:
            return False


@final
class FileCache:
    """Keep FileHandles open.

    Parameters:
        size: Maximum number of files to keep open. The default is 8.
        lock: Reentrant lock to synchronize reads and writes.

    """

    __slots__ = ('files', 'keep', 'past', 'lock', 'size')

    size: int
    """Maximum number of files to keep open."""

    files: dict[FileHandle, int]
    """Reference counts of opened files."""

    keep: set[FileHandle]
    """Set of files to keep open."""

    past: list[FileHandle]
    """FIFO list of opened files."""

    lock: threading.RLock | NullContext
    """Reentrant lock to synchronize reads and writes."""

    def __init__(
        self,
        size: int | None = None,
        *,
        lock: threading.RLock | NullContext | None = None,
    ) -> None:
        self.past = []
        self.files = {}
        self.keep = set()
        self.size = 8 if size is None else int(size)
        self.lock = NullContext() if lock is None else lock

    def open(self, fh: FileHandle, /) -> None:
        """Open file, re-open if necessary."""
        with self.lock:
            if fh in self.files:
                self.files[fh] += 1
            elif fh.closed:
                fh.open()
                self.files[fh] = 1
                self.past.append(fh)
            else:
                self.files[fh] = 2
                self.keep.add(fh)
                self.past.append(fh)

    def close(self, fh: FileHandle, /) -> None:
        """Close least recently used open files."""
        with self.lock:
            if fh in self.files:
                self.files[fh] -= 1
            self._trim()

    def clear(self) -> None:
        """Close all opened files if not in use when opened first."""
        with self.lock:
            for fh, refcount in list(self.files.items()):
                if fh not in self.keep:
                    fh.close()
                    del self.files[fh]
                    del self.past[self.past.index(fh)]

    def read(
        self,
        fh: FileHandle,
        /,
        offset: int,
        bytecount: int,
        whence: int = 0,
    ) -> bytes:
        """Return bytes read from binary file.

        Parameters:
            fh:
                File handle to read from.
            offset:
                Position in file to start reading from relative to the
                position indicated by `whence`.
            bytecount:
                Number of bytes to read.
            whence:
                Relative position of offset.
                0 (`os.SEEK_SET`) beginning of file (default).
                1 (`os.SEEK_CUR`) current position.
                2 (`os.SEEK_END`) end of file.

        """
        # this function is more efficient than
        # filecache.open(fh)
        # with lock:
        #     fh.seek()
        #     data = fh.read()
        # filecache.close(fh)
        with self.lock:
            b = fh not in self.files
            if b:
                if fh.closed:
                    fh.open()
                    self.files[fh] = 0
                else:
                    self.files[fh] = 1
                    self.keep.add(fh)
                self.past.append(fh)
            fh.seek(offset, whence)
            data = fh.read(bytecount)
            if b:
                self._trim()
        return data

    def write(
        self,
        fh: FileHandle,
        /,
        offset: int,
        data: bytes,
        whence: int = 0,
    ) -> int:
        """Write bytes to binary file.

        Parameters:
            fh:
                File handle to write to.
            offset:
                Position in file to start writing from relative to the
                position indicated by `whence`.
            value:
                Bytes to write.
            whence:
                Relative position of offset.
                0 (`os.SEEK_SET`) beginning of file (default).
                1 (`os.SEEK_CUR`) current position.
                2 (`os.SEEK_END`) end of file.

        """
        with self.lock:
            b = fh not in self.files
            if b:
                if fh.closed:
                    fh.open()
                    self.files[fh] = 0
                else:
                    self.files[fh] = 1
                    self.keep.add(fh)
                self.past.append(fh)
            fh.seek(offset, whence)
            written = fh.write(data)
            if b:
                self._trim()
        return written

    def _trim(self) -> None:
        """Trim file cache."""
        index = 0
        size = len(self.past)
        while index < size > self.size:
            fh = self.past[index]
            if fh not in self.keep and self.files[fh] <= 0:
                fh.close()
                del self.files[fh]
                del self.past[index]
                size -= 1
            else:
                index += 1

    def __len__(self) -> int:
        """Return number of open files."""
        return len(self.files)

    def __repr__(self) -> str:
        return f'<tifffile.FileCache @0x{id(self):016X}>'


@final
class StoredShape(Sequence[int]):
    """Normalized shape of image array in TIFF pages.

    Parameters:
        frames:
            Number of TIFF pages.
        separate_samples:
            Number of separate samples.
        depth:
            Image depth.
        length:
            Image length (height).
        width:
            Image width.
        contig_samples:
            Number of contiguous samples.
        extrasamples:
            Number of extra samples.

    """

    __slots__ = (
        'frames',
        'separate_samples',
        'depth',
        'length',
        'width',
        'contig_samples',
        'extrasamples',
    )

    frames: int
    """Number of TIFF pages."""

    separate_samples: int
    """Number of separate samples."""

    depth: int
    """Image depth. Value of ImageDepth tag or 1."""

    length: int
    """Image length (height). Value of ImageLength tag."""

    width: int
    """Image width. Value of ImageWidth tag."""

    contig_samples: int
    """Number of contiguous samples."""

    extrasamples: int
    """Number of extra samples. Count of ExtraSamples tag or 0."""

    def __init__(
        self,
        frames: int = 1,
        separate_samples: int = 1,
        depth: int = 1,
        length: int = 1,
        width: int = 1,
        contig_samples: int = 1,
        extrasamples: int = 0,
    ) -> None:
        if separate_samples != 1 and contig_samples != 1:
            raise ValueError('invalid samples')

        self.frames = int(frames)
        self.separate_samples = int(separate_samples)
        self.depth = int(depth)
        self.length = int(length)
        self.width = int(width)
        self.contig_samples = int(contig_samples)
        self.extrasamples = int(extrasamples)

    @property
    def size(self) -> int:
        """Product of all dimensions."""
        return (
            abs(self.frames)
            * self.separate_samples
            * self.depth
            * self.length
            * self.width
            * self.contig_samples
        )

    @property
    def samples(self) -> int:
        """Number of samples. Count of SamplesPerPixel tag."""
        assert self.separate_samples == 1 or self.contig_samples == 1
        samples = (
            self.separate_samples
            if self.separate_samples > 1
            else self.contig_samples
        )
        assert self.extrasamples < samples
        return samples

    @property
    def photometric_samples(self) -> int:
        """Number of photometric samples."""
        return self.samples - self.extrasamples

    @property
    def shape(self) -> tuple[int, int, int, int, int, int]:
        """Normalized 6D shape of image array in all pages."""
        return (
            self.frames,
            self.separate_samples,
            self.depth,
            self.length,
            self.width,
            self.contig_samples,
        )

    @property
    def page_shape(self) -> tuple[int, int, int, int, int]:
        """Normalized 5D shape of image array in single page."""
        return (
            self.separate_samples,
            self.depth,
            self.length,
            self.width,
            self.contig_samples,
        )

    @property
    def page_size(self) -> int:
        """Product of dimensions in single page."""
        return (
            self.separate_samples
            * self.depth
            * self.length
            * self.width
            * self.contig_samples
        )

    @property
    def squeezed(self) -> tuple[int, ...]:
        """Shape with length-1 removed, except for length and width."""
        shape = [self.length, self.width]
        if self.separate_samples > 1:
            shape.insert(0, self.separate_samples)
        elif self.contig_samples > 1:
            shape.append(self.contig_samples)
        if self.frames > 1:
            shape.insert(0, self.frames)
        return tuple(shape)

    @property
    def is_valid(self) -> bool:
        """Shape is valid."""
        return (
            self.frames >= 1
            and self.depth >= 1
            and self.length >= 1
            and self.width >= 1
            and (self.separate_samples == 1 or self.contig_samples == 1)
            and (
                self.contig_samples
                if self.contig_samples > 1
                else self.separate_samples
            )
            > self.extrasamples
        )

    @property
    def is_planar(self) -> bool:
        """Shape contains planar samples."""
        return self.separate_samples > 1

    @property
    def planarconfig(self) -> int | None:
        """Value of PlanarConfiguration tag."""
        if self.separate_samples > 1:
            return 2  # PLANARCONFIG.SEPARATE
        if self.contig_samples > 1:
            return 1  # PLANARCONFIG.CONTIG
        return None

    def __len__(self) -> int:
        return 6

    @overload
    def __getitem__(self, key: int, /) -> int: ...

    @overload
    def __getitem__(self, key: slice, /) -> tuple[int, ...]: ...

    def __getitem__(self, key: int | slice, /) -> int | tuple[int, ...]:
        return (
            self.frames,
            self.separate_samples,
            self.depth,
            self.length,
            self.width,
            self.contig_samples,
        )[key]

    def __eq__(self, other: object, /) -> bool:
        return (
            isinstance(other, StoredShape)
            and self.frames == other.frames
            and self.separate_samples == other.separate_samples
            and self.depth == other.depth
            and self.length == other.length
            and self.width == other.width
            and self.contig_samples == other.contig_samples
        )

    def __repr__(self) -> str:
        return (
            '<StoredShape('
            f'frames={self.frames}, '
            f'separate_samples={self.separate_samples}, '
            f'depth={self.depth}, '
            f'length={self.length}, '
            f'width={self.width}, '
            f'contig_samples={self.contig_samples}, '
            f'extrasamples={self.extrasamples}'
            ')>'
        )


@final
class NullContext:
    """Null context manager. Can be used as a dummy reentrant lock.

    >>> with NullContext():
    ...     pass
    ...

    """

    __slots__ = ()

    def __enter__(self) -> NullContext:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        pass

    def __repr__(self) -> str:
        return 'NullContext()'


@final
class Timer:
    """Stopwatch for timing execution speed.

    Parameters:
        message:
            Message to print.
        end:
            End of print statement.
        started:
            Value of performance counter when started.
            The default is the current performance counter.

    Examples:
        >>> import time
        >>> with Timer('sleep:'):
        ...     time.sleep(1.05)
        sleep: 1.0... s

    """

    __slots__ = ('started', 'stopped', 'duration')

    started: float
    """Value of performance counter when started."""

    stopped: float
    """Value of performance counter when stopped."""

    duration: float
    """Duration between `started` and `stopped` in seconds."""

    def __init__(
        self,
        message: str | None = None,
        *,
        end: str = ' ',
        started: float | None = None,
    ) -> None:
        if message is not None:
            print(message, end=end, flush=True)
        self.duration = 0.0
        if started is None:
            started = time.perf_counter()
        self.started = self.stopped = started

    def start(self, message: str | None = None, *, end: str = ' ') -> float:
        """Start timer and return current time."""
        if message is not None:
            print(message, end=end, flush=True)
        self.duration = 0.0
        self.started = self.stopped = time.perf_counter()
        return self.started

    def stop(self, message: str | None = None, *, end: str = ' ') -> float:
        """Return duration of timer till start.

        Parameters:
            message: Message to print.
            end: End of print statement.

        """
        self.stopped = time.perf_counter()
        if message is not None:
            print(message, end=end, flush=True)
        self.duration = self.stopped - self.started
        return self.duration

    def print(
        self, message: str | None = None, *, end: str | None = None
    ) -> None:
        """Print duration from timer start till last stop or now.

        Parameters:
            message: Message to print.
            end: End of print statement.

        """
        msg = str(self)
        if message is not None:
            print(message, end=' ')
        print(msg, end=end, flush=True)

    @staticmethod
    def clock() -> float:
        """Return value of performance counter."""
        return time.perf_counter()

    def __str__(self) -> str:
        """Return duration from timer start till last stop or now."""
        if self.duration <= 0.0:
            # not stopped
            duration = time.perf_counter() - self.started
        else:
            duration = self.duration
        s = str(TimeDelta(seconds=duration))
        i = 0
        while i < len(s) and s[i : i + 2] in '0:0010203040506070809':
            i += 1
        if s[i : i + 1] == ':':
            i += 1
        return f'{s[i:]} s'

    def __repr__(self) -> str:
        return f'Timer(started={self.started})'

    def __enter__(self) -> Timer:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.print()


class OmeXmlError(Exception):
    """Exception to indicate invalid OME-XML or unsupported cases."""


@final
class OmeXml:
    """Create OME-TIFF XML metadata.

    Parameters:
        **metadata:
            Additional OME-XML attributes or elements to be stored.

            Creator:
                Name of creating application. The default is 'tifffile'.
            UUID:
                Unique identifier.

    Examples:
        >>> omexml = OmeXml()
        >>> omexml.addimage(
        ...     dtype='uint16',
        ...     shape=(32, 256, 256),
        ...     storedshape=(32, 1, 1, 256, 256, 1),
        ...     axes='CYX',
        ...     Name='First Image',
        ...     PhysicalSizeX=2.0,
        ...     MapAnnotation={'key': 'value'},
        ...     Dataset={'Name': 'FirstDataset'},
        ... )
        >>> xml = omexml.tostring()
        >>> xml
        '<OME ...<Image ID="Image:0" Name="First Image">...</Image>...</OME>'
        >>> OmeXml.validate(xml)
        True

    """

    images: list[str]
    """OME-XML Image elements."""

    annotations: list[str]
    """OME-XML Annotation elements."""

    datasets: list[str]
    """OME-XML Dataset elements."""

    _xml: str
    _ifd: int

    def __init__(self, **metadata: Any) -> None:
        metadata = metadata.get('OME', metadata)

        self._ifd = 0
        self.images = []
        self.annotations = []
        self.datasets = []
        # TODO: parse other OME elements from metadata
        #   Project
        #   Folder
        #   Experiment
        #   Plate
        #   Screen
        #   Experimenter
        #   ExperimenterGroup
        #   Instrument
        #   ROI
        if 'UUID' in metadata:
            uuid = metadata['UUID'].split(':')[-1]
        else:
            from uuid import uuid1

            uuid = str(uuid1())
        creator = OmeXml._attribute(
            metadata, 'Creator', default=f'tifffile.py {__version__}'
        )
        schema = 'http://www.openmicroscopy.org/Schemas/OME/2016-06'
        self._xml = (
            '{declaration}'
            f'<OME xmlns="{schema}" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            f'xsi:schemaLocation="{schema} {schema}/ome.xsd" '
            f'UUID="urn:uuid:{uuid}"{creator}>'
            '{datasets}'
            '{images}'
            '{annotations}'
            '</OME>'
        )

    def addimage(
        self,
        dtype: DTypeLike,
        shape: Sequence[int],
        storedshape: tuple[int, int, int, int, int, int],
        *,
        axes: str | None = None,
        **metadata: Any,
    ) -> None:
        """Add image to OME-XML.

        The OME model can handle up to 9 dimensional images for selected
        axes orders. Refer to the OME-XML specification for details.
        Non-TZCYXS (modulo) dimensions must be after a TZC dimension or
        require an unused TZC dimension.

        Parameters:
            dtype:
                Data type of image array.
            shape:
                Shape of image array.
            storedshape:
                Normalized shape describing how image array is stored in
                TIFF file as (pages, separate_samples, depth, length, width,
                contig_samples).
            axes:
                Character codes for dimensions in `shape`.
                By default, `axes` is determined from the DimensionOrder
                metadata attribute or matched to the `shape` in reverse order
                of TZC(S)YX(S) based on `storedshape`.
                The following codes are supported: 'S' sample, 'X' width,
                'Y' length, 'Z' depth, 'C' channel, 'T' time, 'A' angle,
                'P' phase, 'R' tile, 'H' lifetime, 'E' lambda, 'Q' other.
            **metadata:
                Additional OME-XML attributes or elements to be stored.

                Image/Pixels:
                    Name, Description,
                    DimensionOrder, TypeDescription,
                    PhysicalSizeX, PhysicalSizeXUnit,
                    PhysicalSizeY, PhysicalSizeYUnit,
                    PhysicalSizeZ, PhysicalSizeZUnit,
                    TimeIncrement, TimeIncrementUnit,
                    StructuredAnnotations, BooleanAnnotation, DoubleAnnotation,
                    LongAnnotation, CommentAnnotation, MapAnnotation,
                    Dataset
                Per Plane:
                    DeltaT, DeltaTUnit,
                    ExposureTime, ExposureTimeUnit,
                    PositionX, PositionXUnit,
                    PositionY, PositionYUnit,
                    PositionZ, PositionZUnit.
                Per Channel:
                    Name, AcquisitionMode, Color, ContrastMethod,
                    EmissionWavelength, EmissionWavelengthUnit,
                    ExcitationWavelength, ExcitationWavelengthUnit,
                    Fluor, IlluminationType, NDFilter,
                    PinholeSize, PinholeSizeUnit, PockelCellSetting.

        Raises:
            OmeXmlError: Image format not supported.

        """
        index = len(self.images)
        annotation_refs = []

        # get Image and Pixels metadata
        metadata = metadata.get('OME', metadata)
        metadata = metadata.get('Image', metadata)
        if isinstance(metadata, (list, tuple)):
            # multiple images
            metadata = metadata[index]
        if 'Pixels' in metadata:
            # merge with Image
            import copy

            metadata = copy.deepcopy(metadata)
            if 'ID' in metadata['Pixels']:
                del metadata['Pixels']['ID']
            metadata.update(metadata['Pixels'])
            del metadata['Pixels']

        try:
            dtype = numpy.dtype(dtype).name
            dtype = {
                'int8': 'int8',
                'int16': 'int16',
                'int32': 'int32',
                'uint8': 'uint8',
                'uint16': 'uint16',
                'uint32': 'uint32',
                'float32': 'float',
                'float64': 'double',
                'complex64': 'complex',
                'complex128': 'double-complex',
                'bool': 'bit',
            }[dtype]
        except KeyError as exc:
            raise OmeXmlError(f'data type {dtype!r} not supported') from exc

        if metadata.get('Type', dtype) != dtype:
            raise OmeXmlError(
                f'metadata Pixels Type {metadata["Type"]!r} '
                f'does not match array dtype {dtype!r}'
            )

        samples = 1
        planecount, separate, depth, length, width, contig = storedshape
        if depth != 1:
            raise OmeXmlError('ImageDepth not supported')
        if not (separate == 1 or contig == 1):
            raise ValueError('invalid stored shape')

        shape = tuple(int(i) for i in shape)
        ndim = len(shape)
        if ndim < 1 or product(shape) <= 0:
            raise OmeXmlError('empty arrays not supported')

        if axes is None:
            # get axes from shape, stored shape, and DimensionOrder
            if contig != 1 or shape[-3:] == (length, width, 1):
                axes = 'YXS'
                samples = contig
            elif separate != 1 or (
                ndim == 6 and shape[-3:] == (1, length, width)
            ):
                axes = 'SYX'
                samples = separate
            else:
                axes = 'YX'
            if not len(axes) <= ndim <= (6 if 'S' in axes else 5):
                raise OmeXmlError(f'{ndim} dimensions not supported')
            hiaxes: str = metadata.get('DimensionOrder', 'XYCZT')[:1:-1]
            axes = hiaxes[(6 if 'S' in axes else 5) - ndim :] + axes
            assert len(axes) == len(shape)

        else:
            # validate axes against shape and stored shape
            axes = axes.upper()
            if len(axes) != len(shape):
                raise ValueError('axes do not match shape')
            if not (
                axes.endswith('YX')
                or axes.endswith('YXS')
                or (axes.endswith('YXC') and 'S' not in axes)
            ):
                raise OmeXmlError('dimensions must end with YX or YXS')
            unique = []
            for ax in axes:
                if ax not in 'TZCYXSAPRHEQ':
                    raise OmeXmlError(f'dimension {ax!r} not supported')
                if ax in unique:
                    raise OmeXmlError(f'multiple {ax!r} dimensions')
                unique.append(ax)
            if ndim > (9 if 'S' in axes else 8):
                raise OmeXmlError('more than 8 dimensions not supported')
            if contig != 1:
                samples = contig
                if ndim < 3:
                    raise ValueError('dimensions do not match stored shape')
                if axes[-1] == 'C':
                    # allow C axis instead of S
                    if 'S' in axes:
                        raise ValueError('invalid axes')
                    axes = axes.replace('C', 'S')
                elif axes[-1] != 'S':
                    raise ValueError('axes do not match stored shape')
                if shape[-1] != contig or shape[-2] != width:
                    raise ValueError('shape does not match stored shape')
            elif separate != 1:
                samples = separate
                if ndim < 3:
                    raise ValueError('dimensions do not match stored shape')
                if axes[-3] == 'C':
                    # allow C axis instead of S
                    if 'S' in axes:
                        raise ValueError('invalid axes')
                    axes = axes.replace('C', 'S')
                elif axes[-3] != 'S':
                    raise ValueError('axes do not match stored shape')
                if shape[-3] != separate or shape[-1] != width:
                    raise ValueError('shape does not match stored shape')

        if shape[axes.index('X')] != width or shape[axes.index('Y')] != length:
            raise ValueError('shape does not match stored shape')

        if 'S' in axes:
            hiaxes = axes[: min(axes.index('S'), axes.index('Y'))]
        else:
            hiaxes = axes[: axes.index('Y')]

        if any(ax in 'APRHEQ' for ax in hiaxes):
            # modulo axes
            modulo = {}
            dimorder = ''
            axestype = {
                'A': 'angle',
                'P': 'phase',
                'R': 'tile',
                'H': 'lifetime',
                'E': 'lambda',
                'Q': 'other',
            }
            axestypedescr = metadata.get('TypeDescription', {})
            for i, ax in enumerate(hiaxes):
                if ax in 'APRHEQ':
                    if ax in axestypedescr:
                        typedescr = f'TypeDescription="{axestypedescr[ax]}" '
                    else:
                        typedescr = ''
                    x = hiaxes[i - 1 : i]
                    if x and x in 'TZC':
                        # use previous axis
                        modulo[x] = axestype[ax], shape[i], typedescr
                    else:
                        # use next unused axis
                        for x in 'TZC':
                            if (
                                x not in dimorder
                                and x not in hiaxes
                                and x not in modulo
                            ):
                                modulo[x] = axestype[ax], shape[i], typedescr
                                dimorder += x
                                break
                        else:
                            # TODO: support any order of axes, such as, APRTZC
                            raise OmeXmlError('more than 3 modulo dimensions')
                else:
                    dimorder += ax
            hiaxes = dimorder

            # TODO: use user-specified start, stop, step, or labels
            moduloalong = ''.join(
                f'<ModuloAlong{ax} Type="{axtype}" {typedescr}'
                f'Start="0" End="{size - 1}"/>'
                for ax, (axtype, size, typedescr) in modulo.items()
            )
            annotation_refs.append(
                f'<AnnotationRef ID="Annotation:{len(self.annotations)}"/>'
            )
            self.annotations.append(
                f'<XMLAnnotation ID="Annotation:{len(self.annotations)}" '
                'Namespace="openmicroscopy.org/omero/dimension/modulo">'
                '<Value>'
                '<Modulo namespace='
                '"http://www.openmicroscopy.org/Schemas/Additions/2011-09">'
                f'{moduloalong}'
                '</Modulo>'
                '</Value>'
                '</XMLAnnotation>'
            )
        else:
            modulo = {}
            annotationref = ''

        hiaxes = hiaxes[::-1]
        for dimorder in (
            metadata.get('DimensionOrder', 'XYCZT'),
            'XYCZT',
            'XYZCT',
            'XYZTC',
            'XYCTZ',
            'XYTCZ',
            'XYTZC',
        ):
            if hiaxes in dimorder:
                break
        else:
            raise OmeXmlError(
                f'dimension order {axes!r} not supported ({hiaxes=})'
            )

        dimsizes = []
        for ax in dimorder:
            if ax == 'S':
                continue
            if ax in axes:
                size = shape[axes.index(ax)]
            else:
                size = 1
            if ax == 'C':
                sizec = size
                size *= samples
            if ax in modulo:
                size *= modulo[ax][1]
            dimsizes.append(size)
        sizes = ''.join(
            f' Size{ax}="{size}"' for ax, size in zip(dimorder, dimsizes)
        )

        # verify DimensionOrder in metadata is compatible
        if 'DimensionOrder' in metadata:
            omedimorder = metadata['DimensionOrder']
            omedimorder = ''.join(
                ax for ax in omedimorder if dimsizes[dimorder.index(ax)] > 1
            )
            if hiaxes not in omedimorder:
                raise OmeXmlError(
                    f'metadata DimensionOrder does not match {axes!r}'
                )

        # verify metadata Size values match shape
        for ax, size in zip(dimorder, dimsizes):
            if metadata.get(f'Size{ax}', size) != size:
                raise OmeXmlError(
                    f'metadata Size{ax} does not match {shape!r}'
                )

        dimsizes[dimorder.index('C')] //= samples
        if planecount != product(dimsizes[2:]):
            raise ValueError('shape does not match stored shape')

        plane_list = []
        planeattributes = metadata.get('Plane', '')
        if planeattributes:
            cztorder = tuple(dimorder[2:].index(ax) for ax in 'CZT')
            for p in range(planecount):
                attributes = OmeXml._attributes(
                    planeattributes,
                    p,
                    'DeltaT',
                    'DeltaTUnit',
                    'ExposureTime',
                    'ExposureTimeUnit',
                    'PositionX',
                    'PositionXUnit',
                    'PositionY',
                    'PositionYUnit',
                    'PositionZ',
                    'PositionZUnit',
                )
                unraveled = numpy.unravel_index(p, dimsizes[2:], order='F')
                c, z, t = (int(unraveled[i]) for i in cztorder)
                plane_list.append(
                    f'<Plane TheC="{c}" TheZ="{z}" TheT="{t}"{attributes}/>'
                )
                # TODO: if possible, verify c, z, t match planeattributes
        planes = ''.join(plane_list)

        channel_list = []
        for c in range(sizec):
            lightpath = '<LightPath/>'
            # TODO: use LightPath elements from metadata
            #    'AnnotationRef',
            #    'DichroicRef',
            #    'EmissionFilterRef',
            #    'ExcitationFilterRef'
            attributes = OmeXml._attributes(
                metadata.get('Channel', ''),
                c,
                'Name',
                'AcquisitionMode',
                'Color',
                'ContrastMethod',
                'EmissionWavelength',
                'EmissionWavelengthUnit',
                'ExcitationWavelength',
                'ExcitationWavelengthUnit',
                'Fluor',
                'IlluminationType',
                'NDFilter',
                'PinholeSize',
                'PinholeSizeUnit',
                'PockelCellSetting',
            )
            channel_list.append(
                f'<Channel ID="Channel:{index}:{c}" '
                f'SamplesPerPixel="{samples}"'
                f'{attributes}>'
                f'{lightpath}'
                '</Channel>'
            )
        channels = ''.join(channel_list)

        # TODO: support more Image elements
        elements = OmeXml._elements(metadata, 'AcquisitionDate', 'Description')

        name = OmeXml._attribute(metadata, 'Name', default=f'Image{index}')
        attributes = OmeXml._attributes(
            metadata,
            None,
            'SignificantBits',
            'PhysicalSizeX',
            'PhysicalSizeXUnit',
            'PhysicalSizeY',
            'PhysicalSizeYUnit',
            'PhysicalSizeZ',
            'PhysicalSizeZUnit',
            'TimeIncrement',
            'TimeIncrementUnit',
        )
        if separate > 1 or contig > 1:
            interleaved = 'false' if separate > 1 else 'true'
            interleaved = f' Interleaved="{interleaved}"'
        else:
            interleaved = ''

        self._dataset(
            metadata.get('Dataset', {}), f'<ImageRef ID="Image:{index}"/>'
        )

        self._annotations(
            metadata.get('StructuredAnnotations', metadata), annotation_refs
        )
        annotationref = ''.join(annotation_refs)

        self.images.append(
            f'<Image ID="Image:{index}"{name}>'
            f'{elements}'
            f'<Pixels ID="Pixels:{index}" '
            f'DimensionOrder="{dimorder}" '
            f'Type="{dtype}"'
            f'{sizes}'
            f'{interleaved}'
            f'{attributes}>'
            f'{channels}'
            f'<TiffData IFD="{self._ifd}" PlaneCount="{planecount}"/>'
            f'{planes}'
            '</Pixels>'
            f'{annotationref}'
            '</Image>'
        )
        self._ifd += planecount

    def tostring(self, *, declaration: bool = False) -> str:
        """Return OME-XML string.

        Parameters:
            declaration: Include XML declaration.

        """
        # TODO: support other top-level elements
        datasets = ''.join(self.datasets)
        images = ''.join(self.images)
        annotations = ''.join(self.annotations)
        if annotations:
            annotations = (
                f'<StructuredAnnotations>{annotations}</StructuredAnnotations>'
            )
        if declaration:
            declaration_str = '<?xml version="1.0" encoding="UTF-8"?>'
        else:
            declaration_str = ''
        xml = self._xml.format(
            declaration=declaration_str,
            images=images,
            annotations=annotations,
            datasets=datasets,
        )
        return xml

    def __repr__(self) -> str:
        return f'<tifffile.OmeXml @0x{id(self):016X}>'

    def __str__(self) -> str:
        """Return OME-XML string."""
        xml = self.tostring()
        try:
            from lxml import etree

            parser = etree.XMLParser(remove_blank_text=True)
            tree = etree.fromstring(xml, parser)
            xml = etree.tostring(
                tree, encoding='utf-8', pretty_print=True, xml_declaration=True
            ).decode()
        except ImportError:
            pass
        except Exception as exc:
            warnings.warn(
                f'<tifffile.OmeXml.__str__> {exc.__class__.__name__}: {exc}',
                UserWarning,
            )
        return xml

    @staticmethod
    def _escape(value: object, /) -> str:
        """Return escaped string of value."""
        if not isinstance(value, str):
            value = str(value)
        elif '&amp;' in value or '&gt;' in value or '&lt;' in value:
            return value
        value = value.replace('&', '&amp;')
        value = value.replace('>', '&gt;')
        value = value.replace('<', '&lt;')
        return value

    @staticmethod
    def _element(
        metadata: dict[str, Any], name: str, default: str | None = None
    ) -> str:
        """Return XML formatted element if name in metadata."""
        value = metadata.get(name, default)
        if value is None:
            return ''
        return f'<{name}>{OmeXml._escape(value)}</{name}>'

    @staticmethod
    def _elements(metadata: dict[str, Any], /, *names: str) -> str:
        """Return XML formatted elements."""
        if not metadata:
            return ''
        elements = (OmeXml._element(metadata, name) for name in names)
        return ''.join(e for e in elements if e)

    @staticmethod
    def _attribute(
        metadata: dict[str, Any],
        name: str,
        /,
        index: int | None = None,
        default: Any = None,
    ) -> str:
        """Return XML formatted attribute if name in metadata."""
        value = metadata.get(name, default)
        if value is None:
            return ''
        if index is not None:
            if isinstance(value, (list, tuple)):
                try:
                    value = value[index]
                except IndexError as exc:
                    raise IndexError(
                        f'list index out of range for attribute {name!r}'
                    ) from exc
            elif index > 0:
                raise TypeError(
                    f'{type(value).__name__!r} is not a list or tuple'
                )
        return f' {name}="{OmeXml._escape(value)}"'

    @staticmethod
    def _attributes(
        metadata: dict[str, Any],
        index_: int | None,
        /,
        *names: str,
    ) -> str:
        """Return XML formatted attributes."""
        if not metadata:
            return ''
        if index_ is None:
            attributes = (OmeXml._attribute(metadata, name) for name in names)
        elif isinstance(metadata, (list, tuple)):
            metadata = metadata[index_]
            attributes = (OmeXml._attribute(metadata, name) for name in names)
        elif isinstance(metadata, dict):
            attributes = (
                OmeXml._attribute(metadata, name, index_) for name in names
            )
        return ''.join(a for a in attributes if a)

    def _dataset(self, metadata: dict[str, Any] | None, imageref: str) -> None:
        """Add Dataset element to self.datasets."""
        index = len(self.datasets)
        if metadata is None:
            # dataset explicitly disabled
            return None
        if not metadata and index == 0:
            # no dataset provided yet
            return None
        if not metadata:
            # use previous dataset
            index -= 1
            if '<AnnotationRef' in self.datasets[index]:
                self.datasets[index] = self.datasets[index].replace(
                    '<AnnotationRef', f'{imageref}<AnnotationRef'
                )
            else:
                self.datasets[index] = self.datasets[index].replace(
                    '</Dataset>', f'{imageref}</Dataset>'
                )
            return None

        # new dataset
        name = metadata.get('Name', '')
        if name:
            name = f' Name="{OmeXml._escape(name)}"'

        description = metadata.get('Description', '')
        if description:
            description = (
                f'<Description>{OmeXml._escape(description)}</Description>'
            )

        annotation_refs: list[str] = []
        self._annotations(metadata, annotation_refs)
        annotationref = ''.join(annotation_refs)

        self.datasets.append(
            f'<Dataset ID="Dataset:{index}"{name}>'
            f'{description}'
            f'{imageref}'
            f'{annotationref}'
            '</Dataset>'
        )
        return None  # f'<DatasetRef ID="Dataset:{index}"/>'

    def _annotations(
        self, metadata: dict[str, Any], annotation_refs: list[str]
    ) -> None:
        """Add annotations to self.annotations and annotation_refs."""
        values: Any
        for name, values in metadata.items():
            if name not in {
                'BooleanAnnotation',
                'DoubleAnnotation',
                'LongAnnotation',
                'CommentAnnotation',
                'MapAnnotation',
                # 'FileAnnotation',
                # 'ListAnnotation',
                # 'TimestampAnnotation,
                # 'XmlAnnotation',
            }:
                continue
            if not values:
                continue
            if not isinstance(values, (list, tuple)):
                values = [values]
            for value in values:
                namespace = ''
                description = ''
                if isinstance(value, dict):
                    value = value.copy()
                    description = value.pop('Description', '')
                    if description:
                        description = (
                            '<Description>'
                            f'{OmeXml._escape(description)}'
                            '</Description>'
                        )
                    namespace = value.pop('Namespace', '')
                    if namespace:
                        namespace = f' Namespace="{OmeXml._escape(namespace)}"'
                    value = value.pop('Value', value)
                if name == 'MapAnnotation':
                    if not isinstance(value, dict):
                        raise ValueError('MapAnnotation is not a dict')
                    values = [
                        f'<M K="{OmeXml._escape(k)}">{OmeXml._escape(v)}</M>'
                        for k, v in value.items()
                    ]
                elif name == 'BooleanAnnotation':
                    values = [f'{bool(value)}'.lower()]
                else:
                    values = [OmeXml._escape(str(value))]
                annotation_refs.append(
                    f'<AnnotationRef ID="Annotation:{len(self.annotations)}"/>'
                )
                self.annotations.append(
                    ''.join(
                        (
                            f'<{name} '
                            f'ID="Annotation:{len(self.annotations)}"'
                            f'{namespace}>',
                            description,
                            '<Value>',
                            ''.join(values),
                            '</Value>',
                            f'</{name}>',
                        )
                    )
                )

    @staticmethod
    def validate(
        omexml: str,
        /,
        omexsd: bytes | None = None,
        assert_: bool = True,
        *,
        _schema: list[Any] = [],  # etree.XMLSchema
    ) -> bool | None:
        r"""Return if OME-XML is valid according to XMLSchema.

        Parameters:
            omexml:
                OME-XML string to validate.
            omexsd:
                Content of OME-XSD schema to validate against.
                By default, the 2016-06 OME XMLSchema is downloaded on first
                run.
            assert\_:
                Raise AssertionError if validation fails.
            _schema:
                Internal use.

        Raises:
            AssertionError:
                Validation failed and `assert\_` is *True*.

        """
        from lxml import etree

        if not _schema:
            if omexsd is None:
                omexsd_path = os.path.join(
                    os.path.dirname(__file__), 'ome.xsd'
                )
                if os.path.exists(omexsd_path):
                    with open(omexsd_path, 'rb') as fh:
                        omexsd = fh.read()
                else:
                    import urllib.request

                    with urllib.request.urlopen(
                        'https://www.openmicroscopy.org/'
                        'Schemas/OME/2016-06/ome.xsd'
                    ) as fh:
                        omexsd = fh.read()
            if omexsd.startswith(b'<?xml'):
                omexsd = omexsd.split(b'>', 1)[-1]
            try:
                _schema.append(
                    etree.XMLSchema(etree.fromstring(omexsd.decode()))
                )
            except Exception:
                # raise
                _schema.append(None)
        if _schema and _schema[0] is not None:
            if omexml.startswith('<?xml'):
                omexml = omexml.split('>', 1)[-1]
            tree = etree.fromstring(omexml)
            if assert_:
                _schema[0].assert_(tree)
                return True
            return bool(_schema[0].validate(tree))
        return None


@final
class CompressionCodec(Mapping[int, Callable[..., object]]):
    """Map :py:class:`COMPRESSION` value to encode or decode function.

    Parameters:
        encode: If *True*, return encode functions, else decode functions.

    """

    _codecs: dict[int, Callable[..., Any]]
    _encode: bool

    def __init__(self, encode: bool) -> None:
        self._codecs = {1: identityfunc}
        self._encode = bool(encode)

    def __getitem__(self, key: int, /) -> Callable[..., Any]:
        if key in self._codecs:
            return self._codecs[key]
        codec: Callable[..., Any]
        try:
            # TODO: enable CCITTRLE decoder for future imagecodecs
            # if key == 2:
            #     if self._encode:
            #         codec = imagecodecs.ccittrle_encode
            #     else:
            #         codec = imagecodecs.ccittrle_decode
            if key == 5:
                if self._encode:
                    codec = imagecodecs.lzw_encode
                else:
                    codec = imagecodecs.lzw_decode
            elif key in {6, 7, 33007}:
                if self._encode:
                    if key in {6, 33007}:
                        raise NotImplementedError
                    codec = imagecodecs.jpeg_encode
                else:
                    codec = imagecodecs.jpeg_decode
            elif key in {8, 32946, 50013}:
                if (
                    hasattr(imagecodecs, 'DEFLATE')
                    and imagecodecs.DEFLATE.available
                ):
                    # imagecodecs built with deflate
                    if self._encode:
                        codec = imagecodecs.deflate_encode
                    else:
                        codec = imagecodecs.deflate_decode
                elif (
                    hasattr(imagecodecs, 'ZLIB') and imagecodecs.ZLIB.available
                ):
                    if self._encode:
                        codec = imagecodecs.zlib_encode
                    else:
                        codec = imagecodecs.zlib_decode
                else:
                    # imagecodecs built without zlib
                    try:
                        from . import _imagecodecs
                    except ImportError:
                        import _imagecodecs  # type: ignore[no-redef]

                    if self._encode:
                        codec = _imagecodecs.zlib_encode
                    else:
                        codec = _imagecodecs.zlib_decode
            elif key == 32773:
                if self._encode:
                    codec = imagecodecs.packbits_encode
                else:
                    codec = imagecodecs.packbits_decode
            elif key in {33003, 33004, 33005, 34712}:
                if self._encode:
                    codec = imagecodecs.jpeg2k_encode
                else:
                    codec = imagecodecs.jpeg2k_decode
            elif key == 34887:
                if self._encode:
                    codec = imagecodecs.lerc_encode
                else:
                    codec = imagecodecs.lerc_decode
            elif key == 34892:
                # DNG lossy
                if self._encode:
                    codec = imagecodecs.jpeg8_encode
                else:
                    codec = imagecodecs.jpeg8_decode
            elif key == 34925:
                if hasattr(imagecodecs, 'LZMA') and imagecodecs.LZMA.available:
                    if self._encode:
                        codec = imagecodecs.lzma_encode
                    else:
                        codec = imagecodecs.lzma_decode
                else:
                    # imagecodecs built without lzma
                    try:
                        from . import _imagecodecs
                    except ImportError:
                        import _imagecodecs  # type: ignore[no-redef]

                    if self._encode:
                        codec = _imagecodecs.lzma_encode
                    else:
                        codec = _imagecodecs.lzma_decode
            elif key == 34933:
                if self._encode:
                    codec = imagecodecs.png_encode
                else:
                    codec = imagecodecs.png_decode
            elif key in {34934, 22610}:
                if self._encode:
                    codec = imagecodecs.jpegxr_encode
                else:
                    codec = imagecodecs.jpegxr_decode
            elif key == 48124:
                if self._encode:
                    codec = imagecodecs.jetraw_encode
                else:
                    codec = imagecodecs.jetraw_decode
            elif key in {50000, 34926}:  # 34926 deprecated
                if hasattr(imagecodecs, 'ZSTD') and imagecodecs.ZSTD.available:
                    if self._encode:
                        codec = imagecodecs.zstd_encode
                    else:
                        codec = imagecodecs.zstd_decode
                else:
                    # imagecodecs built without zstd
                    try:
                        from . import _imagecodecs
                    except ImportError:
                        import _imagecodecs  # type: ignore[no-redef]

                    if self._encode:
                        codec = _imagecodecs.zstd_encode
                    else:
                        codec = _imagecodecs.zstd_decode
            elif key in {50001, 34927}:  # 34927 deprecated
                if self._encode:
                    codec = imagecodecs.webp_encode
                else:
                    codec = imagecodecs.webp_decode
            elif key in {65000, 65001, 65002} and not self._encode:
                codec = imagecodecs.eer_decode
            elif key in {50002, 52546}:
                if self._encode:
                    codec = imagecodecs.jpegxl_encode
                else:
                    codec = imagecodecs.jpegxl_decode
            else:
                try:
                    msg = f'{COMPRESSION(key)!r} not supported'
                except ValueError:
                    msg = f'{key} is not a known COMPRESSION'
                raise KeyError(msg)
        except (AttributeError, ImportError) as exc:
            raise KeyError(
                f'{COMPRESSION(key)!r} ' "requires the 'imagecodecs' package"
            ) from exc
        except NotImplementedError as exc:
            raise KeyError(f'{COMPRESSION(key)!r} not implemented') from exc
        self._codecs[key] = codec
        return codec

    def __contains__(self, key: Any, /) -> bool:
        try:
            self[key]
        except KeyError:
            return False
        return True

    def __iter__(self) -> Iterator[int]:
        yield 1  # dummy

    def __len__(self) -> int:
        return 1  # dummy


@final
class PredictorCodec(Mapping[int, Callable[..., object]]):
    """Map :py:class:`PREDICTOR` value to encode or decode function.

    Parameters:
        encode: If *True*, return encode functions, else decode functions.

    """

    _codecs: dict[int, Callable[..., Any]]
    _encode: bool

    def __init__(self, encode: bool) -> None:
        self._codecs = {1: identityfunc}
        self._encode = bool(encode)

    def __getitem__(self, key: int, /) -> Callable[..., Any]:
        if key in self._codecs:
            return self._codecs[key]
        codec: Callable[..., Any]
        try:
            if key == 2:
                if self._encode:
                    codec = imagecodecs.delta_encode
                else:
                    codec = imagecodecs.delta_decode
            elif key == 3:
                if self._encode:
                    codec = imagecodecs.floatpred_encode
                else:
                    codec = imagecodecs.floatpred_decode
            elif key == 34892:
                if self._encode:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.delta_encode(
                            data, axis=axis, out=out, dist=2
                        )

                else:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.delta_decode(
                            data, axis=axis, out=out, dist=2
                        )

            elif key == 34893:
                if self._encode:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.delta_encode(
                            data, axis=axis, out=out, dist=4
                        )

                else:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.delta_decode(
                            data, axis=axis, out=out, dist=4
                        )

            elif key == 34894:
                if self._encode:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.floatpred_encode(
                            data, axis=axis, out=out, dist=2
                        )

                else:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.floatpred_decode(
                            data, axis=axis, out=out, dist=2
                        )

            elif key == 34895:
                if self._encode:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.floatpred_encode(
                            data, axis=axis, out=out, dist=4
                        )

                else:

                    def codec(data, axis=-1, out=None):
                        return imagecodecs.floatpred_decode(
                            data, axis=axis, out=out, dist=4
                        )

            else:
                raise KeyError(f'{key} is not a known PREDICTOR')
        except AttributeError as exc:
            raise KeyError(
                f'{PREDICTOR(key)!r}' " requires the 'imagecodecs' package"
            ) from exc
        except NotImplementedError as exc:
            raise KeyError(f'{PREDICTOR(key)!r} not implemented') from exc
        self._codecs[key] = codec
        return codec

    def __contains__(self, key: Any, /) -> bool:
        try:
            self[key]
        except KeyError:
            return False
        return True

    def __iter__(self) -> Iterator[int]:
        yield 1  # dummy

    def __len__(self) -> int:
        return 1  # dummy


class DATATYPE(enum.IntEnum):
    """TIFF tag data types."""

    BYTE = 1
    """8-bit unsigned integer."""
    ASCII = 2
    """8-bit byte with last byte null, containing 7-bit ASCII code."""
    SHORT = 3
    """16-bit unsigned integer."""
    LONG = 4
    """32-bit unsigned integer."""
    RATIONAL = 5
    """Two 32-bit unsigned integers, numerator and denominator of fraction."""
    SBYTE = 6
    """8-bit signed integer."""
    UNDEFINED = 7
    """8-bit byte that may contain anything."""
    SSHORT = 8
    """16-bit signed integer."""
    SLONG = 9
    """32-bit signed integer."""
    SRATIONAL = 10
    """Two 32-bit signed integers, numerator and denominator of fraction."""
    FLOAT = 11
    """Single precision (4-byte) IEEE format."""
    DOUBLE = 12
    """Double precision (8-byte) IEEE format."""
    IFD = 13
    """Unsigned 4 byte IFD offset."""
    UNICODE = 14
    """UTF-16 (2-byte) unicode string."""
    COMPLEX = 15
    """Single precision (8-byte) complex number."""
    LONG8 = 16
    """Unsigned 8 byte integer (BigTIFF)."""
    SLONG8 = 17
    """Signed 8 byte integer (BigTIFF)."""
    IFD8 = 18
    """Unsigned 8 byte IFD offset (BigTIFF)."""


class COMPRESSION(enum.IntEnum):
    """Values of Compression tag.

    Compression scheme used on image data.

    """

    NONE = 1
    """No compression (default)."""
    CCITTRLE = 2  # CCITT 1D
    CCITT_T4 = 3  # T4/Group 3 Fax
    CCITT_T6 = 4  # T6/Group 4 Fax
    LZW = 5
    """Lempel-Ziv-Welch."""
    OJPEG = 6  # old-style JPEG
    JPEG = 7
    """New style JPEG."""
    ADOBE_DEFLATE = 8
    """Deflate, aka ZLIB."""
    JBIG_BW = 9  # VC5
    JBIG_COLOR = 10
    JPEG_99 = 99  # Leaf MOS lossless JPEG
    IMPACJ = 103  # Pegasus Imaging Corporation DCT
    KODAK_262 = 262
    JPEGXR_NDPI = 22610
    """JPEG XR (Hammatsu NDPI)."""
    NEXT = 32766
    SONY_ARW = 32767
    PACKED_RAW = 32769
    SAMSUNG_SRW = 32770
    CCIRLEW = 32771  # Word-aligned 1D Huffman compression
    SAMSUNG_SRW2 = 32772
    PACKBITS = 32773
    """PackBits, aka Macintosh RLE."""
    THUNDERSCAN = 32809
    IT8CTPAD = 32895  # TIFF/IT
    IT8LW = 32896  # TIFF/IT
    IT8MP = 32897  # TIFF/IT
    IT8BL = 32898  # TIFF/IT
    PIXARFILM = 32908
    PIXARLOG = 32909
    DEFLATE = 32946
    DCS = 32947
    APERIO_JP2000_YCBC = 33003  # Matrox libraries
    """JPEG 2000 YCbCr (Leica Aperio)."""
    JPEG_2000_LOSSY = 33004
    """Lossy JPEG 2000 (Bio-Formats)."""
    APERIO_JP2000_RGB = 33005  # Kakadu libraries
    """JPEG 2000 RGB (Leica Aperio)."""
    ALT_JPEG = 33007
    """JPEG (Bio-Formats)."""
    # PANASONIC_RAW1 = 34316
    # PANASONIC_RAW2 = 34826
    # PANASONIC_RAW3 = 34828
    # PANASONIC_RAW4 = 34830
    JBIG = 34661
    SGILOG = 34676  # LogLuv32
    SGILOG24 = 34677
    LURADOC = 34692  # LuraWave
    JPEG2000 = 34712
    """JPEG 2000."""
    NIKON_NEF = 34713
    JBIG2 = 34715
    MDI_BINARY = 34718  # Microsoft Document Imaging
    MDI_PROGRESSIVE = 34719  # Microsoft Document Imaging
    MDI_VECTOR = 34720  # Microsoft Document Imaging
    LERC = 34887
    """ESRI Limited Error Raster Compression."""
    JPEG_LOSSY = 34892  # DNG
    LZMA = 34925
    """Lempel-Ziv-Markov chain Algorithm."""
    ZSTD_DEPRECATED = 34926
    WEBP_DEPRECATED = 34927
    PNG = 34933  # Objective Pathology Services
    """Portable Network Graphics (Zoomable Image File format)."""
    JPEGXR = 34934
    """JPEG XR (Zoomable Image File format)."""
    JETRAW = 48124
    """Jetraw by Dotphoton."""
    ZSTD = 50000
    """Zstandard."""
    WEBP = 50001
    """WebP."""
    JPEGXL = 50002  # GDAL
    """JPEG XL."""
    PIXTIFF = 50013
    """ZLIB (Atalasoft)."""
    JPEGXL_DNG = 52546
    """JPEG XL (DNG)."""
    EER_V0 = 65000  # FIXED82 Thermo Fisher Scientific
    EER_V1 = 65001  # FIXED72 Thermo Fisher Scientific
    EER_V2 = 65002  # VARIABLE Thermo Fisher Scientific
    # KODAK_DCR = 65000
    # PENTAX_PEF = 65535

    def __bool__(self) -> bool:
        return self > 1


class PREDICTOR(enum.IntEnum):
    """Values of Predictor tag.

    A mathematical operator that is applied to the image data before
    compression.

    """

    NONE = 1
    """No prediction scheme used (default)."""
    HORIZONTAL = 2
    """Horizontal differencing."""
    FLOATINGPOINT = 3
    """Floating-point horizontal differencing."""
    HORIZONTALX2 = 34892  # DNG
    HORIZONTALX4 = 34893
    FLOATINGPOINTX2 = 34894
    FLOATINGPOINTX4 = 34895

    def __bool__(self) -> bool:
        return self > 1


class PHOTOMETRIC(enum.IntEnum):
    """Values of PhotometricInterpretation tag.

    The color space of the image.

    """

    MINISWHITE = 0
    """For bilevel and grayscale images, 0 is imaged as white."""
    MINISBLACK = 1
    """For bilevel and grayscale images, 0 is imaged as black."""
    RGB = 2
    """Chroma components are Red, Green, Blue."""
    PALETTE = 3
    """Single chroma component is index into colormap."""
    MASK = 4
    SEPARATED = 5
    """Chroma components are Cyan, Magenta, Yellow, and Key (black)."""
    YCBCR = 6
    """Chroma components are Luma, blue-difference, and red-difference."""
    CIELAB = 8
    ICCLAB = 9
    ITULAB = 10
    CFA = 32803
    """Color Filter Array."""
    LOGL = 32844
    LOGLUV = 32845
    LINEAR_RAW = 34892
    DEPTH_MAP = 51177  # DNG 1.5
    SEMANTIC_MASK = 52527  # DNG 1.6


class FILETYPE(enum.IntFlag):
    """Values of NewSubfileType tag.

    A general indication of the kind of the image.

    """

    UNDEFINED = 0
    """Image is full-resolution (default)."""
    REDUCEDIMAGE = 1
    """Image is reduced-resolution version of another image."""
    PAGE = 2
    """Image is single page of multi-page image."""
    MASK = 4
    """Image is transparency mask for another image."""
    MACRO = 8  # Aperio SVS, or DNG Depth map
    """Image is MACRO image (SVS) or depth map for another image (DNG)."""
    ENHANCED = 16  # DNG
    """Image contains enhanced image (DNG)."""
    DNG = 65536  # 65537: Alternative, 65540: Semantic mask


class OFILETYPE(enum.IntEnum):
    """Values of deprecated SubfileType tag."""

    UNDEFINED = 0
    IMAGE = 1  # full-resolution image
    REDUCEDIMAGE = 2  # reduced-resolution image
    PAGE = 3  # single page of multi-page image


class FILLORDER(enum.IntEnum):
    """Values of FillOrder tag.

    The logical order of bits within a byte.

    """

    MSB2LSB = 1
    """Pixel values are stored in higher-order bits of byte (default)."""
    LSB2MSB = 2
    """Pixels values are stored in lower-order bits of byte."""


class ORIENTATION(enum.IntEnum):
    """Values of Orientation tag.

    The orientation of the image with respect to the rows and columns.

    """

    TOPLEFT = 1  # default
    TOPRIGHT = 2
    BOTRIGHT = 3
    BOTLEFT = 4
    LEFTTOP = 5
    RIGHTTOP = 6
    RIGHTBOT = 7
    LEFTBOT = 8


class PLANARCONFIG(enum.IntEnum):
    """Values of PlanarConfiguration tag.

    Specifies how components of each pixel are stored.

    """

    CONTIG = 1
    """Chunky, component values are stored contiguously (default)."""
    SEPARATE = 2
    """Planar, component values are stored in separate planes."""


class RESUNIT(enum.IntEnum):
    """Values of ResolutionUnit tag.

    The unit of measurement for XResolution and YResolution.

    """

    NONE = 1
    """No absolute unit of measurement."""
    INCH = 2
    """Inch (default)."""
    CENTIMETER = 3
    """Centimeter."""
    MILLIMETER = 4
    """Millimeter (DNG)."""
    MICROMETER = 5
    """Micrometer (DNG)."""

    def __bool__(self) -> bool:
        return self > 1


class EXTRASAMPLE(enum.IntEnum):
    """Values of ExtraSamples tag.

    Interpretation of extra components in a pixel.

    """

    UNSPECIFIED = 0
    """Unspecified data."""
    ASSOCALPHA = 1
    """Associated alpha data with premultiplied color."""
    UNASSALPHA = 2
    """Unassociated alpha data."""


class SAMPLEFORMAT(enum.IntEnum):
    """Values of SampleFormat tag.

    Data type of samples in a pixel.

    """

    UINT = 1
    """Unsigned integer."""
    INT = 2
    """Signed integer."""
    IEEEFP = 3
    """IEEE floating-point"""
    VOID = 4
    """Undefined."""
    COMPLEXINT = 5
    """Complex integer."""
    COMPLEXIEEEFP = 6
    """Complex floating-point."""


class CHUNKMODE(enum.IntEnum):
    """ZarrStore chunk modes.

    Specifies how to chunk data in Zarr 2 stores.

    """

    STRILE = 0
    """Chunk is strip or tile."""
    PLANE = 1
    """Chunk is image plane."""
    PAGE = 2
    """Chunk is image in page."""
    FILE = 3
    """Chunk is image in file."""


# class THRESHOLD(enum.IntEnum):
#     BILEVEL = 1
#     HALFTONE = 2
#     ERRORDIFFUSE = 3
#
# class GRAYRESPONSEUNIT(enum.IntEnum):
#     _10S = 1
#     _100S = 2
#     _1000S = 3
#     _10000S = 4
#     _100000S = 5
#
# class COLORRESPONSEUNIT(enum.IntEnum):
#     _10S = 1
#     _100S = 2
#     _1000S = 3
#     _10000S = 4
#     _100000S = 5
#
# class GROUP4OPT(enum.IntEnum):
#     UNCOMPRESSED = 2


class _TIFF:
    """Delay-loaded constants, accessible via :py:attr:`TIFF` instance."""

    @cached_property
    def CLASSIC_LE(self) -> TiffFormat:
        """32-bit little-endian TIFF format."""
        return TiffFormat(
            version=42,
            byteorder='<',
            offsetsize=4,
            offsetformat='<I',
            tagnosize=2,
            tagnoformat='<H',
            tagsize=12,
            tagformat1='<HH',
            tagformat2='<I4s',
            tagoffsetthreshold=4,
        )

    @cached_property
    def CLASSIC_BE(self) -> TiffFormat:
        """32-bit big-endian TIFF format."""
        return TiffFormat(
            version=42,
            byteorder='>',
            offsetsize=4,
            offsetformat='>I',
            tagnosize=2,
            tagnoformat='>H',
            tagsize=12,
            tagformat1='>HH',
            tagformat2='>I4s',
            tagoffsetthreshold=4,
        )

    @cached_property
    def BIG_LE(self) -> TiffFormat:
        """64-bit little-endian TIFF format."""
        return TiffFormat(
            version=43,
            byteorder='<',
            offsetsize=8,
            offsetformat='<Q',
            tagnosize=8,
            tagnoformat='<Q',
            tagsize=20,
            tagformat1='<HH',
            tagformat2='<Q8s',
            tagoffsetthreshold=8,
        )

    @cached_property
    def BIG_BE(self) -> TiffFormat:
        """64-bit big-endian TIFF format."""
        return TiffFormat(
            version=43,
            byteorder='>',
            offsetsize=8,
            offsetformat='>Q',
            tagnosize=8,
            tagnoformat='>Q',
            tagsize=20,
            tagformat1='>HH',
            tagformat2='>Q8s',
            tagoffsetthreshold=8,
        )

    @cached_property
    def NDPI_LE(self) -> TiffFormat:
        """32-bit little-endian TIFF format with 64-bit offsets."""
        return TiffFormat(
            version=42,
            byteorder='<',
            offsetsize=8,  # NDPI uses 8 bytes IFD and tag offsets
            offsetformat='<Q',
            tagnosize=2,
            tagnoformat='<H',
            tagsize=12,  # 16 after patching
            tagformat1='<HH',
            tagformat2='<I8s',  # after patching
            tagoffsetthreshold=4,
        )

    @cached_property
    def TAGS(self) -> TiffTagRegistry:
        """Registry of TIFF tag codes and names from TIFF6, TIFF/EP, EXIF."""
        # TODO: divide into baseline, exif, private, ... tags
        return TiffTagRegistry(
            (
                (11, 'ProcessingSoftware'),
                (254, 'NewSubfileType'),
                (255, 'SubfileType'),
                (256, 'ImageWidth'),
                (257, 'ImageLength'),
                (258, 'BitsPerSample'),
                (259, 'Compression'),
                (262, 'PhotometricInterpretation'),
                (263, 'Thresholding'),
                (264, 'CellWidth'),
                (265, 'CellLength'),
                (266, 'FillOrder'),
                (269, 'DocumentName'),
                (270, 'ImageDescription'),
                (271, 'Make'),
                (272, 'Model'),
                (273, 'StripOffsets'),
                (274, 'Orientation'),
                (277, 'SamplesPerPixel'),
                (278, 'RowsPerStrip'),
                (279, 'StripByteCounts'),
                (280, 'MinSampleValue'),
                (281, 'MaxSampleValue'),
                (282, 'XResolution'),
                (283, 'YResolution'),
                (284, 'PlanarConfiguration'),
                (285, 'PageName'),
                (286, 'XPosition'),
                (287, 'YPosition'),
                (288, 'FreeOffsets'),
                (289, 'FreeByteCounts'),
                (290, 'GrayResponseUnit'),
                (291, 'GrayResponseCurve'),
                (292, 'T4Options'),
                (293, 'T6Options'),
                (296, 'ResolutionUnit'),
                (297, 'PageNumber'),
                (300, 'ColorResponseUnit'),
                (301, 'TransferFunction'),
                (305, 'Software'),
                (306, 'DateTime'),
                (315, 'Artist'),
                (316, 'HostComputer'),
                (317, 'Predictor'),
                (318, 'WhitePoint'),
                (319, 'PrimaryChromaticities'),
                (320, 'ColorMap'),
                (321, 'HalftoneHints'),
                (322, 'TileWidth'),
                (323, 'TileLength'),
                (324, 'TileOffsets'),
                (325, 'TileByteCounts'),
                (326, 'BadFaxLines'),
                (327, 'CleanFaxData'),
                (328, 'ConsecutiveBadFaxLines'),
                (330, 'SubIFDs'),
                (332, 'InkSet'),
                (333, 'InkNames'),
                (334, 'NumberOfInks'),
                (336, 'DotRange'),
                (337, 'TargetPrinter'),
                (338, 'ExtraSamples'),
                (339, 'SampleFormat'),
                (340, 'SMinSampleValue'),
                (341, 'SMaxSampleValue'),
                (342, 'TransferRange'),
                (343, 'ClipPath'),
                (344, 'XClipPathUnits'),
                (345, 'YClipPathUnits'),
                (346, 'Indexed'),
                (347, 'JPEGTables'),
                (351, 'OPIProxy'),
                (400, 'GlobalParametersIFD'),
                (401, 'ProfileType'),
                (402, 'FaxProfile'),
                (403, 'CodingMethods'),
                (404, 'VersionYear'),
                (405, 'ModeNumber'),
                (433, 'Decode'),
                (434, 'DefaultImageColor'),
                (435, 'T82Options'),
                (437, 'JPEGTables'),  # 347
                (512, 'JPEGProc'),
                (513, 'JPEGInterchangeFormat'),
                (514, 'JPEGInterchangeFormatLength'),
                (515, 'JPEGRestartInterval'),
                (517, 'JPEGLosslessPredictors'),
                (518, 'JPEGPointTransforms'),
                (519, 'JPEGQTables'),
                (520, 'JPEGDCTables'),
                (521, 'JPEGACTables'),
                (529, 'YCbCrCoefficients'),
                (530, 'YCbCrSubSampling'),
                (531, 'YCbCrPositioning'),
                (532, 'ReferenceBlackWhite'),
                (559, 'StripRowCounts'),
                (700, 'XMP'),  # XMLPacket
                (769, 'GDIGamma'),  # GDI+
                (770, 'ICCProfileDescriptor'),  # GDI+
                (771, 'SRGBRenderingIntent'),  # GDI+
                (800, 'ImageTitle'),  # GDI+
                (907, 'SiffCompress'),  # https://github.com/MaimonLab/SiffPy
                (999, 'USPTO_Miscellaneous'),
                (4864, 'AndorId'),  # TODO, Andor Technology 4864 - 5030
                (4869, 'AndorTemperature'),
                (4876, 'AndorExposureTime'),
                (4878, 'AndorKineticCycleTime'),
                (4879, 'AndorAccumulations'),
                (4881, 'AndorAcquisitionCycleTime'),
                (4882, 'AndorReadoutTime'),
                (4884, 'AndorPhotonCounting'),
                (4885, 'AndorEmDacLevel'),
                (4890, 'AndorFrames'),
                (4896, 'AndorHorizontalFlip'),
                (4897, 'AndorVerticalFlip'),
                (4898, 'AndorClockwise'),
                (4899, 'AndorCounterClockwise'),
                (4904, 'AndorVerticalClockVoltage'),
                (4905, 'AndorVerticalShiftSpeed'),
                (4907, 'AndorPreAmpSetting'),
                (4908, 'AndorCameraSerial'),
                (4911, 'AndorActualTemperature'),
                (4912, 'AndorBaselineClamp'),
                (4913, 'AndorPrescans'),
                (4914, 'AndorModel'),
                (4915, 'AndorChipSizeX'),
                (4916, 'AndorChipSizeY'),
                (4944, 'AndorBaselineOffset'),
                (4966, 'AndorSoftwareVersion'),
                (18246, 'Rating'),
                (18247, 'XP_DIP_XML'),
                (18248, 'StitchInfo'),
                (18249, 'RatingPercent'),
                (20481, 'ResolutionXUnit'),  # GDI+
                (20482, 'ResolutionYUnit'),  # GDI+
                (20483, 'ResolutionXLengthUnit'),  # GDI+
                (20484, 'ResolutionYLengthUnit'),  # GDI+
                (20485, 'PrintFlags'),  # GDI+
                (20486, 'PrintFlagsVersion'),  # GDI+
                (20487, 'PrintFlagsCrop'),  # GDI+
                (20488, 'PrintFlagsBleedWidth'),  # GDI+
                (20489, 'PrintFlagsBleedWidthScale'),  # GDI+
                (20490, 'HalftoneLPI'),  # GDI+
                (20491, 'HalftoneLPIUnit'),  # GDI+
                (20492, 'HalftoneDegree'),  # GDI+
                (20493, 'HalftoneShape'),  # GDI+
                (20494, 'HalftoneMisc'),  # GDI+
                (20495, 'HalftoneScreen'),  # GDI+
                (20496, 'JPEGQuality'),  # GDI+
                (20497, 'GridSize'),  # GDI+
                (20498, 'ThumbnailFormat'),  # GDI+
                (20499, 'ThumbnailWidth'),  # GDI+
                (20500, 'ThumbnailHeight'),  # GDI+
                (20501, 'ThumbnailColorDepth'),  # GDI+
                (20502, 'ThumbnailPlanes'),  # GDI+
                (20503, 'ThumbnailRawBytes'),  # GDI+
                (20504, 'ThumbnailSize'),  # GDI+
                (20505, 'ThumbnailCompressedSize'),  # GDI+
                (20506, 'ColorTransferFunction'),  # GDI+
                (20507, 'ThumbnailData'),
                (20512, 'ThumbnailImageWidth'),  # GDI+
                (20513, 'ThumbnailImageHeight'),  # GDI+
                (20514, 'ThumbnailBitsPerSample'),  # GDI+
                (20515, 'ThumbnailCompression'),
                (20516, 'ThumbnailPhotometricInterp'),  # GDI+
                (20517, 'ThumbnailImageDescription'),  # GDI+
                (20518, 'ThumbnailEquipMake'),  # GDI+
                (20519, 'ThumbnailEquipModel'),  # GDI+
                (20520, 'ThumbnailStripOffsets'),  # GDI+
                (20521, 'ThumbnailOrientation'),  # GDI+
                (20522, 'ThumbnailSamplesPerPixel'),  # GDI+
                (20523, 'ThumbnailRowsPerStrip'),  # GDI+
                (20524, 'ThumbnailStripBytesCount'),  # GDI+
                (20525, 'ThumbnailResolutionX'),
                (20526, 'ThumbnailResolutionY'),
                (20527, 'ThumbnailPlanarConfig'),  # GDI+
                (20528, 'ThumbnailResolutionUnit'),
                (20529, 'ThumbnailTransferFunction'),
                (20530, 'ThumbnailSoftwareUsed'),  # GDI+
                (20531, 'ThumbnailDateTime'),  # GDI+
                (20532, 'ThumbnailArtist'),  # GDI+
                (20533, 'ThumbnailWhitePoint'),  # GDI+
                (20534, 'ThumbnailPrimaryChromaticities'),  # GDI+
                (20535, 'ThumbnailYCbCrCoefficients'),  # GDI+
                (20536, 'ThumbnailYCbCrSubsampling'),  # GDI+
                (20537, 'ThumbnailYCbCrPositioning'),
                (20538, 'ThumbnailRefBlackWhite'),  # GDI+
                (20539, 'ThumbnailCopyRight'),  # GDI+
                (20545, 'InteroperabilityIndex'),
                (20546, 'InteroperabilityVersion'),
                (20624, 'LuminanceTable'),
                (20625, 'ChrominanceTable'),
                (20736, 'FrameDelay'),  # GDI+
                (20737, 'LoopCount'),  # GDI+
                (20738, 'GlobalPalette'),  # GDI+
                (20739, 'IndexBackground'),  # GDI+
                (20740, 'IndexTransparent'),  # GDI+
                (20752, 'PixelUnit'),  # GDI+
                (20753, 'PixelPerUnitX'),  # GDI+
                (20754, 'PixelPerUnitY'),  # GDI+
                (20755, 'PaletteHistogram'),  # GDI+
                (28672, 'SonyRawFileType'),  # Sony ARW
                (28722, 'VignettingCorrParams'),  # Sony ARW
                (28725, 'ChromaticAberrationCorrParams'),  # Sony ARW
                (28727, 'DistortionCorrParams'),  # Sony ARW
                # Private tags >= 32768
                (32781, 'ImageID'),
                (32931, 'WangTag1'),
                (32932, 'WangAnnotation'),
                (32933, 'WangTag3'),
                (32934, 'WangTag4'),
                (32953, 'ImageReferencePoints'),
                (32954, 'RegionXformTackPoint'),
                (32955, 'WarpQuadrilateral'),
                (32956, 'AffineTransformMat'),
                (32995, 'Matteing'),
                (32996, 'DataType'),  # use SampleFormat
                (32997, 'ImageDepth'),
                (32998, 'TileDepth'),
                (33300, 'ImageFullWidth'),
                (33301, 'ImageFullLength'),
                (33302, 'TextureFormat'),
                (33303, 'TextureWrapModes'),
                (33304, 'FieldOfViewCotangent'),
                (33305, 'MatrixWorldToScreen'),
                (33306, 'MatrixWorldToCamera'),
                (33405, 'Model2'),
                (33421, 'CFARepeatPatternDim'),
                (33422, 'CFAPattern'),
                (33423, 'BatteryLevel'),
                (33424, 'KodakIFD'),
                (33434, 'ExposureTime'),
                (33437, 'FNumber'),
                (33432, 'Copyright'),
                (33445, 'MDFileTag'),
                (33446, 'MDScalePixel'),
                (33447, 'MDColorTable'),
                (33448, 'MDLabName'),
                (33449, 'MDSampleInfo'),
                (33450, 'MDPrepDate'),
                (33451, 'MDPrepTime'),
                (33452, 'MDFileUnits'),
                (33465, 'NiffRotation'),  # NIFF
                (33466, 'NiffNavyCompression'),  # NIFF
                (33467, 'NiffTileIndex'),  # NIFF
                (33471, 'OlympusINI'),
                (33550, 'ModelPixelScaleTag'),
                (33560, 'OlympusSIS'),  # see also 33471 and 34853
                (33589, 'AdventScale'),
                (33590, 'AdventRevision'),
                (33628, 'UIC1tag'),  # Metamorph  Universal Imaging Corp STK
                (33629, 'UIC2tag'),
                (33630, 'UIC3tag'),
                (33631, 'UIC4tag'),
                (33723, 'IPTCNAA'),
                (33858, 'ExtendedTagsOffset'),  # DEFF points IFD with tags
                (33918, 'IntergraphPacketData'),  # INGRPacketDataTag
                (33919, 'IntergraphFlagRegisters'),  # INGRFlagRegisters
                (33920, 'IntergraphMatrixTag'),  # IrasBTransformationMatrix
                (33921, 'INGRReserved'),
                (33922, 'ModelTiepointTag'),
                (33923, 'LeicaMagic'),
                (34016, 'Site'),  # 34016..34032 ANSI IT8 TIFF/IT
                (34017, 'ColorSequence'),
                (34018, 'IT8Header'),
                (34019, 'RasterPadding'),
                (34020, 'BitsPerRunLength'),
                (34021, 'BitsPerExtendedRunLength'),
                (34022, 'ColorTable'),
                (34023, 'ImageColorIndicator'),
                (34024, 'BackgroundColorIndicator'),
                (34025, 'ImageColorValue'),
                (34026, 'BackgroundColorValue'),
                (34027, 'PixelIntensityRange'),
                (34028, 'TransparencyIndicator'),
                (34029, 'ColorCharacterization'),
                (34030, 'HCUsage'),
                (34031, 'TrapIndicator'),
                (34032, 'CMYKEquivalent'),
                (34118, 'CZ_SEM'),  # Zeiss SEM
                (34152, 'AFCP_IPTC'),
                (34232, 'PixelMagicJBIGOptions'),  # EXIF, also TI FrameCount
                (34263, 'JPLCartoIFD'),
                (34122, 'IPLAB'),  # number of images
                (34264, 'ModelTransformationTag'),
                (34306, 'WB_GRGBLevels'),  # Leaf MOS
                (34310, 'LeafData'),
                (34361, 'MM_Header'),
                (34362, 'MM_Stamp'),
                (34363, 'MM_Unknown'),
                (34377, 'ImageResources'),  # Photoshop
                (34386, 'MM_UserBlock'),
                (34412, 'CZ_LSMINFO'),
                (34665, 'ExifTag'),
                (34675, 'InterColorProfile'),  # ICCProfile
                (34680, 'FEI_SFEG'),  #
                (34682, 'FEI_HELIOS'),  #
                (34683, 'FEI_TITAN'),  #
                (34687, 'FXExtensions'),
                (34688, 'MultiProfiles'),
                (34689, 'SharedData'),
                (34690, 'T88Options'),
                (34710, 'MarCCD'),  # offset to MarCCD header
                (34732, 'ImageLayer'),
                (34735, 'GeoKeyDirectoryTag'),
                (34736, 'GeoDoubleParamsTag'),
                (34737, 'GeoAsciiParamsTag'),
                (34750, 'JBIGOptions'),
                (34821, 'PIXTIFF'),  # ? Pixel Translations Inc
                (34850, 'ExposureProgram'),
                (34852, 'SpectralSensitivity'),
                (34853, 'GPSTag'),  # GPSIFD  also OlympusSIS2
                (34853, 'OlympusSIS2'),
                (34855, 'ISOSpeedRatings'),
                (34855, 'PhotographicSensitivity'),
                (34856, 'OECF'),  # optoelectric conversion factor
                (34857, 'Interlace'),  # TIFF/EP
                (34858, 'TimeZoneOffset'),  # TIFF/EP
                (34859, 'SelfTimerMode'),  # TIFF/EP
                (34864, 'SensitivityType'),
                (34865, 'StandardOutputSensitivity'),
                (34866, 'RecommendedExposureIndex'),
                (34867, 'ISOSpeed'),
                (34868, 'ISOSpeedLatitudeyyy'),
                (34869, 'ISOSpeedLatitudezzz'),
                (34908, 'HylaFAXFaxRecvParams'),
                (34909, 'HylaFAXFaxSubAddress'),
                (34910, 'HylaFAXFaxRecvTime'),
                (34911, 'FaxDcs'),
                (34929, 'FedexEDR'),
                (34954, 'LeafSubIFD'),
                (34959, 'Aphelion1'),
                (34960, 'Aphelion2'),
                (34961, 'AphelionInternal'),  # ADCIS
                (36864, 'ExifVersion'),
                (36867, 'DateTimeOriginal'),
                (36868, 'DateTimeDigitized'),
                (36873, 'GooglePlusUploadCode'),
                (36880, 'OffsetTime'),
                (36881, 'OffsetTimeOriginal'),
                (36882, 'OffsetTimeDigitized'),
                # TODO, Pilatus/CHESS/TV6 36864..37120 conflicting with Exif
                (36864, 'TVX_Unknown'),
                (36865, 'TVX_NumExposure'),
                (36866, 'TVX_NumBackground'),
                (36867, 'TVX_ExposureTime'),
                (36868, 'TVX_BackgroundTime'),
                (36870, 'TVX_Unknown'),
                (36873, 'TVX_SubBpp'),
                (36874, 'TVX_SubWide'),
                (36875, 'TVX_SubHigh'),
                (36876, 'TVX_BlackLevel'),
                (36877, 'TVX_DarkCurrent'),
                (36878, 'TVX_ReadNoise'),
                (36879, 'TVX_DarkCurrentNoise'),
                (36880, 'TVX_BeamMonitor'),
                (37120, 'TVX_UserVariables'),  # A/D values
                (37121, 'ComponentsConfiguration'),
                (37122, 'CompressedBitsPerPixel'),
                (37377, 'ShutterSpeedValue'),
                (37378, 'ApertureValue'),
                (37379, 'BrightnessValue'),
                (37380, 'ExposureBiasValue'),
                (37381, 'MaxApertureValue'),
                (37382, 'SubjectDistance'),
                (37383, 'MeteringMode'),
                (37384, 'LightSource'),
                (37385, 'Flash'),
                (37386, 'FocalLength'),
                (37387, 'FlashEnergy'),  # TIFF/EP
                (37388, 'SpatialFrequencyResponse'),  # TIFF/EP
                (37389, 'Noise'),  # TIFF/EP
                (37390, 'FocalPlaneXResolution'),  # TIFF/EP
                (37391, 'FocalPlaneYResolution'),  # TIFF/EP
                (37392, 'FocalPlaneResolutionUnit'),  # TIFF/EP
                (37393, 'ImageNumber'),  # TIFF/EP
                (37394, 'SecurityClassification'),  # TIFF/EP
                (37395, 'ImageHistory'),  # TIFF/EP
                (37396, 'SubjectLocation'),  # TIFF/EP
                (37397, 'ExposureIndex'),  # TIFF/EP
                (37398, 'TIFFEPStandardID'),  # TIFF/EP
                (37399, 'SensingMethod'),  # TIFF/EP
                (37434, 'CIP3DataFile'),
                (37435, 'CIP3Sheet'),
                (37436, 'CIP3Side'),
                (37439, 'StoNits'),
                (37500, 'MakerNote'),
                (37510, 'UserComment'),
                (37520, 'SubsecTime'),
                (37521, 'SubsecTimeOriginal'),
                (37522, 'SubsecTimeDigitized'),
                (37679, 'MODIText'),  # Microsoft Office Document Imaging
                (37680, 'MODIOLEPropertySetStorage'),
                (37681, 'MODIPositioning'),
                (37701, 'AgilentBinary'),  # private structure
                (37702, 'AgilentString'),  # file description
                (37706, 'TVIPS'),  # offset to TemData structure
                (37707, 'TVIPS1'),
                (37708, 'TVIPS2'),  # same TemData structure as undefined
                (37724, 'ImageSourceData'),  # Photoshop
                (37888, 'Temperature'),
                (37889, 'Humidity'),
                (37890, 'Pressure'),
                (37891, 'WaterDepth'),
                (37892, 'Acceleration'),
                (37893, 'CameraElevationAngle'),
                (40000, 'XPos'),  # Janelia
                (40001, 'YPos'),
                (40002, 'ZPos'),
                (40001, 'MC_IpWinScal'),  # Media Cybernetics
                (40001, 'RecipName'),  # MS FAX
                (40002, 'RecipNumber'),
                (40003, 'SenderName'),
                (40004, 'Routing'),
                (40005, 'CallerId'),
                (40006, 'TSID'),
                (40007, 'CSID'),
                (40008, 'FaxTime'),
                (40100, 'MC_IdOld'),
                (40106, 'MC_Unknown'),
                (40965, 'InteroperabilityTag'),  # InteropOffset
                (40091, 'XPTitle'),
                (40092, 'XPComment'),
                (40093, 'XPAuthor'),
                (40094, 'XPKeywords'),
                (40095, 'XPSubject'),
                (40960, 'FlashpixVersion'),
                (40961, 'ColorSpace'),
                (40962, 'PixelXDimension'),
                (40963, 'PixelYDimension'),
                (40964, 'RelatedSoundFile'),
                (40976, 'SamsungRawPointersOffset'),
                (40977, 'SamsungRawPointersLength'),
                (41217, 'SamsungRawByteOrder'),
                (41218, 'SamsungRawUnknown'),
                (41483, 'FlashEnergy'),
                (41484, 'SpatialFrequencyResponse'),
                (41485, 'Noise'),  # 37389
                (41486, 'FocalPlaneXResolution'),  # 37390
                (41487, 'FocalPlaneYResolution'),  # 37391
                (41488, 'FocalPlaneResolutionUnit'),  # 37392
                (41489, 'ImageNumber'),  # 37393
                (41490, 'SecurityClassification'),  # 37394
                (41491, 'ImageHistory'),  # 37395
                (41492, 'SubjectLocation'),  # 37395
                (41493, 'ExposureIndex '),  # 37397
                (41494, 'TIFF-EPStandardID'),
                (41495, 'SensingMethod'),  # 37399
                (41728, 'FileSource'),
                (41729, 'SceneType'),
                (41730, 'CFAPattern'),  # 33422
                (41985, 'CustomRendered'),
                (41986, 'ExposureMode'),
                (41987, 'WhiteBalance'),
                (41988, 'DigitalZoomRatio'),
                (41989, 'FocalLengthIn35mmFilm'),
                (41990, 'SceneCaptureType'),
                (41991, 'GainControl'),
                (41992, 'Contrast'),
                (41993, 'Saturation'),
                (41994, 'Sharpness'),
                (41995, 'DeviceSettingDescription'),
                (41996, 'SubjectDistanceRange'),
                (42016, 'ImageUniqueID'),
                (42032, 'CameraOwnerName'),
                (42033, 'BodySerialNumber'),
                (42034, 'LensSpecification'),
                (42035, 'LensMake'),
                (42036, 'LensModel'),
                (42037, 'LensSerialNumber'),
                (42080, 'CompositeImage'),
                (42081, 'SourceImageNumberCompositeImage'),
                (42082, 'SourceExposureTimesCompositeImage'),
                (42112, 'GDAL_METADATA'),
                (42113, 'GDAL_NODATA'),
                (42240, 'Gamma'),
                (43314, 'NIHImageHeader'),
                (44992, 'ExpandSoftware'),
                (44993, 'ExpandLens'),
                (44994, 'ExpandFilm'),
                (44995, 'ExpandFilterLens'),
                (44996, 'ExpandScanner'),
                (44997, 'ExpandFlashLamp'),
                (48129, 'PixelFormat'),  # HDP and WDP
                (48130, 'Transformation'),
                (48131, 'Uncompressed'),
                (48132, 'ImageType'),
                (48256, 'ImageWidth'),  # 256
                (48257, 'ImageHeight'),
                (48258, 'WidthResolution'),
                (48259, 'HeightResolution'),
                (48320, 'ImageOffset'),
                (48321, 'ImageByteCount'),
                (48322, 'AlphaOffset'),
                (48323, 'AlphaByteCount'),
                (48324, 'ImageDataDiscard'),
                (48325, 'AlphaDataDiscard'),
                (50003, 'KodakAPP3'),
                (50215, 'OceScanjobDescription'),
                (50216, 'OceApplicationSelector'),
                (50217, 'OceIdentificationNumber'),
                (50218, 'OceImageLogicCharacteristics'),
                (50255, 'Annotations'),
                (50288, 'MC_Id'),  # Media Cybernetics
                (50289, 'MC_XYPosition'),
                (50290, 'MC_ZPosition'),
                (50291, 'MC_XYCalibration'),
                (50292, 'MC_LensCharacteristics'),
                (50293, 'MC_ChannelName'),
                (50294, 'MC_ExcitationWavelength'),
                (50295, 'MC_TimeStamp'),
                (50296, 'MC_FrameProperties'),
                (50341, 'PrintImageMatching'),
                (50495, 'PCO_RAW'),  # TODO, PCO CamWare
                (50547, 'OriginalFileName'),
                (50560, 'USPTO_OriginalContentType'),  # US Patent Office
                (50561, 'USPTO_RotationCode'),
                (50648, 'CR2Unknown1'),
                (50649, 'CR2Unknown2'),
                (50656, 'CR2CFAPattern'),
                (50674, 'LercParameters'),  # ESGI 50674 .. 50677
                (50706, 'DNGVersion'),  # DNG 50706 .. 51114
                (50707, 'DNGBackwardVersion'),
                (50708, 'UniqueCameraModel'),
                (50709, 'LocalizedCameraModel'),
                (50710, 'CFAPlaneColor'),
                (50711, 'CFALayout'),
                (50712, 'LinearizationTable'),
                (50713, 'BlackLevelRepeatDim'),
                (50714, 'BlackLevel'),
                (50715, 'BlackLevelDeltaH'),
                (50716, 'BlackLevelDeltaV'),
                (50717, 'WhiteLevel'),
                (50718, 'DefaultScale'),
                (50719, 'DefaultCropOrigin'),
                (50720, 'DefaultCropSize'),
                (50721, 'ColorMatrix1'),
                (50722, 'ColorMatrix2'),
                (50723, 'CameraCalibration1'),
                (50724, 'CameraCalibration2'),
                (50725, 'ReductionMatrix1'),
                (50726, 'ReductionMatrix2'),
                (50727, 'AnalogBalance'),
                (50728, 'AsShotNeutral'),
                (50729, 'AsShotWhiteXY'),
                (50730, 'BaselineExposure'),
                (50731, 'BaselineNoise'),
                (50732, 'BaselineSharpness'),
                (50733, 'BayerGreenSplit'),
                (50734, 'LinearResponseLimit'),
                (50735, 'CameraSerialNumber'),
                (50736, 'LensInfo'),
                (50737, 'ChromaBlurRadius'),
                (50738, 'AntiAliasStrength'),
                (50739, 'ShadowScale'),
                (50740, 'DNGPrivateData'),
                (50741, 'MakerNoteSafety'),
                (50752, 'RawImageSegmentation'),
                (50778, 'CalibrationIlluminant1'),
                (50779, 'CalibrationIlluminant2'),
                (50780, 'BestQualityScale'),
                (50781, 'RawDataUniqueID'),
                (50784, 'AliasLayerMetadata'),
                (50827, 'OriginalRawFileName'),
                (50828, 'OriginalRawFileData'),
                (50829, 'ActiveArea'),
                (50830, 'MaskedAreas'),
                (50831, 'AsShotICCProfile'),
                (50832, 'AsShotPreProfileMatrix'),
                (50833, 'CurrentICCProfile'),
                (50834, 'CurrentPreProfileMatrix'),
                (50838, 'IJMetadataByteCounts'),
                (50839, 'IJMetadata'),
                (50844, 'RPCCoefficientTag'),
                (50879, 'ColorimetricReference'),
                (50885, 'SRawType'),
                (50898, 'PanasonicTitle'),
                (50899, 'PanasonicTitle2'),
                (50908, 'RSID'),  # DGIWG
                (50909, 'GEO_METADATA'),  # DGIWG XML
                (50931, 'CameraCalibrationSignature'),
                (50932, 'ProfileCalibrationSignature'),
                (50933, 'ProfileIFD'),  # EXTRACAMERAPROFILES
                (50934, 'AsShotProfileName'),
                (50935, 'NoiseReductionApplied'),
                (50936, 'ProfileName'),
                (50937, 'ProfileHueSatMapDims'),
                (50938, 'ProfileHueSatMapData1'),
                (50939, 'ProfileHueSatMapData2'),
                (50940, 'ProfileToneCurve'),
                (50941, 'ProfileEmbedPolicy'),
                (50942, 'ProfileCopyright'),
                (50964, 'ForwardMatrix1'),
                (50965, 'ForwardMatrix2'),
                (50966, 'PreviewApplicationName'),
                (50967, 'PreviewApplicationVersion'),
                (50968, 'PreviewSettingsName'),
                (50969, 'PreviewSettingsDigest'),
                (50970, 'PreviewColorSpace'),
                (50971, 'PreviewDateTime'),
                (50972, 'RawImageDigest'),
                (50973, 'OriginalRawFileDigest'),
                (50974, 'SubTileBlockSize'),
                (50975, 'RowInterleaveFactor'),
                (50981, 'ProfileLookTableDims'),
                (50982, 'ProfileLookTableData'),
                (51008, 'OpcodeList1'),
                (51009, 'OpcodeList2'),
                (51022, 'OpcodeList3'),
                (51023, 'FibicsXML'),  #
                (51041, 'NoiseProfile'),
                (51043, 'TimeCodes'),
                (51044, 'FrameRate'),
                (51058, 'TStop'),
                (51081, 'ReelName'),
                (51089, 'OriginalDefaultFinalSize'),
                (51090, 'OriginalBestQualitySize'),
                (51091, 'OriginalDefaultCropSize'),
                (51105, 'CameraLabel'),
                (51107, 'ProfileHueSatMapEncoding'),
                (51108, 'ProfileLookTableEncoding'),
                (51109, 'BaselineExposureOffset'),
                (51110, 'DefaultBlackRender'),
                (51111, 'NewRawImageDigest'),
                (51112, 'RawToPreviewGain'),
                (51113, 'CacheBlob'),
                (51114, 'CacheVersion'),
                (51123, 'MicroManagerMetadata'),
                (51125, 'DefaultUserCrop'),
                (51159, 'ZIFmetadata'),  # Objective Pathology Services
                (51160, 'ZIFannotations'),  # Objective Pathology Services
                (51177, 'DepthFormat'),
                (51178, 'DepthNear'),
                (51179, 'DepthFar'),
                (51180, 'DepthUnits'),
                (51181, 'DepthMeasureType'),
                (51182, 'EnhanceParams'),
                (52525, 'ProfileGainTableMap'),  # DNG 1.6
                (52526, 'SemanticName'),  # DNG 1.6
                (52528, 'SemanticInstanceID'),  # DNG 1.6
                (52536, 'MaskSubArea'),  # DNG 1.6
                (52543, 'RGBTables'),  # DNG 1.6
                (52529, 'CalibrationIlluminant3'),  # DNG 1.6
                (52531, 'ColorMatrix3'),  # DNG 1.6
                (52530, 'CameraCalibration3'),  # DNG 1.6
                (52538, 'ReductionMatrix3'),  # DNG 1.6
                (52537, 'ProfileHueSatMapData3'),  # DNG 1.6
                (52532, 'ForwardMatrix3'),  # DNG 1.6
                (52533, 'IlluminantData1'),  # DNG 1.6
                (52534, 'IlluminantData2'),  # DNG 1.6
                (53535, 'IlluminantData3'),  # DNG 1.6
                (52544, 'ProfileGainTableMap2'),  # DNG 1.7
                (52547, 'ColumnInterleaveFactor'),  # DNG 1.7
                (52548, 'ImageSequenceInfo'),  # DNG 1.7
                (52550, 'ImageStats'),  # DNG 1.7
                (52551, 'ProfileDynamicRange'),  # DNG 1.7
                (52552, 'ProfileGroupName'),  # DNG 1.7
                (52553, 'JXLDistance'),  # DNG 1.7
                (52554, 'JXLEffort'),  # DNG 1.7
                (52555, 'JXLDecodeSpeed'),  # DNG 1.7
                (55000, 'AperioUnknown55000'),
                (55001, 'AperioMagnification'),
                (55002, 'AperioMPP'),
                (55003, 'AperioScanScopeID'),
                (55004, 'AperioDate'),
                (59932, 'Padding'),
                (59933, 'OffsetSchema'),
                # Reusable Tags 65000-65535
                # (65000, 'DimapDocumentXML'),
                # EER metadata:
                # (65001, 'AcquisitionMetadata'),
                # (65002, 'FrameMetadata'),
                # (65006, 'ImageMetadata'),
                # (65007, 'PosSkipBits'),
                # (65008, 'HorzSubBits'),
                # (65009, 'VertSubBits'),
                # Photoshop Camera RAW EXIF tags:
                # (65000, 'OwnerName'),
                # (65001, 'SerialNumber'),
                # (65002, 'Lens'),
                # (65024, 'KodakKDCPrivateIFD'),
                # (65100, 'RawFile'),
                # (65101, 'Converter'),
                # (65102, 'WhiteBalance'),
                # (65105, 'Exposure'),
                # (65106, 'Shadows'),
                # (65107, 'Brightness'),
                # (65108, 'Contrast'),
                # (65109, 'Saturation'),
                # (65110, 'Sharpness'),
                # (65111, 'Smoothness'),
                # (65112, 'MoireFilter'),
                (65200, 'FlexXML'),
            )
        )

    @cached_property
    def TAG_READERS(
        self,
    ) -> dict[int, Callable[[FileHandle, ByteOrder, int, int, int], Any]]:
        # map tag codes to import functions
        return {
            301: read_colormap,
            320: read_colormap,
            # 700: read_bytes,  # read_utf8,
            # 34377: read_bytes,
            33723: read_bytes,
            # 34675: read_bytes,
            33628: read_uic1tag,  # Universal Imaging Corp STK
            33629: read_uic2tag,
            33630: read_uic3tag,
            33631: read_uic4tag,
            34118: read_cz_sem,  # Carl Zeiss SEM
            34361: read_mm_header,  # Olympus FluoView
            34362: read_mm_stamp,
            34363: read_numpy,  # MM_Unknown
            34386: read_numpy,  # MM_UserBlock
            34412: read_cz_lsminfo,  # Carl Zeiss LSM
            34680: read_fei_metadata,  # S-FEG
            34682: read_fei_metadata,  # Helios NanoLab
            37706: read_tvips_header,  # TVIPS EMMENU
            37724: read_bytes,  # ImageSourceData
            33923: read_bytes,  # read_leica_magic
            43314: read_nih_image_header,
            # 40001: read_bytes,
            40100: read_bytes,
            50288: read_bytes,
            50296: read_bytes,
            50839: read_bytes,
            51123: read_json,
            33471: read_sis_ini,
            33560: read_sis,
            34665: read_exif_ifd,
            34853: read_gps_ifd,  # conflicts with OlympusSIS
            40965: read_interoperability_ifd,
            65426: read_numpy,  # NDPI McuStarts
            65432: read_numpy,  # NDPI McuStartsHighBytes
            65439: read_numpy,  # NDPI unknown
            65459: read_bytes,  # NDPI bytes, not string
        }

    @cached_property
    def TAG_LOAD(self) -> frozenset[int]:
        # tags whose values are not delay loaded
        return frozenset(
            (
                258,  # BitsPerSample
                270,  # ImageDescription
                273,  # StripOffsets
                277,  # SamplesPerPixel
                279,  # StripByteCounts
                282,  # XResolution
                283,  # YResolution
                # 301,  # TransferFunction
                305,  # Software
                # 306,  # DateTime
                # 320,  # ColorMap
                324,  # TileOffsets
                325,  # TileByteCounts
                330,  # SubIFDs
                338,  # ExtraSamples
                339,  # SampleFormat
                347,  # JPEGTables
                513,  # JPEGInterchangeFormat
                514,  # JPEGInterchangeFormatLength
                530,  # YCbCrSubSampling
                33628,  # UIC1tag
                42113,  # GDAL_NODATA
                50838,  # IJMetadataByteCounts
                50839,  # IJMetadata
            )
        )

    @cached_property
    def TAG_FILTERED(self) -> frozenset[int]:
        # tags filtered from extratags in :py:meth:`TiffWriter.write`
        return frozenset(
            (
                256,  # ImageWidth
                257,  # ImageLength
                258,  # BitsPerSample
                259,  # Compression
                262,  # PhotometricInterpretation
                266,  # FillOrder
                273,  # StripOffsets
                277,  # SamplesPerPixel
                278,  # RowsPerStrip
                279,  # StripByteCounts
                284,  # PlanarConfiguration
                317,  # Predictor
                322,  # TileWidth
                323,  # TileLength
                324,  # TileOffsets
                325,  # TileByteCounts
                330,  # SubIFDs,
                338,  # ExtraSamples
                339,  # SampleFormat
                400,  # GlobalParametersIFD
                32997,  # ImageDepth
                32998,  # TileDepth
                34665,  # ExifTag
                34853,  # GPSTag
                40965,  # InteroperabilityTag
            )
        )

    @cached_property
    def TAG_TUPLE(self) -> frozenset[int]:
        # tags whose values must be stored as tuples
        return frozenset(
            (
                273,
                279,
                282,
                283,
                324,
                325,
                330,
                338,
                513,
                514,
                530,
                531,
                34736,
                50838,
            )
        )

    @cached_property
    def TAG_ATTRIBUTES(self) -> dict[int, str]:
        # map tag codes to TiffPage attribute names
        return {
            254: 'subfiletype',
            256: 'imagewidth',
            257: 'imagelength',
            # 258: 'bitspersample',  # set manually
            259: 'compression',
            262: 'photometric',
            266: 'fillorder',
            270: 'description',
            277: 'samplesperpixel',
            278: 'rowsperstrip',
            284: 'planarconfig',
            # 301: 'transferfunction',  # delay load
            305: 'software',
            # 320: 'colormap',  # delay load
            317: 'predictor',
            322: 'tilewidth',
            323: 'tilelength',
            330: 'subifds',
            338: 'extrasamples',
            # 339: 'sampleformat',  # set manually
            347: 'jpegtables',
            530: 'subsampling',
            32997: 'imagedepth',
            32998: 'tiledepth',
        }

    @cached_property
    def TAG_ENUM(self) -> dict[int, type[enum.Enum]]:
        # map tag codes to Enums
        return {
            254: FILETYPE,
            255: OFILETYPE,
            259: COMPRESSION,
            262: PHOTOMETRIC,
            # 263: THRESHOLD,
            266: FILLORDER,
            274: ORIENTATION,
            284: PLANARCONFIG,
            # 290: GRAYRESPONSEUNIT,
            # 292: GROUP3OPT
            # 293: GROUP4OPT
            296: RESUNIT,
            # 300: COLORRESPONSEUNIT,
            317: PREDICTOR,
            338: EXTRASAMPLE,
            339: SAMPLEFORMAT,
            # 512: JPEGPROC
            # 531: YCBCRPOSITION
        }

    @cached_property
    def EXIF_TAGS(self) -> TiffTagRegistry:
        """Registry of EXIF tags, including private Photoshop Camera RAW."""
        # 65000 - 65112  Photoshop Camera RAW EXIF tags
        tags = TiffTagRegistry(
            (
                (65000, 'OwnerName'),
                (65001, 'SerialNumber'),
                (65002, 'Lens'),
                (65100, 'RawFile'),
                (65101, 'Converter'),
                (65102, 'WhiteBalance'),
                (65105, 'Exposure'),
                (65106, 'Shadows'),
                (65107, 'Brightness'),
                (65108, 'Contrast'),
                (65109, 'Saturation'),
                (65110, 'Sharpness'),
                (65111, 'Smoothness'),
                (65112, 'MoireFilter'),
            )
        )
        tags.update(TIFF.TAGS)
        return tags

    @cached_property
    def NDPI_TAGS(self) -> TiffTagRegistry:
        """Registry of private TIFF tags for Hamamatsu NDPI (65420-65458)."""
        # TODO: obtain specification
        return TiffTagRegistry(
            (
                (65324, 'OffsetHighBytes'),
                (65325, 'ByteCountHighBytes'),
                (65420, 'FileFormat'),
                (65421, 'Magnification'),  # SourceLens
                (65422, 'XOffsetFromSlideCenter'),
                (65423, 'YOffsetFromSlideCenter'),
                (65424, 'ZOffsetFromSlideCenter'),  # FocalPlane
                (65425, 'TissueIndex'),
                (65426, 'McuStarts'),
                (65427, 'SlideLabel'),
                (65428, 'AuthCode'),  # ?
                (65429, '65429'),
                (65430, '65430'),
                (65431, '65431'),
                (65432, 'McuStartsHighBytes'),
                (65433, '65433'),
                (65434, 'Fluorescence'),  # FilterSetName, Channel
                (65435, 'ExposureRatio'),
                (65436, 'RedMultiplier'),
                (65437, 'GreenMultiplier'),
                (65438, 'BlueMultiplier'),
                (65439, 'FocusPoints'),
                (65440, 'FocusPointRegions'),
                (65441, 'CaptureMode'),
                (65442, 'ScannerSerialNumber'),
                (65443, '65443'),
                (65444, 'JpegQuality'),
                (65445, 'RefocusInterval'),
                (65446, 'FocusOffset'),
                (65447, 'BlankLines'),
                (65448, 'FirmwareVersion'),
                (65449, 'Comments'),  # PropertyMap, CalibrationInfo
                (65450, 'LabelObscured'),
                (65451, 'Wavelength'),
                (65452, '65452'),
                (65453, 'LampAge'),
                (65454, 'ExposureTime'),
                (65455, 'FocusTime'),
                (65456, 'ScanTime'),
                (65457, 'WriteTime'),
                (65458, 'FullyAutoFocus'),
                (65500, 'DefaultGamma'),
            )
        )

    @cached_property
    def GPS_TAGS(self) -> TiffTagRegistry:
        """Registry of GPS IFD tags."""
        return TiffTagRegistry(
            (
                (0, 'GPSVersionID'),
                (1, 'GPSLatitudeRef'),
                (2, 'GPSLatitude'),
                (3, 'GPSLongitudeRef'),
                (4, 'GPSLongitude'),
                (5, 'GPSAltitudeRef'),
                (6, 'GPSAltitude'),
                (7, 'GPSTimeStamp'),
                (8, 'GPSSatellites'),
                (9, 'GPSStatus'),
                (10, 'GPSMeasureMode'),
                (11, 'GPSDOP'),
                (12, 'GPSSpeedRef'),
                (13, 'GPSSpeed'),
                (14, 'GPSTrackRef'),
                (15, 'GPSTrack'),
                (16, 'GPSImgDirectionRef'),
                (17, 'GPSImgDirection'),
                (18, 'GPSMapDatum'),
                (19, 'GPSDestLatitudeRef'),
                (20, 'GPSDestLatitude'),
                (21, 'GPSDestLongitudeRef'),
                (22, 'GPSDestLongitude'),
                (23, 'GPSDestBearingRef'),
                (24, 'GPSDestBearing'),
                (25, 'GPSDestDistanceRef'),
                (26, 'GPSDestDistance'),
                (27, 'GPSProcessingMethod'),
                (28, 'GPSAreaInformation'),
                (29, 'GPSDateStamp'),
                (30, 'GPSDifferential'),
                (31, 'GPSHPositioningError'),
            )
        )

    @cached_property
    def IOP_TAGS(self) -> TiffTagRegistry:
        """Registry of Interoperability IFD tags."""
        return TiffTagRegistry(
            (
                (1, 'InteroperabilityIndex'),
                (2, 'InteroperabilityVersion'),
                (4096, 'RelatedImageFileFormat'),
                (4097, 'RelatedImageWidth'),
                (4098, 'RelatedImageLength'),
            )
        )

    @cached_property
    def PHOTOMETRIC_SAMPLES(self) -> dict[int, int]:
        """Map :py:class:`PHOTOMETRIC` to number of photometric samples."""
        return {
            0: 1,  # MINISWHITE
            1: 1,  # MINISBLACK
            2: 3,  # RGB
            3: 1,  # PALETTE
            4: 1,  # MASK
            5: 4,  # SEPARATED
            6: 3,  # YCBCR
            8: 3,  # CIELAB
            9: 3,  # ICCLAB
            10: 3,  # ITULAB
            32803: 1,  # CFA
            32844: 1,  # LOGL ?
            32845: 3,  # LOGLUV
            34892: 3,  # LINEAR_RAW ?
            51177: 1,  # DEPTH_MAP ?
            52527: 1,  # SEMANTIC_MASK ?
        }

    @cached_property
    def DATA_FORMATS(self) -> dict[int, str]:
        """Map :py:class:`DATATYPE` to Python struct formats."""
        return {
            1: '1B',
            2: '1s',
            3: '1H',
            4: '1I',
            5: '2I',
            6: '1b',
            7: '1B',
            8: '1h',
            9: '1i',
            10: '2i',
            11: '1f',
            12: '1d',
            13: '1I',
            # 14: '',
            # 15: '',
            16: '1Q',
            17: '1q',
            18: '1Q',
        }

    @cached_property
    def DATA_DTYPES(self) -> dict[str, int]:
        """Map NumPy dtype to :py:class:`DATATYPE`."""
        return {
            'B': 1,
            's': 2,
            'H': 3,
            'I': 4,
            '2I': 5,
            'b': 6,
            'h': 8,
            'i': 9,
            '2i': 10,
            'f': 11,
            'd': 12,
            'Q': 16,
            'q': 17,
        }

    @cached_property
    def SAMPLE_DTYPES(self) -> dict[tuple[int, int | tuple[int, ...]], str]:
        """Map :py:class:`SAMPLEFORMAT` and BitsPerSample to NumPy dtype."""
        return {
            # UINT
            (1, 1): '?',  # bitmap
            (1, 2): 'B',
            (1, 3): 'B',
            (1, 4): 'B',
            (1, 5): 'B',
            (1, 6): 'B',
            (1, 7): 'B',
            (1, 8): 'B',
            (1, 9): 'H',
            (1, 10): 'H',
            (1, 11): 'H',
            (1, 12): 'H',
            (1, 13): 'H',
            (1, 14): 'H',
            (1, 15): 'H',
            (1, 16): 'H',
            (1, 17): 'I',
            (1, 18): 'I',
            (1, 19): 'I',
            (1, 20): 'I',
            (1, 21): 'I',
            (1, 22): 'I',
            (1, 23): 'I',
            (1, 24): 'I',
            (1, 25): 'I',
            (1, 26): 'I',
            (1, 27): 'I',
            (1, 28): 'I',
            (1, 29): 'I',
            (1, 30): 'I',
            (1, 31): 'I',
            (1, 32): 'I',
            (1, 64): 'Q',
            # VOID : treat as UINT
            (4, 1): '?',  # bitmap
            (4, 2): 'B',
            (4, 3): 'B',
            (4, 4): 'B',
            (4, 5): 'B',
            (4, 6): 'B',
            (4, 7): 'B',
            (4, 8): 'B',
            (4, 9): 'H',
            (4, 10): 'H',
            (4, 11): 'H',
            (4, 12): 'H',
            (4, 13): 'H',
            (4, 14): 'H',
            (4, 15): 'H',
            (4, 16): 'H',
            (4, 17): 'I',
            (4, 18): 'I',
            (4, 19): 'I',
            (4, 20): 'I',
            (4, 21): 'I',
            (4, 22): 'I',
            (4, 23): 'I',
            (4, 24): 'I',
            (4, 25): 'I',
            (4, 26): 'I',
            (4, 27): 'I',
            (4, 28): 'I',
            (4, 29): 'I',
            (4, 30): 'I',
            (4, 31): 'I',
            (4, 32): 'I',
            (4, 64): 'Q',
            # INT
            (2, 8): 'b',
            (2, 16): 'h',
            (2, 32): 'i',
            (2, 64): 'q',
            # IEEEFP
            (3, 16): 'e',
            (3, 24): 'f',  # float24 bit not supported by numpy
            (3, 32): 'f',
            (3, 64): 'd',
            # COMPLEXIEEEFP
            (6, 64): 'F',
            (6, 128): 'D',
            # RGB565
            (1, (5, 6, 5)): 'B',
            # COMPLEXINT : not supported by numpy
            (5, 16): 'E',
            (5, 32): 'F',
            (5, 64): 'D',
        }

    @cached_property
    def PREDICTORS(self) -> Mapping[int, Callable[..., Any]]:
        """Map :py:class:`PREDICTOR` value to encode function."""
        return PredictorCodec(True)

    @cached_property
    def UNPREDICTORS(self) -> Mapping[int, Callable[..., Any]]:
        """Map :py:class:`PREDICTOR` value to decode function."""
        return PredictorCodec(False)

    @cached_property
    def COMPRESSORS(self) -> Mapping[int, Callable[..., Any]]:
        """Map :py:class:`COMPRESSION` value to compress function."""
        return CompressionCodec(True)

    @cached_property
    def DECOMPRESSORS(self) -> Mapping[int, Callable[..., Any]]:
        """Map :py:class:`COMPRESSION` value to decompress function."""
        return CompressionCodec(False)

    @cached_property
    def IMAGE_COMPRESSIONS(self) -> set[int]:
        # set of compression to encode/decode images
        # encode/decode preserves shape and dtype
        # cannot be used with predictors or fillorder
        return {
            6,  # jpeg
            7,  # jpeg
            22610,  # jpegxr
            33003,  # jpeg2k
            33004,  # jpeg2k
            33005,  # jpeg2k
            33007,  # alt_jpeg
            34712,  # jpeg2k
            34892,  # jpeg
            34933,  # png
            34934,  # jpegxr ZIF
            48124,  # jetraw
            50001,  # webp
            50002,  # jpegxl
            52546,  # jpegxl DNG
        }

    @cached_property
    def AXES_NAMES(self) -> dict[str, str]:
        """Map axes character codes to dimension names.

        - **X : width** (image width)
        - **Y : height** (image length)
        - **Z : depth** (image depth)
        - **S : sample** (color space and extra samples)
        - **I : sequence** (generic sequence of images, frames, planes, pages)
        - **T : time** (time series)
        - **C : channel** (acquisition path or emission wavelength)
        - **A : angle** (OME)
        - **P : phase** (OME. In LSM, **P** maps to **position**)
        - **R : tile** (OME. Region, position, or mosaic)
        - **H : lifetime** (OME. Histogram)
        - **E : lambda** (OME. Excitation wavelength)
        - **Q : other** (OME)
        - **L : exposure** (FluoView)
        - **V : event** (FluoView)
        - **M : mosaic** (LSM 6)
        - **J : column** (NDTiff)
        - **K : row** (NDTiff)

        There is no universal standard for dimension codes or names.
        This mapping mainly follows TIFF, OME-TIFF, ImageJ, LSM, and FluoView
        conventions.

        """
        return {
            'X': 'width',
            'Y': 'height',
            'Z': 'depth',
            'S': 'sample',
            'I': 'sequence',
            # 'F': 'file',
            'T': 'time',
            'C': 'channel',
            'A': 'angle',
            'P': 'phase',
            'R': 'tile',
            'H': 'lifetime',
            'E': 'lambda',
            'L': 'exposure',
            'V': 'event',
            'M': 'mosaic',
            'Q': 'other',
            'J': 'column',
            'K': 'row',
        }

    @cached_property
    def AXES_CODES(self) -> dict[str, str]:
        """Map dimension names to axes character codes.

        Reverse mapping of :py:attr:`AXES_NAMES`.

        """
        codes = {name: code for code, name in TIFF.AXES_NAMES.items()}
        codes['z'] = 'Z'  # NDTiff
        codes['position'] = 'R'  # NDTiff
        return codes

    @cached_property
    def GEO_KEYS(self) -> type[enum.IntEnum]:
        """:py:class:`geodb.GeoKeys`."""
        try:
            from .geodb import GeoKeys
        except ImportError:

            class GeoKeys(enum.IntEnum):  # type: ignore[no-redef]
                pass

        return GeoKeys

    @cached_property
    def GEO_CODES(self) -> dict[int, type[enum.IntEnum]]:
        """Map :py:class:`geodb.GeoKeys` to GeoTIFF codes."""
        try:
            from .geodb import GEO_CODES
        except ImportError:
            GEO_CODES = {}
        return GEO_CODES

    @cached_property
    def PAGE_FLAGS(self) -> set[str]:
        # TiffFile and TiffPage 'is_\*' attributes
        exclude = {
            'reduced',
            'mask',
            'final',
            'memmappable',
            'contiguous',
            'tiled',
            'subsampled',
            'jfif',
        }
        return {
            a[3:]
            for a in dir(TiffPage)
            if a[:3] == 'is_' and a[3:] not in exclude
        }

    @cached_property
    def FILE_FLAGS(self) -> set[str]:
        # TiffFile 'is_\*' attributes
        exclude = {'bigtiff', 'appendable'}
        return {
            a[3:]
            for a in dir(TiffFile)
            if a[:3] == 'is_' and a[3:] not in exclude
        }.union(TIFF.PAGE_FLAGS)

    @property
    def FILE_PATTERNS(self) -> dict[str, str]:
        # predefined FileSequence patterns
        return {
            'axes': r"""(?ix)
                # matches Olympus OIF and Leica TIFF series
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                _?(?:(q|l|p|a|c|t|x|y|z|ch|tp)(\d{1,4}))?
                """
        }

    @property
    def FILE_EXTENSIONS(self) -> tuple[str, ...]:
        """Known TIFF file extensions."""
        return (
            'tif',
            'tiff',
            'ome.tif',
            'lsm',
            'stk',
            'qpi',
            'pcoraw',
            'qptiff',
            'ptiff',
            'ptif',
            'gel',
            'seq',
            'svs',
            'avs',
            'scn',
            'zif',
            'ndpi',
            'bif',
            'tf8',
            'tf2',
            'btf',
            'eer',
        )

    @property
    def FILEOPEN_FILTER(self) -> list[tuple[str, str]]:
        # string for use in Windows File Open box
        return [
            (f'{ext.upper()} files', f'*.{ext}')
            for ext in TIFF.FILE_EXTENSIONS
        ] + [('All files', '*')]

    @property
    def CZ_LSMINFO(self) -> list[tuple[str, str]]:
        # numpy data type of LSMINFO structure
        return [
            ('MagicNumber', 'u4'),
            ('StructureSize', 'i4'),
            ('DimensionX', 'i4'),
            ('DimensionY', 'i4'),
            ('DimensionZ', 'i4'),
            ('DimensionChannels', 'i4'),
            ('DimensionTime', 'i4'),
            ('DataType', 'i4'),  # DATATYPES
            ('ThumbnailX', 'i4'),
            ('ThumbnailY', 'i4'),
            ('VoxelSizeX', 'f8'),
            ('VoxelSizeY', 'f8'),
            ('VoxelSizeZ', 'f8'),
            ('OriginX', 'f8'),
            ('OriginY', 'f8'),
            ('OriginZ', 'f8'),
            ('ScanType', 'u2'),
            ('SpectralScan', 'u2'),
            ('TypeOfData', 'u4'),  # TYPEOFDATA
            ('OffsetVectorOverlay', 'u4'),
            ('OffsetInputLut', 'u4'),
            ('OffsetOutputLut', 'u4'),
            ('OffsetChannelColors', 'u4'),
            ('TimeIntervall', 'f8'),
            ('OffsetChannelDataTypes', 'u4'),
            ('OffsetScanInformation', 'u4'),  # SCANINFO
            ('OffsetKsData', 'u4'),
            ('OffsetTimeStamps', 'u4'),
            ('OffsetEventList', 'u4'),
            ('OffsetRoi', 'u4'),
            ('OffsetBleachRoi', 'u4'),
            ('OffsetNextRecording', 'u4'),
            # LSM 2.0 ends here
            ('DisplayAspectX', 'f8'),
            ('DisplayAspectY', 'f8'),
            ('DisplayAspectZ', 'f8'),
            ('DisplayAspectTime', 'f8'),
            ('OffsetMeanOfRoisOverlay', 'u4'),
            ('OffsetTopoIsolineOverlay', 'u4'),
            ('OffsetTopoProfileOverlay', 'u4'),
            ('OffsetLinescanOverlay', 'u4'),
            ('ToolbarFlags', 'u4'),
            ('OffsetChannelWavelength', 'u4'),
            ('OffsetChannelFactors', 'u4'),
            ('ObjectiveSphereCorrection', 'f8'),
            ('OffsetUnmixParameters', 'u4'),
            # LSM 3.2, 4.0 end here
            ('OffsetAcquisitionParameters', 'u4'),
            ('OffsetCharacteristics', 'u4'),
            ('OffsetPalette', 'u4'),
            ('TimeDifferenceX', 'f8'),
            ('TimeDifferenceY', 'f8'),
            ('TimeDifferenceZ', 'f8'),
            ('InternalUse1', 'u4'),
            ('DimensionP', 'i4'),
            ('DimensionM', 'i4'),
            ('DimensionsReserved', '16i4'),
            ('OffsetTilePositions', 'u4'),
            ('', '9u4'),  # Reserved
            ('OffsetPositions', 'u4'),
            # ('', '21u4'),  # must be 0
        ]

    @property
    def CZ_LSMINFO_READERS(
        self,
    ) -> dict[str, Callable[[FileHandle], Any] | None]:
        # import functions for CZ_LSMINFO sub-records
        # TODO: read more CZ_LSMINFO sub-records
        return {
            'ScanInformation': read_lsm_scaninfo,
            'TimeStamps': read_lsm_timestamps,
            'EventList': read_lsm_eventlist,
            'ChannelColors': read_lsm_channelcolors,
            'Positions': read_lsm_positions,
            'TilePositions': read_lsm_positions,
            'VectorOverlay': None,
            'InputLut': read_lsm_lookuptable,
            'OutputLut': read_lsm_lookuptable,
            'TimeIntervall': None,
            'ChannelDataTypes': read_lsm_channeldatatypes,
            'KsData': None,
            'Roi': None,
            'BleachRoi': None,
            'NextRecording': None,  # read with TiffFile(fh, offset=)
            'MeanOfRoisOverlay': None,
            'TopoIsolineOverlay': None,
            'TopoProfileOverlay': None,
            'ChannelWavelength': read_lsm_channelwavelength,
            'SphereCorrection': None,
            'ChannelFactors': None,
            'UnmixParameters': None,
            'AcquisitionParameters': None,
            'Characteristics': None,
        }

    @property
    def CZ_LSMINFO_SCANTYPE(self) -> dict[int, str]:
        # map CZ_LSMINFO.ScanType to dimension order
        return {
            0: 'ZCYX',  # Stack, normal x-y-z-scan
            1: 'CZX',  # Z-Scan, x-z-plane
            2: 'CTX',  # Line or Time Series Line
            3: 'TCYX',  # Time Series Plane, x-y
            4: 'TCZX',  # Time Series z-Scan, x-z
            5: 'CTX',  # Time Series Mean-of-ROIs
            6: 'TZCYX',  # Time Series Stack, x-y-z
            7: 'TZCYX',  # TODO: Spline Scan
            8: 'CZX',  # Spline Plane, x-z
            9: 'TCZX',  # Time Series Spline Plane, x-z
            10: 'CTX',  # Point or Time Series Point
        }

    @property
    def CZ_LSMINFO_DIMENSIONS(self) -> dict[str, str]:
        # map dimension codes to CZ_LSMINFO attribute
        return {
            'X': 'DimensionX',
            'Y': 'DimensionY',
            'Z': 'DimensionZ',
            'C': 'DimensionChannels',
            'T': 'DimensionTime',
            'P': 'DimensionP',
            'M': 'DimensionM',
        }

    @property
    def CZ_LSMINFO_DATATYPES(self) -> dict[int, str]:
        # description of CZ_LSMINFO.DataType
        return {
            0: 'varying data types',
            1: '8 bit unsigned integer',
            2: '12 bit unsigned integer',
            5: '32 bit float',
        }

    @property
    def CZ_LSMINFO_TYPEOFDATA(self) -> dict[int, str]:
        # description of CZ_LSMINFO.TypeOfData
        return {
            0: 'Original scan data',
            1: 'Calculated data',
            2: '3D reconstruction',
            3: 'Topography height map',
        }

    @property
    def CZ_LSMINFO_SCANINFO_ARRAYS(self) -> dict[int, str]:
        return {
            0x20000000: 'Tracks',
            0x30000000: 'Lasers',
            0x60000000: 'DetectionChannels',
            0x80000000: 'IlluminationChannels',
            0xA0000000: 'BeamSplitters',
            0xC0000000: 'DataChannels',
            0x11000000: 'Timers',
            0x13000000: 'Markers',
        }

    @property
    def CZ_LSMINFO_SCANINFO_STRUCTS(self) -> dict[int, str]:
        return {
            # 0x10000000: 'Recording',
            0x40000000: 'Track',
            0x50000000: 'Laser',
            0x70000000: 'DetectionChannel',
            0x90000000: 'IlluminationChannel',
            0xB0000000: 'BeamSplitter',
            0xD0000000: 'DataChannel',
            0x12000000: 'Timer',
            0x14000000: 'Marker',
        }

    @property
    def CZ_LSMINFO_SCANINFO_ATTRIBUTES(self) -> dict[int, str]:
        return {
            # Recording
            0x10000001: 'Name',
            0x10000002: 'Description',
            0x10000003: 'Notes',
            0x10000004: 'Objective',
            0x10000005: 'ProcessingSummary',
            0x10000006: 'SpecialScanMode',
            0x10000007: 'ScanType',
            0x10000008: 'ScanMode',
            0x10000009: 'NumberOfStacks',
            0x1000000A: 'LinesPerPlane',
            0x1000000B: 'SamplesPerLine',
            0x1000000C: 'PlanesPerVolume',
            0x1000000D: 'ImagesWidth',
            0x1000000E: 'ImagesHeight',
            0x1000000F: 'ImagesNumberPlanes',
            0x10000010: 'ImagesNumberStacks',
            0x10000011: 'ImagesNumberChannels',
            0x10000012: 'LinscanXySize',
            0x10000013: 'ScanDirection',
            0x10000014: 'TimeSeries',
            0x10000015: 'OriginalScanData',
            0x10000016: 'ZoomX',
            0x10000017: 'ZoomY',
            0x10000018: 'ZoomZ',
            0x10000019: 'Sample0X',
            0x1000001A: 'Sample0Y',
            0x1000001B: 'Sample0Z',
            0x1000001C: 'SampleSpacing',
            0x1000001D: 'LineSpacing',
            0x1000001E: 'PlaneSpacing',
            0x1000001F: 'PlaneWidth',
            0x10000020: 'PlaneHeight',
            0x10000021: 'VolumeDepth',
            0x10000023: 'Nutation',
            0x10000034: 'Rotation',
            0x10000035: 'Precession',
            0x10000036: 'Sample0time',
            0x10000037: 'StartScanTriggerIn',
            0x10000038: 'StartScanTriggerOut',
            0x10000039: 'StartScanEvent',
            0x10000040: 'StartScanTime',
            0x10000041: 'StopScanTriggerIn',
            0x10000042: 'StopScanTriggerOut',
            0x10000043: 'StopScanEvent',
            0x10000044: 'StopScanTime',
            0x10000045: 'UseRois',
            0x10000046: 'UseReducedMemoryRois',
            0x10000047: 'User',
            0x10000048: 'UseBcCorrection',
            0x10000049: 'PositionBcCorrection1',
            0x10000050: 'PositionBcCorrection2',
            0x10000051: 'InterpolationY',
            0x10000052: 'CameraBinning',
            0x10000053: 'CameraSupersampling',
            0x10000054: 'CameraFrameWidth',
            0x10000055: 'CameraFrameHeight',
            0x10000056: 'CameraOffsetX',
            0x10000057: 'CameraOffsetY',
            0x10000059: 'RtBinning',
            0x1000005A: 'RtFrameWidth',
            0x1000005B: 'RtFrameHeight',
            0x1000005C: 'RtRegionWidth',
            0x1000005D: 'RtRegionHeight',
            0x1000005E: 'RtOffsetX',
            0x1000005F: 'RtOffsetY',
            0x10000060: 'RtZoom',
            0x10000061: 'RtLinePeriod',
            0x10000062: 'Prescan',
            0x10000063: 'ScanDirectionZ',
            # Track
            0x40000001: 'MultiplexType',  # 0 After Line; 1 After Frame
            0x40000002: 'MultiplexOrder',
            0x40000003: 'SamplingMode',  # 0 Sample; 1 Line Avg; 2 Frame Avg
            0x40000004: 'SamplingMethod',  # 1 Mean; 2 Sum
            0x40000005: 'SamplingNumber',
            0x40000006: 'Acquire',
            0x40000007: 'SampleObservationTime',
            0x4000000B: 'TimeBetweenStacks',
            0x4000000C: 'Name',
            0x4000000D: 'Collimator1Name',
            0x4000000E: 'Collimator1Position',
            0x4000000F: 'Collimator2Name',
            0x40000010: 'Collimator2Position',
            0x40000011: 'IsBleachTrack',
            0x40000012: 'IsBleachAfterScanNumber',
            0x40000013: 'BleachScanNumber',
            0x40000014: 'TriggerIn',
            0x40000015: 'TriggerOut',
            0x40000016: 'IsRatioTrack',
            0x40000017: 'BleachCount',
            0x40000018: 'SpiCenterWavelength',
            0x40000019: 'PixelTime',
            0x40000021: 'CondensorFrontlens',
            0x40000023: 'FieldStopValue',
            0x40000024: 'IdCondensorAperture',
            0x40000025: 'CondensorAperture',
            0x40000026: 'IdCondensorRevolver',
            0x40000027: 'CondensorFilter',
            0x40000028: 'IdTransmissionFilter1',
            0x40000029: 'IdTransmission1',
            0x40000030: 'IdTransmissionFilter2',
            0x40000031: 'IdTransmission2',
            0x40000032: 'RepeatBleach',
            0x40000033: 'EnableSpotBleachPos',
            0x40000034: 'SpotBleachPosx',
            0x40000035: 'SpotBleachPosy',
            0x40000036: 'SpotBleachPosz',
            0x40000037: 'IdTubelens',
            0x40000038: 'IdTubelensPosition',
            0x40000039: 'TransmittedLight',
            0x4000003A: 'ReflectedLight',
            0x4000003B: 'SimultanGrabAndBleach',
            0x4000003C: 'BleachPixelTime',
            # Laser
            0x50000001: 'Name',
            0x50000002: 'Acquire',
            0x50000003: 'Power',
            # DetectionChannel
            0x70000001: 'IntegrationMode',
            0x70000002: 'SpecialMode',
            0x70000003: 'DetectorGainFirst',
            0x70000004: 'DetectorGainLast',
            0x70000005: 'AmplifierGainFirst',
            0x70000006: 'AmplifierGainLast',
            0x70000007: 'AmplifierOffsFirst',
            0x70000008: 'AmplifierOffsLast',
            0x70000009: 'PinholeDiameter',
            0x7000000A: 'CountingTrigger',
            0x7000000B: 'Acquire',
            0x7000000C: 'PointDetectorName',
            0x7000000D: 'AmplifierName',
            0x7000000E: 'PinholeName',
            0x7000000F: 'FilterSetName',
            0x70000010: 'FilterName',
            0x70000013: 'IntegratorName',
            0x70000014: 'ChannelName',
            0x70000015: 'DetectorGainBc1',
            0x70000016: 'DetectorGainBc2',
            0x70000017: 'AmplifierGainBc1',
            0x70000018: 'AmplifierGainBc2',
            0x70000019: 'AmplifierOffsetBc1',
            0x70000020: 'AmplifierOffsetBc2',
            0x70000021: 'SpectralScanChannels',
            0x70000022: 'SpiWavelengthStart',
            0x70000023: 'SpiWavelengthStop',
            0x70000026: 'DyeName',
            0x70000027: 'DyeFolder',
            # IlluminationChannel
            0x90000001: 'Name',
            0x90000002: 'Power',
            0x90000003: 'Wavelength',
            0x90000004: 'Aquire',
            0x90000005: 'DetchannelName',
            0x90000006: 'PowerBc1',
            0x90000007: 'PowerBc2',
            # BeamSplitter
            0xB0000001: 'FilterSet',
            0xB0000002: 'Filter',
            0xB0000003: 'Name',
            # DataChannel
            0xD0000001: 'Name',
            0xD0000003: 'Acquire',
            0xD0000004: 'Color',
            0xD0000005: 'SampleType',
            0xD0000006: 'BitsPerSample',
            0xD0000007: 'RatioType',
            0xD0000008: 'RatioTrack1',
            0xD0000009: 'RatioTrack2',
            0xD000000A: 'RatioChannel1',
            0xD000000B: 'RatioChannel2',
            0xD000000C: 'RatioConst1',
            0xD000000D: 'RatioConst2',
            0xD000000E: 'RatioConst3',
            0xD000000F: 'RatioConst4',
            0xD0000010: 'RatioConst5',
            0xD0000011: 'RatioConst6',
            0xD0000012: 'RatioFirstImages1',
            0xD0000013: 'RatioFirstImages2',
            0xD0000014: 'DyeName',
            0xD0000015: 'DyeFolder',
            0xD0000016: 'Spectrum',
            0xD0000017: 'Acquire',
            # Timer
            0x12000001: 'Name',
            0x12000002: 'Description',
            0x12000003: 'Interval',
            0x12000004: 'TriggerIn',
            0x12000005: 'TriggerOut',
            0x12000006: 'ActivationTime',
            0x12000007: 'ActivationNumber',
            # Marker
            0x14000001: 'Name',
            0x14000002: 'Description',
            0x14000003: 'TriggerIn',
            0x14000004: 'TriggerOut',
        }

    @cached_property
    def CZ_LSM_LUTTYPE(self):  # TODO: type this
        class CZ_LSM_LUTTYPE(enum.IntEnum):
            NORMAL = 0
            ORIGINAL = 1
            RAMP = 2
            POLYLINE = 3
            SPLINE = 4
            GAMMA = 5

        return CZ_LSM_LUTTYPE

    @cached_property
    def CZ_LSM_SUBBLOCK_TYPE(self):  # TODO: type this
        class CZ_LSM_SUBBLOCK_TYPE(enum.IntEnum):
            END = 0
            GAMMA = 1
            BRIGHTNESS = 2
            CONTRAST = 3
            RAMP = 4
            KNOTS = 5
            PALETTE_12_TO_12 = 6

        return CZ_LSM_SUBBLOCK_TYPE

    @property
    def NIH_IMAGE_HEADER(self):  # TODO: type this
        return [
            ('FileID', 'S8'),
            ('nLines', 'i2'),
            ('PixelsPerLine', 'i2'),
            ('Version', 'i2'),
            ('OldLutMode', 'i2'),
            ('OldnColors', 'i2'),
            ('Colors', 'u1', (3, 32)),
            ('OldColorStart', 'i2'),
            ('ColorWidth', 'i2'),
            ('ExtraColors', 'u2', (6, 3)),
            ('nExtraColors', 'i2'),
            ('ForegroundIndex', 'i2'),
            ('BackgroundIndex', 'i2'),
            ('XScale', 'f8'),
            ('Unused2', 'i2'),
            ('Unused3', 'i2'),
            ('UnitsID', 'i2'),  # NIH_UNITS_TYPE
            ('p1', [('x', 'i2'), ('y', 'i2')]),
            ('p2', [('x', 'i2'), ('y', 'i2')]),
            ('CurveFitType', 'i2'),  # NIH_CURVEFIT_TYPE
            ('nCoefficients', 'i2'),
            ('Coeff', 'f8', 6),
            ('UMsize', 'u1'),
            ('UM', 'S15'),
            ('UnusedBoolean', 'u1'),
            ('BinaryPic', 'b1'),
            ('SliceStart', 'i2'),
            ('SliceEnd', 'i2'),
            ('ScaleMagnification', 'f4'),
            ('nSlices', 'i2'),
            ('SliceSpacing', 'f4'),
            ('CurrentSlice', 'i2'),
            ('FrameInterval', 'f4'),
            ('PixelAspectRatio', 'f4'),
            ('ColorStart', 'i2'),
            ('ColorEnd', 'i2'),
            ('nColors', 'i2'),
            ('Fill1', '3u2'),
            ('Fill2', '3u2'),
            ('Table', 'u1'),  # NIH_COLORTABLE_TYPE
            ('LutMode', 'u1'),  # NIH_LUTMODE_TYPE
            ('InvertedTable', 'b1'),
            ('ZeroClip', 'b1'),
            ('XUnitSize', 'u1'),
            ('XUnit', 'S11'),
            ('StackType', 'i2'),  # NIH_STACKTYPE_TYPE
            # ('UnusedBytes', 'u1', 200)
        ]

    @property
    def NIH_COLORTABLE_TYPE(self) -> tuple[str, ...]:
        return (
            'CustomTable',
            'AppleDefault',
            'Pseudo20',
            'Pseudo32',
            'Rainbow',
            'Fire1',
            'Fire2',
            'Ice',
            'Grays',
            'Spectrum',
        )

    @property
    def NIH_LUTMODE_TYPE(self) -> tuple[str, ...]:
        return (
            'PseudoColor',
            'OldAppleDefault',
            'OldSpectrum',
            'GrayScale',
            'ColorLut',
            'CustomGrayscale',
        )

    @property
    def NIH_CURVEFIT_TYPE(self) -> tuple[str, ...]:
        return (
            'StraightLine',
            'Poly2',
            'Poly3',
            'Poly4',
            'Poly5',
            'ExpoFit',
            'PowerFit',
            'LogFit',
            'RodbardFit',
            'SpareFit1',
            'Uncalibrated',
            'UncalibratedOD',
        )

    @property
    def NIH_UNITS_TYPE(self) -> tuple[str, ...]:
        return (
            'Nanometers',
            'Micrometers',
            'Millimeters',
            'Centimeters',
            'Meters',
            'Kilometers',
            'Inches',
            'Feet',
            'Miles',
            'Pixels',
            'OtherUnits',
        )

    @property
    def TVIPS_HEADER_V1(self) -> list[tuple[str, str]]:
        # TVIPS TemData structure from EMMENU Help file
        return [
            ('Version', 'i4'),
            ('CommentV1', 'S80'),
            ('HighTension', 'i4'),
            ('SphericalAberration', 'i4'),
            ('IlluminationAperture', 'i4'),
            ('Magnification', 'i4'),
            ('PostMagnification', 'i4'),
            ('FocalLength', 'i4'),
            ('Defocus', 'i4'),
            ('Astigmatism', 'i4'),
            ('AstigmatismDirection', 'i4'),
            ('BiprismVoltage', 'i4'),
            ('SpecimenTiltAngle', 'i4'),
            ('SpecimenTiltDirection', 'i4'),
            ('IlluminationTiltDirection', 'i4'),
            ('IlluminationTiltAngle', 'i4'),
            ('ImageMode', 'i4'),
            ('EnergySpread', 'i4'),
            ('ChromaticAberration', 'i4'),
            ('ShutterType', 'i4'),
            ('DefocusSpread', 'i4'),
            ('CcdNumber', 'i4'),
            ('CcdSize', 'i4'),
            ('OffsetXV1', 'i4'),
            ('OffsetYV1', 'i4'),
            ('PhysicalPixelSize', 'i4'),
            ('Binning', 'i4'),
            ('ReadoutSpeed', 'i4'),
            ('GainV1', 'i4'),
            ('SensitivityV1', 'i4'),
            ('ExposureTimeV1', 'i4'),
            ('FlatCorrected', 'i4'),
            ('DeadPxCorrected', 'i4'),
            ('ImageMean', 'i4'),
            ('ImageStd', 'i4'),
            ('DisplacementX', 'i4'),
            ('DisplacementY', 'i4'),
            ('DateV1', 'i4'),
            ('TimeV1', 'i4'),
            ('ImageMin', 'i4'),
            ('ImageMax', 'i4'),
            ('ImageStatisticsQuality', 'i4'),
        ]

    @property
    def TVIPS_HEADER_V2(self) -> list[tuple[str, str]]:
        return [
            ('ImageName', 'V160'),  # utf16
            ('ImageFolder', 'V160'),
            ('ImageSizeX', 'i4'),
            ('ImageSizeY', 'i4'),
            ('ImageSizeZ', 'i4'),
            ('ImageSizeE', 'i4'),
            ('ImageDataType', 'i4'),
            ('Date', 'i4'),
            ('Time', 'i4'),
            ('Comment', 'V1024'),
            ('ImageHistory', 'V1024'),
            ('Scaling', '16f4'),
            ('ImageStatistics', '16c16'),
            ('ImageType', 'i4'),
            ('ImageDisplaType', 'i4'),
            ('PixelSizeX', 'f4'),  # distance between two px in x, [nm]
            ('PixelSizeY', 'f4'),  # distance between two px in y, [nm]
            ('ImageDistanceZ', 'f4'),
            ('ImageDistanceE', 'f4'),
            ('ImageMisc', '32f4'),
            ('TemType', 'V160'),
            ('TemHighTension', 'f4'),
            ('TemAberrations', '32f4'),
            ('TemEnergy', '32f4'),
            ('TemMode', 'i4'),
            ('TemMagnification', 'f4'),
            ('TemMagnificationCorrection', 'f4'),
            ('PostMagnification', 'f4'),
            ('TemStageType', 'i4'),
            ('TemStagePosition', '5f4'),  # x, y, z, a, b
            ('TemImageShift', '2f4'),
            ('TemBeamShift', '2f4'),
            ('TemBeamTilt', '2f4'),
            ('TilingParameters', '7f4'),  # 0: tiling? 1:x 2:y 3: max x
            #                               4: max y 5: overlap x 6: overlap y
            ('TemIllumination', '3f4'),  # 0: spotsize 1: intensity
            ('TemShutter', 'i4'),
            ('TemMisc', '32f4'),
            ('CameraType', 'V160'),
            ('PhysicalPixelSizeX', 'f4'),
            ('PhysicalPixelSizeY', 'f4'),
            ('OffsetX', 'i4'),
            ('OffsetY', 'i4'),
            ('BinningX', 'i4'),
            ('BinningY', 'i4'),
            ('ExposureTime', 'f4'),
            ('Gain', 'f4'),
            ('ReadoutRate', 'f4'),
            ('FlatfieldDescription', 'V160'),
            ('Sensitivity', 'f4'),
            ('Dose', 'f4'),
            ('CamMisc', '32f4'),
            ('FeiMicroscopeInformation', 'V1024'),
            ('FeiSpecimenInformation', 'V1024'),
            ('Magic', 'u4'),
        ]

    @property
    def MM_HEADER(self) -> list[tuple[Any, ...]]:
        # Olympus FluoView MM_Header
        MM_DIMENSION = [
            ('Name', 'S16'),
            ('Size', 'i4'),
            ('Origin', 'f8'),
            ('Resolution', 'f8'),
            ('Unit', 'S64'),
        ]
        return [
            ('HeaderFlag', 'i2'),
            ('ImageType', 'u1'),
            ('ImageName', 'S257'),
            ('OffsetData', 'u4'),
            ('PaletteSize', 'i4'),
            ('OffsetPalette0', 'u4'),
            ('OffsetPalette1', 'u4'),
            ('CommentSize', 'i4'),
            ('OffsetComment', 'u4'),
            ('Dimensions', MM_DIMENSION, 10),
            ('OffsetPosition', 'u4'),
            ('MapType', 'i2'),
            ('MapMin', 'f8'),
            ('MapMax', 'f8'),
            ('MinValue', 'f8'),
            ('MaxValue', 'f8'),
            ('OffsetMap', 'u4'),
            ('Gamma', 'f8'),
            ('Offset', 'f8'),
            ('GrayChannel', MM_DIMENSION),
            ('OffsetThumbnail', 'u4'),
            ('VoiceField', 'i4'),
            ('OffsetVoiceField', 'u4'),
        ]

    @property
    def MM_DIMENSIONS(self) -> dict[str, str]:
        # map FluoView MM_Header.Dimensions to axes characters
        return {
            'X': 'X',
            'Y': 'Y',
            'Z': 'Z',
            'T': 'T',
            'CH': 'C',
            'WAVELENGTH': 'C',
            'TIME': 'T',
            'XY': 'R',
            'EVENT': 'V',
            'EXPOSURE': 'L',
        }

    @property
    def UIC_TAGS(self) -> list[tuple[str, Any]]:
        # map Universal Imaging Corporation MetaMorph internal tag ids to
        # name and type
        from fractions import Fraction

        return [
            ('AutoScale', int),
            ('MinScale', int),
            ('MaxScale', int),
            ('SpatialCalibration', int),
            ('XCalibration', Fraction),
            ('YCalibration', Fraction),
            ('CalibrationUnits', str),
            ('Name', str),
            ('ThreshState', int),
            ('ThreshStateRed', int),
            ('tagid_10', None),  # undefined
            ('ThreshStateGreen', int),
            ('ThreshStateBlue', int),
            ('ThreshStateLo', int),
            ('ThreshStateHi', int),
            ('Zoom', int),
            ('CreateTime', julian_datetime),
            ('LastSavedTime', julian_datetime),
            ('currentBuffer', int),
            ('grayFit', None),
            ('grayPointCount', None),
            ('grayX', Fraction),
            ('grayY', Fraction),
            ('grayMin', Fraction),
            ('grayMax', Fraction),
            ('grayUnitName', str),
            ('StandardLUT', int),
            ('wavelength', int),
            ('StagePosition', '(%i,2,2)u4'),  # N xy positions as fract
            ('CameraChipOffset', '(%i,2,2)u4'),  # N xy offsets as fract
            ('OverlayMask', None),
            ('OverlayCompress', None),
            ('Overlay', None),
            ('SpecialOverlayMask', None),
            ('SpecialOverlayCompress', None),
            ('SpecialOverlay', None),
            ('ImageProperty', read_uic_property),
            ('StageLabel', '%ip'),  # N str
            ('AutoScaleLoInfo', Fraction),
            ('AutoScaleHiInfo', Fraction),
            ('AbsoluteZ', '(%i,2)u4'),  # N fractions
            ('AbsoluteZValid', '(%i,)u4'),  # N long
            ('Gamma', 'I'),  # 'I' uses offset
            ('GammaRed', 'I'),
            ('GammaGreen', 'I'),
            ('GammaBlue', 'I'),
            ('CameraBin', '2I'),
            ('NewLUT', int),
            ('ImagePropertyEx', None),
            ('PlaneProperty', int),
            ('UserLutTable', '(256,3)u1'),
            ('RedAutoScaleInfo', int),
            ('RedAutoScaleLoInfo', Fraction),
            ('RedAutoScaleHiInfo', Fraction),
            ('RedMinScaleInfo', int),
            ('RedMaxScaleInfo', int),
            ('GreenAutoScaleInfo', int),
            ('GreenAutoScaleLoInfo', Fraction),
            ('GreenAutoScaleHiInfo', Fraction),
            ('GreenMinScaleInfo', int),
            ('GreenMaxScaleInfo', int),
            ('BlueAutoScaleInfo', int),
            ('BlueAutoScaleLoInfo', Fraction),
            ('BlueAutoScaleHiInfo', Fraction),
            ('BlueMinScaleInfo', int),
            ('BlueMaxScaleInfo', int),
            # ('OverlayPlaneColor', read_uic_overlay_plane_color),
        ]

    @property
    def PILATUS_HEADER(self) -> dict[str, Any]:
        # PILATUS CBF Header Specification, Version 1.4
        # map key to [value_indices], type
        return {
            'Detector': ([slice(1, None)], str),
            'Pixel_size': ([1, 4], float),
            'Silicon': ([3], float),
            'Exposure_time': ([1], float),
            'Exposure_period': ([1], float),
            'Tau': ([1], float),
            'Count_cutoff': ([1], int),
            'Threshold_setting': ([1], float),
            'Gain_setting': ([1, 2], str),
            'N_excluded_pixels': ([1], int),
            'Excluded_pixels': ([1], str),
            'Flat_field': ([1], str),
            'Trim_file': ([1], str),
            'Image_path': ([1], str),
            # optional
            'Wavelength': ([1], float),
            'Energy_range': ([1, 2], float),
            'Detector_distance': ([1], float),
            'Detector_Voffset': ([1], float),
            'Beam_xy': ([1, 2], float),
            'Flux': ([1], str),
            'Filter_transmission': ([1], float),
            'Start_angle': ([1], float),
            'Angle_increment': ([1], float),
            'Detector_2theta': ([1], float),
            'Polarization': ([1], float),
            'Alpha': ([1], float),
            'Kappa': ([1], float),
            'Phi': ([1], float),
            'Phi_increment': ([1], float),
            'Chi': ([1], float),
            'Chi_increment': ([1], float),
            'Oscillation_axis': ([slice(1, None)], str),
            'N_oscillations': ([1], int),
            'Start_position': ([1], float),
            'Position_increment': ([1], float),
            'Shutter_time': ([1], float),
            'Omega': ([1], float),
            'Omega_increment': ([1], float),
        }

    @cached_property
    def ALLOCATIONGRANULARITY(self) -> int:
        # alignment for writing contiguous data to TIFF
        import mmap

        return mmap.ALLOCATIONGRANULARITY

    @cached_property
    def MAXWORKERS(self) -> int:
        """Default maximum number of threads for de/compressing segments.

        The value of the ``TIFFFILE_NUM_THREADS`` environment variable if set,
        else half the CPU cores up to 32.

        """
        if 'TIFFFILE_NUM_THREADS' in os.environ:
            return max(1, int(os.environ['TIFFFILE_NUM_THREADS']))
        cpu_count: int | None
        try:
            cpu_count = len(
                os.sched_getaffinity(0)  # type: ignore[attr-defined]
            )
        except AttributeError:
            cpu_count = os.cpu_count()
        if cpu_count is None:
            return 1
        return min(32, max(1, cpu_count // 2))

    @cached_property
    def MAXIOWORKERS(self) -> int:
        """Default maximum number of I/O threads for reading file sequences.

        The value of the ``TIFFFILE_NUM_IOTHREADS`` environment variable if
        set, else 4 more than the number of CPU cores up to 32.

        """
        if 'TIFFFILE_NUM_IOTHREADS' in os.environ:
            return max(1, int(os.environ['TIFFFILE_NUM_IOTHREADS']))
        cpu_count: int | None
        try:
            cpu_count = len(
                os.sched_getaffinity(0)  # type: ignore[attr-defined]
            )
        except AttributeError:
            cpu_count = os.cpu_count()
        if cpu_count is None:
            return 5
        return min(32, cpu_count + 4)

    BUFFERSIZE: int = 268435456
    """Default number of bytes to read or encode in one pass (256 MB)."""


TIFF = _TIFF()


def read_tags(
    fh: FileHandle,
    /,
    byteorder: ByteOrder,
    offsetsize: int,
    tagnames: TiffTagRegistry,
    *,
    maxifds: int | None = None,
    customtags: (
        dict[int, Callable[[FileHandle, ByteOrder, int, int, int], Any]] | None
    ) = None,
) -> list[dict[str, Any]]:
    """Read tag values from chain of IFDs.

    Parameters:
        fh:
            Binary file handle to read from.
            The file handle position must be at a valid IFD header.
        byteorder:
            Byte order of TIFF file.
        offsetsize:
            Size of offsets in TIFF file (8 for BigTIFF, else 4).
        tagnames:
            Map of tag codes to names.
            For example, :py:class:`_TIFF.GPS_TAGS` or
            :py:class:`_TIFF.IOP_TAGS`.
        maxifds:
            Maximum number of IFDs to read.
            By default, read the whole IFD chain.
        customtags:
            Mapping of tag codes to functions reading special tag value from
            file.

    Raises:
        TiffFileError: Invalid TIFF structure.

    Notes:
        This implementation does not support 64-bit NDPI files.

    """
    code: int
    dtype: int
    count: int
    valuebytes: bytes
    valueoffset: int

    if offsetsize == 4:
        offsetformat = byteorder + 'I'
        tagnosize = 2
        tagnoformat = byteorder + 'H'
        tagsize = 12
        tagformat1 = byteorder + 'HH'
        tagformat2 = byteorder + 'I4s'
    elif offsetsize == 8:
        offsetformat = byteorder + 'Q'
        tagnosize = 8
        tagnoformat = byteorder + 'Q'
        tagsize = 20
        tagformat1 = byteorder + 'HH'
        tagformat2 = byteorder + 'Q8s'
    else:
        raise ValueError('invalid offset size')

    if customtags is None:
        customtags = {}
    if maxifds is None:
        maxifds = 2**32

    result: list[dict[str, Any]] = []
    unpack = struct.unpack
    offset = fh.tell()
    while len(result) < maxifds:
        # loop over IFDs
        try:
            tagno = unpack(tagnoformat, fh.read(tagnosize))[0]
            if tagno > 4096:
                raise TiffFileError(f'suspicious number of tags {tagno}')
        except Exception as exc:
            logger().error(
                f'<tifffile.read_tags> corrupted tag list @{offset} '
                f'raised {exc!r:.128}'
            )
            break

        tags = {}
        data = fh.read(tagsize * tagno)
        pos = fh.tell()
        index = 0

        for _ in range(tagno):
            code, dtype = unpack(tagformat1, data[index : index + 4])
            count, valuebytes = unpack(
                tagformat2, data[index + 4 : index + tagsize]
            )
            index += tagsize
            name = tagnames.get(code, str(code))
            try:
                valueformat = TIFF.DATA_FORMATS[dtype]
            except KeyError:
                logger().error(f'invalid data type {dtype!r} for tag #{code}')
                continue

            valuesize = count * struct.calcsize(valueformat)
            if valuesize > offsetsize or code in customtags:
                valueoffset = unpack(offsetformat, valuebytes)[0]
                if valueoffset < 8 or valueoffset + valuesize > fh.size:
                    logger().error(
                        f'invalid value offset {valueoffset} for tag #{code}'
                    )
                    continue
                fh.seek(valueoffset)
                if code in customtags:
                    readfunc = customtags[code]
                    value = readfunc(fh, byteorder, dtype, count, offsetsize)
                elif dtype in {1, 2, 7}:
                    # BYTES, ASCII, UNDEFINED
                    value = fh.read(valuesize)
                    if len(value) != valuesize:
                        logger().warning(
                            '<tifffile.read_tags> '
                            f'could not read all values for tag #{code}'
                        )
                elif code in tagnames:
                    fmt = (
                        f'{byteorder}'
                        f'{count * int(valueformat[0])}'
                        f'{valueformat[1]}'
                    )
                    value = unpack(fmt, fh.read(valuesize))
                else:
                    value = read_numpy(fh, byteorder, dtype, count, offsetsize)
            elif dtype in {1, 2, 7}:
                # BYTES, ASCII, UNDEFINED
                value = valuebytes[:valuesize]
            else:
                fmt = (
                    f'{byteorder}'
                    f'{count * int(valueformat[0])}'
                    f'{valueformat[1]}'
                )
                value = unpack(fmt, valuebytes[:valuesize])

            process = (
                code not in customtags
                and code not in TIFF.TAG_TUPLE
                and dtype != 7  # UNDEFINED
            )
            if process and dtype == 2:
                # TIFF ASCII fields can contain multiple strings,
                #   each terminated with a NUL
                value = value.rstrip(b'\x00')
                try:
                    value = value.decode('utf-8').strip()
                except UnicodeDecodeError:
                    try:
                        value = value.decode('cp1252').strip()
                    except UnicodeDecodeError as exc:
                        logger().warning(
                            '<tifffile.read_tags> coercing invalid ASCII to '
                            f'bytes for tag #{code}, due to {exc!r:.128}'
                        )
            else:
                if code in TIFF.TAG_ENUM:
                    t = TIFF.TAG_ENUM[code]
                    try:
                        value = tuple(t(v) for v in value)
                    except ValueError as exc:
                        if code not in {259, 317}:
                            # ignore compression/predictor
                            logger().warning(
                                f'<tifffile.read_tags> tag #{code} '
                                f'raised {exc!r:.128}'
                            )
                if process and len(value) == 1:
                    value = value[0]
            tags[name] = value

        result.append(tags)

        # read offset to next page
        fh.seek(pos)
        offset = unpack(offsetformat, fh.read(offsetsize))[0]
        if offset == 0:
            break
        if offset >= fh.size:
            logger().error(f'<tifffile.read_tags> invalid next page {offset=}')
            break
        fh.seek(offset)

    return result


def read_exif_ifd(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read EXIF tags from file."""
    exif = read_tags(fh, byteorder, offsetsize, TIFF.EXIF_TAGS, maxifds=1)[0]
    for name in ('ExifVersion', 'FlashpixVersion'):
        try:
            exif[name] = bytes2str(exif[name])
        except Exception:
            pass
    if 'UserComment' in exif:
        idcode = exif['UserComment'][:8]
        try:
            if idcode == b'ASCII\x00\x00\x00':
                exif['UserComment'] = bytes2str(exif['UserComment'][8:])
            elif idcode == b'UNICODE\x00':
                exif['UserComment'] = exif['UserComment'][8:].decode('utf-16')
        except Exception:
            pass
    return exif


def read_gps_ifd(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read GPS tags from file."""
    return read_tags(fh, byteorder, offsetsize, TIFF.GPS_TAGS, maxifds=1)[0]


def read_interoperability_ifd(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read Interoperability tags from file."""
    return read_tags(fh, byteorder, offsetsize, TIFF.IOP_TAGS, maxifds=1)[0]


def read_bytes(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> bytes:
    """Read tag data from file."""
    count *= numpy.dtype(
        'B' if dtype == 2 else byteorder + TIFF.DATA_FORMATS[dtype][-1]
    ).itemsize
    data = fh.read(count)
    if len(data) != count:
        logger().warning(
            '<tifffile.read_bytes> '
            f'failed to read {count} bytes, got {len(data)})'
        )
    return data


def read_utf8(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> str:
    """Read unicode tag value from file."""
    return fh.read(count).decode()


def read_numpy(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> NDArray[Any]:
    """Read NumPy array tag value from file."""
    return fh.read_array(
        'b' if dtype == 2 else byteorder + TIFF.DATA_FORMATS[dtype][-1], count
    )


def read_colormap(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> NDArray[Any]:
    """Read ColorMap or TransferFunction tag value from file."""
    cmap = fh.read_array(byteorder + TIFF.DATA_FORMATS[dtype][-1], count)
    if count % 3 == 0:
        cmap.shape = (3, -1)
    return cmap


def read_json(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> Any:
    """Read JSON tag value from file."""
    data = fh.read(count)
    try:
        return json.loads(bytes2str(data, 'utf-8'))
    except ValueError as exc:
        logger().warning(f'<tifffile.read_json> raised {exc!r:.128}')
    return None


def read_mm_header(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read FluoView mm_header tag value from file."""
    meta = recarray2dict(
        fh.read_record(numpy.dtype(TIFF.MM_HEADER), byteorder=byteorder)
    )
    meta['Dimensions'] = [
        (bytes2str(d[0]).strip(), d[1], d[2], d[3], bytes2str(d[4]).strip())
        for d in meta['Dimensions']
    ]
    d = meta['GrayChannel']
    meta['GrayChannel'] = (
        bytes2str(d[0]).strip(),
        d[1],
        d[2],
        d[3],
        bytes2str(d[4]).strip(),
    )
    return meta


def read_mm_stamp(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> NDArray[Any]:
    """Read FluoView mm_stamp tag value from file."""
    return fh.read_array(byteorder + 'f8', 8)


def read_uic1tag(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
    planecount: int = 0,
) -> dict[str, Any]:
    """Read MetaMorph STK UIC1Tag value from file.

    Return empty dictionary if planecount is unknown.

    """
    if dtype not in {4, 5} or byteorder != '<':
        raise ValueError(f'invalid UIC1Tag {byteorder}{dtype}')
    result = {}
    if dtype == 5:
        # pre MetaMorph 2.5 (not tested)
        values = fh.read_array('<u4', 2 * count).reshape(count, 2)
        result = {'ZDistance': values[:, 0] / values[:, 1]}
    else:
        for _ in range(count):
            tagid = struct.unpack('<I', fh.read(4))[0]
            if planecount > 1 and tagid in {28, 29, 37, 40, 41}:
                # silently skip unexpected tags
                fh.read(4)
                continue
            name, value = read_uic_tag(fh, tagid, planecount, True)
            if name == 'PlaneProperty':
                pos = fh.tell()
                fh.seek(value + 4)
                result.setdefault(name, []).append(read_uic_property(fh))
                fh.seek(pos)
            else:
                result[name] = value
    return result


def read_uic2tag(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, NDArray[Any]]:
    """Read MetaMorph STK UIC2Tag value from file."""
    if dtype != 5 or byteorder != '<':
        raise ValueError('invalid UIC2Tag')
    values = fh.read_array('<u4', 6 * count).reshape(count, 6)
    return {
        'ZDistance': values[:, 0] / values[:, 1],
        'DateCreated': values[:, 2],  # julian days
        'TimeCreated': values[:, 3],  # milliseconds
        'DateModified': values[:, 4],  # julian days
        'TimeModified': values[:, 5],  # milliseconds
    }


def read_uic3tag(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, NDArray[Any]]:
    """Read MetaMorph STK UIC3Tag value from file."""
    if dtype != 5 or byteorder != '<':
        raise ValueError('invalid UIC3Tag')
    values = fh.read_array('<u4', 2 * count).reshape(count, 2)
    return {'Wavelengths': values[:, 0] / values[:, 1]}


def read_uic4tag(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, NDArray[Any]]:
    """Read MetaMorph STK UIC4Tag value from file."""
    if dtype != 4 or byteorder != '<':
        raise ValueError('invalid UIC4Tag')
    result = {}
    while True:
        tagid: int = struct.unpack('<H', fh.read(2))[0]
        if tagid == 0:
            break
        name, value = read_uic_tag(fh, tagid, count, False)
        result[name] = value
    return result


def read_uic_tag(
    fh: FileHandle, tagid: int, planecount: int, offset: bool, /
) -> tuple[str, Any]:
    """Read single UIC tag value from file and return tag name and value.

    UIC1Tags use an offset.

    """

    def read_int() -> int:
        return int(struct.unpack('<I', fh.read(4))[0])

    def read_int2() -> tuple[int, int]:
        value = struct.unpack('<2I', fh.read(8))
        return int(value[0]), (value[1])

    try:
        name, dtype = TIFF.UIC_TAGS[tagid]
    except IndexError:
        # unknown tag
        return f'_TagId{tagid}', read_int()

    Fraction = TIFF.UIC_TAGS[4][1]

    if offset:
        pos = fh.tell()
        if dtype not in {int, None}:
            off = read_int()
            if off < 8:
                # undocumented cases, or invalid offset
                if dtype is str:
                    return name, ''
                if tagid == 41:  # AbsoluteZValid
                    return name, off
                logger().warning(
                    '<tifffile.read_uic_tag> '
                    f'invalid offset for tag {name!r} @{off}'
                )
                return name, off
            fh.seek(off)

    value: Any

    if dtype is None:
        # skip
        name = '_' + name
        value = read_int()
    elif dtype is int:
        # int
        value = read_int()
    elif dtype is Fraction:
        # fraction
        value = read_int2()
        value = value[0] / value[1]
    elif dtype is julian_datetime:
        # datetime
        value = read_int2()
        try:
            value = julian_datetime(*value)
        except Exception as exc:
            value = None
            logger().warning(
                f'<tifffile.read_uic_tag> reading {name} raised {exc!r:.128}'
            )
    elif dtype is read_uic_property:
        # ImagePropertyEx
        value = read_uic_property(fh)
    elif dtype is str:
        # pascal string
        size = read_int()
        if 0 <= size < 2**10:
            value = struct.unpack(f'{size}s', fh.read(size))[0][:-1]
            value = bytes2str(value)
        elif offset:
            value = ''
            logger().warning(
                f'<tifffile.read_uic_tag> invalid string in tag {name!r}'
            )
        else:
            raise ValueError(f'invalid string size {size}')
    elif planecount == 0:
        value = None
    elif dtype == '%ip':
        # sequence of pascal strings
        value = []
        for _ in range(planecount):
            size = read_int()
            if 0 <= size < 2**10:
                string = struct.unpack(f'{size}s', fh.read(size))[0][:-1]
                value.append(bytes2str(string))
            elif offset:
                logger().warning(
                    f'<tifffile.read_uic_tag> invalid string in tag {name!r}'
                )
            else:
                raise ValueError(f'invalid string size: {size}')
    else:
        # struct or numpy type
        dtype = '<' + dtype
        if '%i' in dtype:
            dtype = dtype % planecount
        if '(' in dtype:
            # numpy type
            value = fh.read_array(dtype, 1)[0]
            if value.shape[-1] == 2:
                # assume fractions
                value = value[..., 0] / value[..., 1]
        else:
            # struct format
            value = struct.unpack(dtype, fh.read(struct.calcsize(dtype)))
            if len(value) == 1:
                value = value[0]

    if offset:
        fh.seek(pos + 4)

    return name, value


def read_uic_property(fh: FileHandle, /) -> dict[str, Any]:
    """Read UIC ImagePropertyEx or PlaneProperty tag from file."""
    size = struct.unpack('B', fh.read(1))[0]
    name = bytes2str(struct.unpack(f'{size}s', fh.read(size))[0])
    flags, prop = struct.unpack('<IB', fh.read(5))
    if prop == 1:
        value = struct.unpack('II', fh.read(8))
        value = value[0] / value[1]
    else:
        size = struct.unpack('B', fh.read(1))[0]
        value = bytes2str(
            struct.unpack(f'{size}s', fh.read(size))[0]
        )  # type: ignore[assignment]
    return {'name': name, 'flags': flags, 'value': value}


def read_cz_lsminfo(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read CZ_LSMINFO tag value from file."""
    if byteorder != '<':
        raise ValueError('invalid CZ_LSMINFO structure')
    magic_number, structure_size = struct.unpack('<II', fh.read(8))
    if magic_number not in {50350412, 67127628}:
        raise ValueError('invalid CZ_LSMINFO structure')
    fh.seek(-8, os.SEEK_CUR)
    CZ_LSMINFO = TIFF.CZ_LSMINFO

    if structure_size < numpy.dtype(CZ_LSMINFO).itemsize:
        # adjust structure according to structure_size
        lsminfo: list[tuple[str, str]] = []
        size = 0
        for name, typestr in CZ_LSMINFO:
            size += numpy.dtype(typestr).itemsize
            if size > structure_size:
                break
            lsminfo.append((name, typestr))
    else:
        lsminfo = CZ_LSMINFO

    result = recarray2dict(
        fh.read_record(numpy.dtype(lsminfo), byteorder=byteorder)
    )

    # read LSM info subrecords at offsets
    for name, reader in TIFF.CZ_LSMINFO_READERS.items():
        if reader is None:
            continue
        offset = result.get('Offset' + name, 0)
        if offset < 8:
            continue
        fh.seek(offset)
        try:
            result[name] = reader(fh)
        except ValueError:
            pass
    return result


def read_lsm_channeldatatypes(fh: FileHandle, /) -> NDArray[Any]:
    """Read LSM channel data type from file."""
    size = struct.unpack('<I', fh.read(4))[0]
    return fh.read_array('<u4', count=size)


def read_lsm_channelwavelength(fh: FileHandle, /) -> NDArray[Any]:
    """Read LSM channel wavelength ranges from file."""
    size = struct.unpack('<i', fh.read(4))[0]
    return fh.read_array('<2f8', count=size)


def read_lsm_positions(fh: FileHandle, /) -> NDArray[Any]:
    """Read LSM positions from file."""
    size = struct.unpack('<I', fh.read(4))[0]
    return fh.read_array('<3f8', count=size)


def read_lsm_timestamps(fh: FileHandle, /) -> NDArray[Any]:
    """Read LSM time stamps from file."""
    size, count = struct.unpack('<ii', fh.read(8))
    if size != (8 + 8 * count):
        logger().warning(
            '<tifffile.read_lsm_timestamps> invalid LSM TimeStamps block'
        )
        return numpy.empty((0,), '<f8')
    # return struct.unpack(f'<{count}d', fh.read(8 * count))
    return fh.read_array('<f8', count=count)


def read_lsm_eventlist(fh: FileHandle, /) -> list[tuple[float, int, str]]:
    """Read LSM events from file and return as list of (time, type, text)."""
    count = struct.unpack('<II', fh.read(8))[1]
    events = []
    while count > 0:
        esize, etime, etype = struct.unpack('<IdI', fh.read(16))
        etext = bytes2str(fh.read(esize - 16))
        events.append((etime, etype, etext))
        count -= 1
    return events


def read_lsm_channelcolors(fh: FileHandle, /) -> dict[str, Any]:
    """Read LSM ChannelColors structure from file."""
    result = {'Mono': False, 'Colors': [], 'ColorNames': []}
    pos = fh.tell()
    (size, ncolors, nnames, coffset, noffset, mono) = struct.unpack(
        '<IIIIII', fh.read(24)
    )
    if ncolors != nnames:
        logger().warning(
            '<tifffile.read_lsm_channelcolors> '
            'invalid LSM ChannelColors structure'
        )
        return result
    result['Mono'] = bool(mono)
    # Colors
    fh.seek(pos + coffset)
    colors = fh.read_array(numpy.uint8, count=ncolors * 4)
    colors = colors.reshape((ncolors, 4))
    result['Colors'] = colors.tolist()
    # ColorNames
    fh.seek(pos + noffset)
    buffer = fh.read(size - noffset)
    names = []
    while len(buffer) > 4:
        size = struct.unpack('<I', buffer[:4])[0]
        names.append(bytes2str(buffer[4 : 3 + size]))
        buffer = buffer[4 + size :]
    result['ColorNames'] = names
    return result


def read_lsm_lookuptable(fh: FileHandle, /) -> dict[str, Any]:
    """Read LSM lookup tables from file."""
    result: dict[str, Any] = {}
    (
        size,
        nsubblocks,
        nchannels,
        luttype,
        advanced,
        currentchannel,
    ) = struct.unpack('<iiiiii', fh.read(24))
    if size < 60:
        logger().warning(
            '<tifffile.read_lsm_lookuptable> '
            'invalid LSM LookupTables structure'
        )
        return result
    fh.read(9 * 4)  # reserved
    result['LutType'] = TIFF.CZ_LSM_LUTTYPE(luttype)
    result['Advanced'] = advanced
    result['NumberChannels'] = nchannels
    result['CurrentChannel'] = currentchannel
    result['SubBlocks'] = subblocks = []
    for _ in range(nsubblocks):
        sbtype = struct.unpack('<i', fh.read(4))[0]
        if sbtype <= 0:
            break
        size = struct.unpack('<i', fh.read(4))[0] - 8
        if sbtype == 1:
            data = fh.read_array('<f8', count=nchannels)
        elif sbtype == 2:
            data = fh.read_array('<f8', count=nchannels)
        elif sbtype == 3:
            data = fh.read_array('<f8', count=nchannels)
        elif sbtype == 4:
            # the data type is wrongly documented as f8
            data = fh.read_array('<i4', count=nchannels * 4)
            data = data.reshape((-1, 2, 2))
        elif sbtype == 5:
            # the data type is wrongly documented as f8
            nknots = struct.unpack('<i', fh.read(4))[0]  # undocumented
            data = fh.read_array('<i4', count=nchannels * nknots * 2)
            data = data.reshape((nchannels, nknots, 2))
        elif sbtype == 6:
            data = fh.read_array('<i2', count=nchannels * 4096)
            data = data.reshape((-1, 4096))
        else:
            logger().warning(
                '<tifffile.read_lsm_lookuptable> '
                f'invalid LSM SubBlock type {sbtype}'
            )
            break
        subblocks.append(
            {'Type': TIFF.CZ_LSM_SUBBLOCK_TYPE(sbtype), 'Data': data}
        )
    return result


def read_lsm_scaninfo(fh: FileHandle, /) -> dict[str, Any]:
    """Read LSM ScanInfo structure from file."""
    value: Any
    block: dict[str, Any] = {}
    blocks = [block]
    unpack = struct.unpack
    if struct.unpack('<I', fh.read(4))[0] != 0x10000000:
        # not a Recording sub block
        logger().warning(
            '<tifffile.read_lsm_scaninfo> invalid LSM ScanInfo structure'
        )
        return block
    fh.read(8)
    while True:
        entry, dtype, size = unpack('<III', fh.read(12))
        if dtype == 2:
            # ascii
            value = bytes2str(fh.read(size))
        elif dtype == 4:
            # long
            value = unpack('<i', fh.read(4))[0]
        elif dtype == 5:
            # rational
            value = unpack('<d', fh.read(8))[0]
        else:
            value = 0
        if entry in TIFF.CZ_LSMINFO_SCANINFO_ARRAYS:
            blocks.append(block)
            name = TIFF.CZ_LSMINFO_SCANINFO_ARRAYS[entry]
            newlist: list[dict[str, Any]] = []
            block[name] = newlist
            # TODO: fix types
            block = newlist  # type: ignore[assignment]
        elif entry in TIFF.CZ_LSMINFO_SCANINFO_STRUCTS:
            blocks.append(block)
            newdict: dict[str, Any] = {}
            # TODO: fix types
            block.append(newdict)  # type: ignore[attr-defined]
            block = newdict
        elif entry in TIFF.CZ_LSMINFO_SCANINFO_ATTRIBUTES:
            block[TIFF.CZ_LSMINFO_SCANINFO_ATTRIBUTES[entry]] = value
        elif entry == 0xFFFFFFFF:
            # end sub block
            block = blocks.pop()
        else:
            # unknown entry
            block[f'Entry0x{entry:x}'] = value
        if not blocks:
            break
    return block


def read_sis(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read OlympusSIS structure from file.

    No specification is available. Only few fields are known.

    """
    result: dict[str, Any] = {}

    (magic, minute, hour, day, month, year, name, tagcount) = struct.unpack(
        '<4s6xhhhhh6x32sh', fh.read(60)
    )

    if magic != b'SIS0':
        raise ValueError('invalid OlympusSIS structure')

    result['name'] = bytes2str(name)
    try:
        result['datetime'] = DateTime(
            1900 + year, month + 1, day, hour, minute
        )
    except ValueError:
        pass

    data = fh.read(8 * tagcount)
    for i in range(0, tagcount * 8, 8):
        tagtype, count, offset = struct.unpack('<hhI', data[i : i + 8])
        fh.seek(offset)
        if tagtype == 1:
            # general data
            (lenexp, xcal, ycal, mag, camname, pictype) = struct.unpack(
                '<10xhdd8xd2x34s32s', fh.read(112)  # 220
            )
            m = math.pow(10, lenexp)
            result['pixelsizex'] = xcal * m
            result['pixelsizey'] = ycal * m
            result['magnification'] = mag
            result['cameraname'] = bytes2str(camname)
            result['picturetype'] = bytes2str(pictype)
        elif tagtype == 10:
            # channel data
            continue
            # TODO: does not seem to work?
            # (length, _, exptime, emv, _, camname, _, mictype,
            #  ) = struct.unpack('<h22sId4s32s48s32s', fh.read(152))  # 720
            # result['exposuretime'] = exptime
            # result['emvoltage'] = emv
            # result['cameraname2'] = bytes2str(camname)
            # result['microscopename'] = bytes2str(mictype)

    return result


def read_sis_ini(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read OlympusSIS INI string from file."""
    try:
        return olympusini_metadata(bytes2str(fh.read(count)))
    except Exception as exc:
        logger().warning(f'<tifffile.olympusini_metadata> raised {exc!r:.128}')
        return {}


def read_tvips_header(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read TVIPS EM-MENU headers from file."""
    result: dict[str, Any] = {}
    header_v1 = TIFF.TVIPS_HEADER_V1
    header = fh.read_record(numpy.dtype(header_v1), byteorder=byteorder)
    for name, typestr in header_v1:
        result[name] = header[name].tolist()
    if header['Version'] == 2:
        header_v2 = TIFF.TVIPS_HEADER_V2
        header = fh.read_record(numpy.dtype(header_v2), byteorder=byteorder)
        if header['Magic'] != 0xAAAAAAAA:
            logger().warning(
                '<tifffile.read_tvips_header> invalid TVIPS v2 magic number'
            )
            return {}
        # decode utf16 strings
        for name, typestr in header_v2:
            if typestr.startswith('V'):
                result[name] = bytes2str(
                    header[name].tobytes(), 'utf-16', 'ignore'
                )
            else:
                result[name] = header[name].tolist()
        # convert nm to m
        for axis in 'XY':
            header['PhysicalPixelSize' + axis] /= 1e9
            header['PixelSize' + axis] /= 1e9
    elif header.version != 1:
        logger().warning(
            '<tifffile.read_tvips_header> unknown TVIPS header version'
        )
        return {}
    return result


def read_fei_metadata(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read FEI SFEG/HELIOS headers from file."""
    result: dict[str, Any] = {}
    section: dict[str, Any] = {}
    data = bytes2str(fh.read(count))
    for line in data.splitlines():
        line = line.strip()
        if line.startswith('['):
            section = {}
            result[line[1:-1]] = section
            continue
        try:
            key, value = line.split('=')
        except ValueError:
            continue
        section[key] = astype(value)
    return result


def read_cz_sem(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read Zeiss SEM tag from file.

    See https://sourceforge.net/p/gwyddion/mailman/message/29275000/ for
    unnamed values.

    """
    result: dict[str, Any] = {'': ()}
    value: Any
    key = None
    data = bytes2str(fh.read(count))
    for line in data.splitlines():
        if line.isupper():
            key = line.lower()
        elif key:
            try:
                name, value = line.split('=')
            except ValueError:
                try:
                    name, value = line.split(':', 1)
                except Exception:
                    continue
            value = value.strip()
            unit = ''
            try:
                v, u = value.split()
                number = astype(v, (int, float))
                if number != v:
                    value = number
                    unit = u
            except Exception:
                number = astype(value, (int, float))
                if number != value:
                    value = number
                if value in {'No', 'Off'}:
                    value = False
                elif value in {'Yes', 'On'}:
                    value = True
            result[key] = (name.strip(), value)
            if unit:
                result[key] += (unit,)
            key = None
        else:
            result[''] += (astype(line, (int, float)),)
    return result


def read_nih_image_header(
    fh: FileHandle,
    byteorder: ByteOrder,
    dtype: int,
    count: int,
    offsetsize: int,
    /,
) -> dict[str, Any]:
    """Read NIH_IMAGE_HEADER tag value from file."""
    arr = fh.read_record(TIFF.NIH_IMAGE_HEADER, byteorder=byteorder)
    arr = arr.view(arr.dtype.newbyteorder(byteorder))
    result = recarray2dict(arr)
    result['XUnit'] = result['XUnit'][: result['XUnitSize']]
    result['UM'] = result['UM'][: result['UMsize']]
    return result


def read_scanimage_metadata(
    fh: FileHandle, /
) -> tuple[dict[str, Any], dict[str, Any], int]:
    """Read ScanImage BigTIFF v3 or v4 static and ROI metadata from file.

    The settings can be used to read image and metadata without parsing
    the TIFF file.

    Frame data and ROI groups can alternatively be obtained from the Software
    and Artist tags of any TIFF page.

    Parameters:
        fh: Binary file handle to read from.

    Returns:
        - Non-varying frame data, parsed with :py:func:`matlabstr2py`.
        - ROI group data, parsed from JSON.
        - Version of metadata (3 or 4).

    Raises:
        ValueError: File does not contain valid ScanImage metadata.

    """
    fh.seek(0)
    try:
        byteorder, version = struct.unpack('<2sH', fh.read(4))
        if byteorder != b'II' or version != 43:
            raise ValueError('not a BigTIFF file')
        fh.seek(16)
        magic, version, size0, size1 = struct.unpack('<IIII', fh.read(16))
        if magic != 117637889 or version not in {3, 4}:
            raise ValueError(
                f'invalid magic {magic} or version {version} number'
            )
    except UnicodeDecodeError as exc:
        raise ValueError('file must be opened in binary mode') from exc
    except Exception as exc:
        raise ValueError('not a ScanImage BigTIFF v3 or v4 file') from exc

    frame_data = matlabstr2py(bytes2str(fh.read(size0)[:-1]))
    roi_data = read_json(fh, '<', 0, size1, 0) if size1 > 1 else {}
    return frame_data, roi_data, version


def read_micromanager_metadata(
    fh: FileHandle | IO[bytes], /, keys: Container[str] | None = None
) -> dict[str, Any]:
    """Return Micro-Manager non-TIFF settings from file.

    The settings can be used to read image data without parsing any TIFF
    structures.

    Parameters:
        fh: Open file handle to Micro-Manager TIFF file.
        keys: Name of keys to return in result.

    Returns:
        Micro-Manager non-TIFF settings, which may contain the following keys:

        - 'MajorVersion' (str)
        - 'MinorVersion' (str)
        - 'Summary' (dict):
          Specifies the dataset, such as shape, dimensions, and coordinates.
        - 'IndexMap' (numpy.ndarray):
          (channel, slice, frame, position, ifd_offset) indices of all frames.
        - 'DisplaySettings' (list[dict]):
          Image display settings such as channel contrast and colors.
        - 'Comments' (dict):
          User comments.

    Notes:
        Summary metadata are the same for all files in a dataset.
        DisplaySettings metadata are frequently corrupted, and Comments are
        often empty.
        The Summary and IndexMap metadata are stored at the beginning of
        the file, while DisplaySettings and Comments are towards the end.
        Excluding DisplaySettings and Comments from the results may
        significantly speed up reading metadata of interest.

    References:
        - https://micro-manager.org/Micro-Manager_File_Formats
        - https://github.com/micro-manager/NDTiffStorage

    """
    if keys is None:
        keys = {'Summary', 'IndexMap', 'DisplaySettings', 'Comments'}
    fh.seek(0)
    try:
        byteorder = {b'II': '<', b'MM': '>'}[fh.read(2)]
        fh.seek(8)
        (
            index_header,
            index_offset,
        ) = struct.unpack(byteorder + 'II', fh.read(8))
    except Exception as exc:
        raise ValueError('not a Micro-Manager TIFF file') from exc

    result = {}
    if index_header == 483729:
        # NDTiff >= v2
        result['MajorVersion'] = index_offset
        try:
            (
                summary_header,
                summary_length,
            ) = struct.unpack(byteorder + 'II', fh.read(8))
            if summary_header != 2355492:
                # NDTiff v3
                result['MinorVersion'] = summary_header
                if summary_length != 2355492:
                    raise ValueError(
                        f'invalid summary_length {summary_length}'
                    )
                summary_length = struct.unpack(byteorder + 'I', fh.read(4))[0]
            if 'Summary' in keys:
                data = fh.read(summary_length)
                if len(data) != summary_length:
                    raise ValueError('not enough data')
                result['Summary'] = json.loads(bytes2str(data, 'utf-8'))
        except Exception as exc:
            logger().warning(
                '<tifffile.read_micromanager_metadata> '
                f'failed to read NDTiffv{index_offset} summary settings, '
                f'raised {exc!r:.128}'
            )
        return result

    # Micro-Manager multipage TIFF or NDTiff v1
    try:
        (
            display_header,
            display_offset,
            comments_header,
            comments_offset,
            summary_header,
            summary_length,
        ) = struct.unpack(byteorder + 'IIIIII', fh.read(24))
    except Exception as exc:
        logger().warning(
            '<tifffile.read_micromanager_metadata> failed to read header, '
            f'raised {exc!r:.128}'
        )

    if 'Summary' in keys:
        try:
            if summary_header != 2355492:
                raise ValueError(f'invalid summary_header {summary_header}')
            data = fh.read(summary_length)
            if len(data) != summary_length:
                raise ValueError('not enough data')
            result['Summary'] = json.loads(bytes2str(data, 'utf-8'))
        except Exception as exc:
            logger().warning(
                '<tifffile.read_micromanager_metadata> '
                f'failed to read summary settings, raised {exc!r:.128}'
            )

    if 'IndexMap' in keys:
        try:
            if index_header != 54773648:
                raise ValueError(f'invalid index_header {index_header}')
            fh.seek(index_offset)
            header, count = struct.unpack(byteorder + 'II', fh.read(8))
            if header != 3453623:
                raise ValueError('invalid header')
            data = fh.read(count * 20)
            result['IndexMap'] = numpy.frombuffer(
                data, byteorder + 'u4', count * 5
            ).reshape(-1, 5)
        except Exception as exc:
            logger().warning(
                '<tifffile.read_micromanager_metadata> '
                f'failed to read index map, raised {exc!r:.128}'
            )

    if 'DisplaySettings' in keys:
        try:
            if display_header != 483765892:
                raise ValueError(f'invalid display_header {display_header}')
            fh.seek(display_offset)
            header, count = struct.unpack(byteorder + 'II', fh.read(8))
            if header != 347834724:
                # display_offset might be wrapped at 4 GB
                fh.seek(display_offset + 2**32)
                header, count = struct.unpack(byteorder + 'II', fh.read(8))
                if header != 347834724:
                    raise ValueError('invalid display header')
            data = fh.read(count)
            if len(data) != count:
                raise ValueError('not enough data')
            result['DisplaySettings'] = json.loads(bytes2str(data, 'utf-8'))
        except json.decoder.JSONDecodeError:
            pass  # ignore frequent truncated JSON data
        except Exception as exc:
            logger().warning(
                '<tifffile.read_micromanager_metadata> '
                f'failed to read display settings, raised {exc!r:.128}'
            )

    result['MajorVersion'] = 0
    try:
        if comments_header == 99384722:
            # Micro-Manager multipage TIFF
            if 'Comments' in keys:
                fh.seek(comments_offset)
                header, count = struct.unpack(byteorder + 'II', fh.read(8))
                if header != 84720485:
                    # comments_offset might be wrapped at 4 GB
                    fh.seek(comments_offset + 2**32)
                    header, count = struct.unpack(byteorder + 'II', fh.read(8))
                    if header != 84720485:
                        raise ValueError('invalid comments header')
                data = fh.read(count)
                if len(data) != count:
                    raise ValueError('not enough data')
                result['Comments'] = json.loads(bytes2str(data, 'utf-8'))
        elif comments_header == 483729:
            # NDTiff v1
            result['MajorVersion'] = comments_offset
        elif comments_header == 0 and comments_offset == 0:
            pass
        elif 'Comments' in keys:
            raise ValueError(f'invalid comments_header {comments_header}')
    except Exception as exc:
        logger().warning(
            '<tifffile.read_micromanager_metadata> failed to read comments, '
            f'raised {exc!r:.128}'
        )

    return result


def read_ndtiff_index(
    file: str | os.PathLike[Any], /
) -> Iterator[
    tuple[dict[str, int | str], str, int, int, int, int, int, int, int, int]
]:
    """Return iterator over fields in Micro-Manager NDTiff.index file.

    Parameters:
        file: Path of NDTiff.index file.

    Yields:
        Fields in NDTiff.index file:

        - axes_dict: Axes indices of frame in image.
        - filename: Name of file containing frame and metadata.
        - dataoffset: Offset of frame data in file.
        - width: Width of frame.
        - height: Height of frame.
        - pixeltype: Pixel type.
          0: 8-bit monochrome;
          1: 16-bit monochrome;
          2: 8-bit RGB;
          3: 10-bit monochrome;
          4: 12-bit monochrome;
          5: 14-bit monochrome;
          6: 11-bit monochrome.
        - compression: Pixel compression. 0: Uncompressed.
        - metaoffset: Offset of JSON metadata in file.
        - metabytecount: Length of metadata.
        - metacompression: Metadata compression. 0: Uncompressed.

    """
    with open(file, 'rb') as fh:
        while True:
            b = fh.read(4)
            if len(b) != 4:
                break
            k = struct.unpack('<i', b)[0]
            axes_dict = json.loads(fh.read(k))
            n = struct.unpack('<i', fh.read(4))[0]
            filename = fh.read(n).decode()
            (
                dataoffset,
                width,
                height,
                pixeltype,
                compression,
                metaoffset,
                metabytecount,
                metacompression,
            ) = struct.unpack('<IiiiiIii', fh.read(32))
            yield (
                axes_dict,
                filename,
                dataoffset,
                width,
                height,
                pixeltype,
                compression,
                metaoffset,
                metabytecount,
                metacompression,
            )


def read_gdal_structural_metadata(
    fh: FileHandle | IO[bytes], /
) -> dict[str, str] | None:
    """Read non-TIFF GDAL structural metadata from file.

    Return None if the file does not contain valid GDAL structural metadata.
    The metadata can be used to optimize reading image data from a COG file.

    """
    fh.seek(0)
    try:
        if fh.read(2) not in {b'II', b'MM'}:
            raise ValueError('not a TIFF file')
        fh.seek({b'*': 8, b'+': 16}[fh.read(1)])
        header = fh.read(43).decode()
        if header[:30] != 'GDAL_STRUCTURAL_METADATA_SIZE=':
            return None
        size = int(header[30:36])
        lines = fh.read(size).decode()
    except Exception:
        return None

    result: dict[str, Any] = {}
    try:
        for line in lines.splitlines():
            if '=' in line:
                key, value = line.split('=', 1)
                result[key.strip()] = value.strip()
    except Exception as exc:
        logger().warning(
            f'<tifffile.read_gdal_structural_metadata> raised {exc!r:.128}'
        )
        return None
    return result


def read_metaseries_catalog(fh: FileHandle | IO[bytes], /) -> None:
    """Read MetaSeries non-TIFF hint catalog from file.

    Raise ValueError if the file does not contain a valid hint catalog.

    """
    # TODO: implement read_metaseries_catalog
    raise NotImplementedError


def imagej_metadata_tag(
    metadata: dict[str, Any], byteorder: ByteOrder, /
) -> tuple[
    tuple[int, int, int, bytes, bool], tuple[int, int, int, bytes, bool]
]:
    """Return IJMetadata and IJMetadataByteCounts tags from metadata dict.

    Parameters:
        metadata:
            May contain the following keys and values:

            'Info' (str):
                Human-readable information as string.
            'Labels' (Sequence[str]):
                Human-readable label for each image.
            'Ranges' (Sequence[float]):
                Lower and upper values for each channel.
            'LUTs' (list[numpy.ndarray[(3, 256), 'uint8']]):
                Color palettes for each channel.
            'Plot' (bytes):
                Undocumented ImageJ internal format.
            'ROI', 'Overlays' (bytes):
                Undocumented ImageJ internal region of interest and overlay
                format. Can be created with the
                `roifile <https://pypi.org/project/roifile/>`_ package.
            'Properties' (dict[str, str]):
                Map of key, value items.

        byteorder:
            Byte order of TIFF file.

    Returns:
        IJMetadata and IJMetadataByteCounts tags in :py:meth:`TiffWriter.write`
        `extratags` format.

    """
    if not metadata:
        return ()  # type: ignore[return-value]
    header_list = [{'>': b'IJIJ', '<': b'JIJI'}[byteorder]]
    bytecount_list = [0]
    body_list = []

    def _string(data: str, byteorder: ByteOrder, /) -> bytes:
        return data.encode('utf-16' + {'>': 'be', '<': 'le'}[byteorder])

    def _doubles(data: Sequence[float], byteorder: ByteOrder, /) -> bytes:
        return struct.pack(f'{byteorder}{len(data)}d', *data)

    def _ndarray(data: NDArray[Any], byteorder: ByteOrder, /) -> bytes:
        return data.tobytes()

    def _bytes(data: bytes, byteorder: ByteOrder, /) -> bytes:
        return data

    metadata_types: tuple[
        tuple[str, bytes, Callable[[Any, ByteOrder], bytes]], ...
    ] = (
        ('Info', b'info', _string),
        ('Labels', b'labl', _string),
        ('Ranges', b'rang', _doubles),
        ('LUTs', b'luts', _ndarray),
        ('Plot', b'plot', _bytes),
        ('ROI', b'roi ', _bytes),
        ('Overlays', b'over', _bytes),
        ('Properties', b'prop', _string),
    )

    for key, mtype, func in metadata_types:
        if key.lower() in metadata:
            key = key.lower()
        elif key not in metadata:
            continue
        if byteorder == '<':
            mtype = mtype[::-1]
        values = metadata[key]
        if isinstance(values, dict):
            values = [str(i) for item in values.items() for i in item]
            count = len(values)
        elif isinstance(values, list):
            count = len(values)
        else:
            values = [values]
            count = 1
        header_list.append(mtype + struct.pack(byteorder + 'I', count))
        for value in values:
            data = func(value, byteorder)
            body_list.append(data)
            bytecount_list.append(len(data))

    if not body_list:
        return ()  # type: ignore[return-value]
    body = b''.join(body_list)
    header = b''.join(header_list)
    data = header + body
    bytecount_list[0] = len(header)
    bytecounts = struct.pack(
        byteorder + ('I' * len(bytecount_list)), *bytecount_list
    )
    return (
        (50839, 1, len(data), data, True),
        (50838, 4, len(bytecounts) // 4, bytecounts, True),
    )


def imagej_metadata(
    data: bytes, bytecounts: Sequence[int], byteorder: ByteOrder, /
) -> dict[str, Any]:
    """Return IJMetadata tag value.

    Parameters:
        bytes:
            Encoded value of IJMetadata tag.
        bytecounts:
            Value of IJMetadataByteCounts tag.
        byteorder:
            Byte order of TIFF file.

    Returns:
        Metadata dict with optional items:

            'Info' (str):
                Human-readable information as string.
                Some formats, such as OIF or ScanImage, can be parsed into
                dicts with :py:func:`matlabstr2py` or the
                `oiffile.SettingsFile()` function of the
                `oiffile <https://pypi.org/project/oiffile/>`_  package.
            'Labels' (Sequence[str]):
                Human-readable labels for each channel.
            'Ranges' (Sequence[float]):
                Lower and upper values for each channel.
            'LUTs' (list[numpy.ndarray[(3, 256), 'uint8']]):
                Color palettes for each channel.
            'Plot' (bytes):
                Undocumented ImageJ internal format.
            'ROI', 'Overlays' (bytes):
                Undocumented ImageJ internal region of interest and overlay
                format. Can be parsed with the
                `roifile <https://pypi.org/project/roifile/>`_  package.
            'Properties' (dict[str, str]):
                Map of key, value items.

    """

    def _string(data: bytes, byteorder: ByteOrder, /) -> str:
        return data.decode('utf-16' + {'>': 'be', '<': 'le'}[byteorder])

    def _doubles(data: bytes, byteorder: ByteOrder, /) -> tuple[float, ...]:
        return struct.unpack(byteorder + ('d' * (len(data) // 8)), data)

    def _lut(data: bytes, byteorder: ByteOrder, /) -> NDArray[numpy.uint8]:
        return numpy.frombuffer(data, numpy.uint8).reshape(-1, 256)

    def _bytes(data: bytes, byteorder: ByteOrder, /) -> bytes:
        return data

    # big-endian
    metadata_types: dict[
        bytes, tuple[str, Callable[[bytes, ByteOrder], Any]]
    ] = {
        b'info': ('Info', _string),
        b'labl': ('Labels', _string),
        b'rang': ('Ranges', _doubles),
        b'luts': ('LUTs', _lut),
        b'plot': ('Plot', _bytes),
        b'roi ': ('ROI', _bytes),
        b'over': ('Overlays', _bytes),
        b'prop': ('Properties', _string),
    }
    # little-endian
    metadata_types.update({k[::-1]: v for k, v in metadata_types.items()})

    if len(bytecounts) == 0:
        raise ValueError('no ImageJ metadata')

    if data[:4] not in {b'IJIJ', b'JIJI'}:
        raise ValueError('invalid ImageJ metadata')

    header_size = bytecounts[0]
    if header_size < 12 or header_size > 804:
        raise ValueError('invalid ImageJ metadata header size')

    ntypes = (header_size - 4) // 8
    header = struct.unpack(
        byteorder + '4sI' * ntypes, data[4 : 4 + ntypes * 8]
    )
    pos = 4 + ntypes * 8
    counter = 0
    result = {}
    for mtype, count in zip(header[::2], header[1::2]):
        values = []
        name, func = metadata_types.get(mtype, (bytes2str(mtype), _bytes))
        for _ in range(count):
            counter += 1
            pos1 = pos + bytecounts[counter]
            values.append(func(data[pos:pos1], byteorder))
            pos = pos1
        result[name.strip()] = values[0] if count == 1 else values
    prop = result.get('Properties')
    if prop and len(prop) % 2 == 0:
        result['Properties'] = dict(
            prop[i : i + 2] for i in range(0, len(prop), 2)
        )
    return result


def imagej_description_metadata(description: str, /) -> dict[str, Any]:
    r"""Return metatata from ImageJ image description.

    Raise ValueError if not a valid ImageJ description.

    >>> description = 'ImageJ=1.11a\nimages=510\nhyperstack=true\n'
    >>> imagej_description_metadata(description)  # doctest: +SKIP
    {'ImageJ': '1.11a', 'images': 510, 'hyperstack': True}

    """

    def _bool(val: str, /) -> bool:
        return {'true': True, 'false': False}[val.lower()]

    result: dict[str, Any] = {}
    for line in description.splitlines():
        try:
            key, val = line.split('=')
        except Exception:
            continue
        key = key.strip()
        val = val.strip()
        for dtype in (int, float, _bool):
            try:
                val = dtype(val)  # type: ignore[assignment]
                break
            except Exception:
                pass
        result[key] = val

    if 'ImageJ' not in result and 'SCIFIO' not in result:
        raise ValueError(f'not an ImageJ image description {result!r}')
    return result


def imagej_description(
    shape: Sequence[int],
    /,
    axes: str | None = None,
    rgb: bool | None = None,
    colormaped: bool = False,
    **metadata: Any,  # TODO: use TypedDict
) -> str:
    """Return ImageJ image description from data shape and metadata.

    Parameters:
        shape:
            Shape of image array.
        axes:
            Character codes for dimensions in `shape`.
            ImageJ can handle up to 6 dimensions in order TZCYXS.
            `Axes` and `shape` are used to determine the images, channels,
            slices, and frames entries of the image description.
        rgb:
            Image is RGB type.
        colormaped:
            Image is indexed color.
        **metadata:
            Additional items to be included in image description:

            hyperstack (bool):
                Image is a hyperstack.
                The default is True unless `colormapped` is true.
            mode (str):
                Display mode: 'grayscale', 'composite', or 'color'.
                The default is 'grayscale' unless `rgb` or `colormaped` are
                true. Ignored if `hyperstack` is false.
            loop (bool):
                Loop frames back and forth. The default is False.
            finterval (float):
                Frame interval in seconds.
            fps (float):
                Frames per seconds. The inverse of `finterval`.
            spacing (float):
                Voxel spacing in `unit` units.
            unit (str):
                Unit for `spacing` and X/YResolution tags.
                Usually 'um' (micrometer) or 'pixel'.
            xorigin, yorigin, zorigin (float):
                X, Y, and Z origins in pixel units.
            version (str):
                ImageJ version string. The default is '1.11a'.
            images, channels, slices, frames (int):
                Ignored.

    Examples:
        >>> imagej_description((51, 5, 2, 196, 171))  # doctest: +SKIP
        ImageJ=1.11a
        images=510
        channels=2
        slices=5
        frames=51
        hyperstack=true
        mode=grayscale
        loop=false

    """
    mode = metadata.pop('mode', None)
    hyperstack = metadata.pop('hyperstack', None)
    loop = metadata.pop('loop', None)
    version = metadata.pop('ImageJ', '1.11a')

    if colormaped:
        hyperstack = False
        rgb = False

    shape = imagej_shape(shape, rgb=rgb, axes=axes)
    rgb = shape[-1] in {3, 4}

    append = []
    result = [f'ImageJ={version}']
    result.append(f'images={product(shape[:-3])}')
    if hyperstack is None:
        hyperstack = True
        append.append('hyperstack=true')
    else:
        append.append(f'hyperstack={bool(hyperstack)}'.lower())
    if shape[2] > 1:
        result.append(f'channels={shape[2]}')
    if mode is None and not rgb and not colormaped:
        mode = 'grayscale'
    if hyperstack and mode:
        append.append(f'mode={mode}')
    if shape[1] > 1:
        result.append(f'slices={shape[1]}')
    if shape[0] > 1:
        result.append(f'frames={shape[0]}')
        if loop is None:
            append.append('loop=false')
    if loop is not None:
        append.append(f'loop={bool(loop)}'.lower())

    for key, value in metadata.items():
        if key not in {'images', 'channels', 'slices', 'frames', 'SCIFIO'}:
            if isinstance(value, bool):
                value = str(value).lower()
            append.append(f'{key.lower()}={value}')

    return '\n'.join(result + append + [''])


def imagej_shape(
    shape: Sequence[int],
    /,
    *,
    rgb: bool | None = None,
    axes: str | None = None,
) -> tuple[int, ...]:
    """Return shape normalized to 6D ImageJ hyperstack TZCYXS.

    Raise ValueError if not a valid ImageJ hyperstack shape or axes order.

    >>> imagej_shape((2, 3, 4, 5, 3), rgb=False)
    (2, 3, 4, 5, 3, 1)

    """
    shape = tuple(int(i) for i in shape)
    ndim = len(shape)
    if 1 > ndim > 6:
        raise ValueError('ImageJ hyperstack must be 2-6 dimensional')

    if axes:
        if len(axes) != ndim:
            raise ValueError('ImageJ hyperstack shape and axes do not match')
        i = 0
        axes = axes.upper()
        for ax in axes:
            j = 'TZCYXS'.find(ax)
            if j < i:
                raise ValueError(
                    'ImageJ hyperstack axes must be in TZCYXS order'
                )
            i = j
        ndims = len(axes)
        newshape = []
        i = 0
        for ax in 'TZCYXS':
            if i < ndims and ax == axes[i]:
                newshape.append(shape[i])
                i += 1
            else:
                newshape.append(1)
        if newshape[-1] not in {1, 3, 4}:
            raise ValueError(
                'ImageJ hyperstack must contain 1, 3, or 4 samples'
            )
        return tuple(newshape)

    if rgb is None:
        rgb = shape[-1] in {3, 4} and ndim > 2
    if rgb and shape[-1] not in {3, 4}:
        raise ValueError('ImageJ hyperstack is not a RGB image')
    if not rgb and ndim == 6 and shape[-1] != 1:
        raise ValueError('ImageJ hyperstack is not a grayscale image')
    if rgb or shape[-1] == 1:
        return (1,) * (6 - ndim) + shape
    return (1,) * (5 - ndim) + shape + (1,)


def jpeg_decode_colorspace(
    photometric: int,
    planarconfig: int,
    extrasamples: tuple[int, ...],
    jfif: bool,
    /,
) -> tuple[int | None, int | str | None]:
    """Return JPEG and output color space for `jpeg_decode` function."""
    colorspace: int | None = None
    outcolorspace: int | str | None = None
    if extrasamples:
        pass
    elif photometric == 6:
        # YCBCR -> RGB
        outcolorspace = 2  # RGB
    elif photometric == 2:
        # RGB -> RGB
        if not jfif:
            # found in Aperio SVS
            colorspace = 2
        outcolorspace = 2
    elif photometric == 5:
        # CMYK
        outcolorspace = 4
    elif photometric > 3:
        outcolorspace = PHOTOMETRIC(photometric).name
    if planarconfig != 1:
        outcolorspace = 1  # decode separate planes to grayscale
    return colorspace, outcolorspace


def jpeg_shape(jpeg: bytes, /) -> tuple[int, int, int, int]:
    """Return bitdepth and shape of JPEG image."""
    i = 0
    while i < len(jpeg):
        marker = struct.unpack('>H', jpeg[i : i + 2])[0]
        i += 2

        if marker == 0xFFD8:
            # start of image
            continue
        if marker == 0xFFD9:
            # end of image
            break
        if 0xFFD0 <= marker <= 0xFFD7:
            # restart marker
            continue
        if marker == 0xFF01:
            # private marker
            continue

        length = struct.unpack('>H', jpeg[i : i + 2])[0]
        i += 2

        if 0xFFC0 <= marker <= 0xFFC3:
            # start of frame
            return struct.unpack('>BHHB', jpeg[i : i + 6])
        if marker == 0xFFDA:
            # start of scan
            break

        # skip to next marker
        i += length - 2

    raise ValueError('no SOF marker found')


def ndpi_jpeg_tile(jpeg: bytes, /) -> tuple[int, int, bytes]:
    """Return tile shape and JPEG header from JPEG with restart markers."""
    marker: int
    length: int
    factor: int
    ncomponents: int
    restartinterval: int = 0
    sofoffset: int = 0
    sosoffset: int = 0
    i: int = 0
    while i < len(jpeg):
        marker = struct.unpack('>H', jpeg[i : i + 2])[0]
        i += 2

        if marker == 0xFFD8:
            # start of image
            continue
        if marker == 0xFFD9:
            # end of image
            break
        if 0xFFD0 <= marker <= 0xFFD7:
            # restart marker
            continue
        if marker == 0xFF01:
            # private marker
            continue

        length = struct.unpack('>H', jpeg[i : i + 2])[0]
        i += 2

        if marker == 0xFFDD:
            # define restart interval
            restartinterval = struct.unpack('>H', jpeg[i : i + 2])[0]

        elif marker == 0xFFC0:
            # start of frame
            sofoffset = i + 1
            precision, imlength, imwidth, ncomponents = struct.unpack(
                '>BHHB', jpeg[i : i + 6]
            )
            i += 6
            mcuwidth = 1
            mcuheight = 1
            for _ in range(ncomponents):
                cid, factor, table = struct.unpack('>BBB', jpeg[i : i + 3])
                i += 3
                if factor >> 4 > mcuwidth:
                    mcuwidth = factor >> 4
                if factor & 0b00001111 > mcuheight:
                    mcuheight = factor & 0b00001111
            mcuwidth *= 8
            mcuheight *= 8
            i = sofoffset - 1

        elif marker == 0xFFDA:
            # start of scan
            sosoffset = i + length - 2
            break

        # skip to next marker
        i += length - 2

    if restartinterval == 0 or sofoffset == 0 or sosoffset == 0:
        raise ValueError('missing required JPEG markers')

    # patch jpeg header for tile size
    tilelength = mcuheight
    tilewidth = restartinterval * mcuwidth
    jpegheader = (
        jpeg[:sofoffset]
        + struct.pack('>HH', tilelength, tilewidth)
        + jpeg[sofoffset + 4 : sosoffset]
    )
    return tilelength, tilewidth, jpegheader


def shaped_description(shape: Sequence[int], /, **metadata: Any) -> str:
    """Return JSON image description from data shape and other metadata.

    Return UTF-8 encoded JSON.

    >>> shaped_description((256, 256, 3), axes='YXS')  # doctest: +SKIP
    '{"shape": [256, 256, 3], "axes": "YXS"}'

    """
    metadata.update(shape=shape)
    return json.dumps(metadata)  # .encode()


def shaped_description_metadata(description: str, /) -> dict[str, Any]:
    """Return metatata from JSON formatted image description.

    Raise ValueError if `description` is of unknown format.

    >>> description = '{"shape": [256, 256, 3], "axes": "YXS"}'
    >>> shaped_description_metadata(description)  # doctest: +SKIP
    {'shape': [256, 256, 3], 'axes': 'YXS'}
    >>> shaped_description_metadata('shape=(256, 256, 3)')
    {'shape': (256, 256, 3)}

    """
    if description[:6] == 'shape=':
        # old-style 'shaped' description; not JSON
        shape = tuple(int(i) for i in description[7:-1].split(','))
        return {'shape': shape}
    if description[:1] == '{' and description[-1:] == '}':
        # JSON description
        return json.loads(description)
    raise ValueError('invalid JSON image description', description)


def fluoview_description_metadata(
    description: str,
    /,
    ignoresections: Container[str] | None = None,
) -> dict[str, Any]:
    r"""Return metatata from FluoView image description.

    The FluoView image description format is unspecified. Expect failures.

    >>> descr = (
    ...     '[Intensity Mapping]\nMap Ch0: Range=00000 to 02047\n'
    ...     '[Intensity Mapping End]'
    ... )
    >>> fluoview_description_metadata(descr)
    {'Intensity Mapping': {'Map Ch0: Range': '00000 to 02047'}}

    """
    if not description.startswith('['):
        raise ValueError('invalid FluoView image description')
    if ignoresections is None:
        ignoresections = {'Region Info (Fields)', 'Protocol Description'}

    section: Any
    result: dict[str, Any] = {}
    sections = [result]
    comment = False
    for line in description.splitlines():
        if not comment:
            line = line.strip()
        if not line:
            continue
        if line[0] == '[':
            if line[-5:] == ' End]':
                # close section
                del sections[-1]
                section = sections[-1]
                name = line[1:-5]
                if comment:
                    section[name] = '\n'.join(section[name])
                if name[:4] == 'LUT ':
                    a = numpy.array(section[name], dtype=numpy.uint8)
                    a.shape = -1, 3
                    section[name] = a
                continue
            # new section
            comment = False
            name = line[1:-1]
            if name[:4] == 'LUT ':
                section = []
            elif name in ignoresections:
                section = []
                comment = True
            else:
                section = {}
            sections.append(section)
            result[name] = section
            continue
        # add entry
        if comment:
            section.append(line)
            continue
        lines = line.split('=', 1)
        if len(line) == 1:
            section[lines[0].strip()] = None
            continue
        key, value = lines
        if key[:4] == 'RGB ':
            section.extend(int(rgb) for rgb in value.split())
        else:
            section[key.strip()] = astype(value.strip())
    return result


def pilatus_description_metadata(description: str, /) -> dict[str, Any]:
    """Return metatata from Pilatus image description.

    Return metadata from Pilatus pixel array detectors by Dectris, created
    by camserver or TVX software.

    >>> pilatus_description_metadata('# Pixel_size 172e-6 m x 172e-6 m')
    {'Pixel_size': (0.000172, 0.000172)}

    """
    result: dict[str, Any] = {}
    values: Any
    if not description.startswith('# '):
        return result
    for c in '#:=,()':
        description = description.replace(c, ' ')
    for lines in description.split('\n'):
        if lines[:2] != '  ':
            continue
        line = lines.split()
        name = line[0]
        if line[0] not in TIFF.PILATUS_HEADER:
            try:
                result['DateTime'] = strptime(
                    ' '.join(line), '%Y-%m-%dT%H %M %S.%f'
                )
            except Exception:
                result[name] = ' '.join(line[1:])
            continue
        indices, dtype = TIFF.PILATUS_HEADER[line[0]]
        if isinstance(indices[0], slice):
            # assumes one slice
            values = line[indices[0]]
        else:
            values = [line[i] for i in indices]
        if dtype is float and values[0] == 'not':
            values = ['NaN']
        values = tuple(dtype(v) for v in values)
        if dtype is str:
            values = ' '.join(values)
        elif len(values) == 1:
            values = values[0]
        result[name] = values
    return result


def svs_description_metadata(description: str, /) -> dict[str, Any]:
    """Return metatata from Aperio image description.

    The Aperio image description format is unspecified. Expect failures.

    >>> svs_description_metadata('Aperio Image Library v1.0')
    {'Header': 'Aperio Image Library v1.0'}

    """
    if not description.startswith('Aperio '):
        raise ValueError('invalid Aperio image description')
    result = {}
    items = description.split('|')
    result['Header'] = items[0]
    for item in items[1:]:
        key, value = item.split('=', maxsplit=1)
        result[key.strip()] = astype(value.strip())
    return result


def stk_description_metadata(description: str, /) -> list[dict[str, Any]]:
    """Return metadata from MetaMorph image description.

    The MetaMorph image description format is unspecified. Expect failures.

    """
    description = description.strip()
    if not description:
        return []
    # try:
    #     description = bytes2str(description)
    # except UnicodeDecodeError as exc:
    #     logger().warning(
    #         '<tifffile.stk_description_metadata> raised {exc!r:.128}'
    #     )
    #     return []
    result = []
    for plane in description.split('\x00'):
        d = {}
        for line in plane.split('\r\n'):
            lines = line.split(':', 1)
            if len(lines) > 1:
                name, value = lines
                d[name.strip()] = astype(value.strip())
            else:
                value = lines[0].strip()
                if value:
                    if '' in d:
                        d[''].append(value)
                    else:
                        d[''] = [value]
        result.append(d)
    return result


def metaseries_description_metadata(description: str, /) -> dict[str, Any]:
    """Return metatata from MetaSeries image description."""
    if not description.startswith('<MetaData>'):
        raise ValueError('invalid MetaSeries image description')

    import uuid
    from xml.etree import ElementTree as etree

    root = etree.fromstring(description)
    types: dict[str, Callable[..., Any]] = {
        'float': float,
        'int': int,
        'bool': lambda x: asbool(x, 'on', 'off'),
        'time': lambda x: strptime(x, '%Y%m%d %H:%M:%S.%f'),
        'guid': uuid.UUID,
        # 'float-array':
        # 'colorref':
    }

    def parse(
        root: etree.Element, result: dict[str, Any], /
    ) -> dict[str, Any]:
        # recursive
        for child in root:
            attrib = child.attrib
            if not attrib:
                result[child.tag] = parse(child, {})
                continue
            if 'id' in attrib:
                i = attrib['id']
                t = attrib['type']
                v = attrib['value']
                if t in types:
                    try:
                        result[i] = types[t](v)
                    except Exception:
                        result[i] = v
                else:
                    result[i] = v
        return result

    adict = parse(root, {})
    if 'Description' in adict:
        adict['Description'] = adict['Description'].replace('&#13;&#10;', '\n')
    return adict


def scanimage_description_metadata(description: str, /) -> Any:
    """Return metatata from ScanImage image description."""
    return matlabstr2py(description)


def scanimage_artist_metadata(artist: str, /) -> dict[str, Any] | None:
    """Return metatata from ScanImage artist tag."""
    try:
        return json.loads(artist)
    except ValueError as exc:
        logger().warning(
            f'<tifffile.scanimage_artist_metadata> raised {exc!r:.128}'
        )
    return None


def olympusini_metadata(inistr: str, /) -> dict[str, Any]:
    """Return OlympusSIS metadata from INI string.

    No specification is available.

    """

    def keyindex(key: str, /) -> tuple[str, int]:
        # split key into name and index
        index = 0
        i = len(key.rstrip('0123456789'))
        if i < len(key):
            index = int(key[i:]) - 1
            key = key[:i]
        return key, index

    result: dict[str, Any] = {}
    bands: list[dict[str, Any]] = []
    value: Any
    zpos: list[Any] | None = None
    tpos: list[Any] | None = None
    for line in inistr.splitlines():
        line = line.strip()
        if line == '' or line[0] == ';':
            continue
        if line[0] == '[' and line[-1] == ']':
            section_name = line[1:-1]
            result[section_name] = section = {}
            if section_name == 'Dimension':
                result['axes'] = axes = []
                result['shape'] = shape = []
            elif section_name == 'ASD':
                result[section_name] = []
            elif section_name == 'Z':
                if 'Dimension' in result:
                    result[section_name]['ZPos'] = zpos = []
            elif section_name == 'Time':
                if 'Dimension' in result:
                    result[section_name]['TimePos'] = tpos = []
            elif section_name == 'Band':
                nbands = result['Dimension']['Band']
                bands = [{'LUT': []} for _ in range(nbands)]
                result[section_name] = bands
                iband = 0
        else:
            key, value = line.split('=')
            if value.strip() == '':
                value = None
            elif ',' in value:
                value = tuple(astype(v) for v in value.split(','))
            else:
                value = astype(value)

            if section_name == 'Dimension':
                section[key] = value
                axes.append(key)
                shape.append(value)
            elif section_name == 'ASD':
                if key == 'Count':
                    result['ASD'] = [{}] * value
                else:
                    key, index = keyindex(key)
                    result['ASD'][index][key] = value
            elif section_name == 'Band':
                if key[:3] == 'LUT':
                    lut = bands[iband]['LUT']
                    value = struct.pack('<I', value)
                    lut.append(
                        [ord(value[0:1]), ord(value[1:2]), ord(value[2:3])]
                    )
                else:
                    key, iband = keyindex(key)
                    bands[iband][key] = value
            elif key[:4] == 'ZPos' and zpos is not None:
                zpos.append(value)
            elif key[:7] == 'TimePos' and tpos is not None:
                tpos.append(value)
            else:
                section[key] = value

    if 'axes' in result:
        sisaxes = {'Band': 'C'}
        axes = []
        shape = []
        for i, x in zip(result['shape'], result['axes']):
            if i > 1:
                axes.append(sisaxes.get(x, x[0].upper()))
                shape.append(i)
        result['axes'] = ''.join(axes)
        result['shape'] = tuple(shape)
    try:
        result['Z']['ZPos'] = numpy.array(
            result['Z']['ZPos'][: result['Dimension']['Z']], numpy.float64
        )
    except Exception:
        pass
    try:
        result['Time']['TimePos'] = numpy.array(
            result['Time']['TimePos'][: result['Dimension']['Time']],
            numpy.int32,
        )
    except Exception:
        pass
    for band in bands:
        band['LUT'] = numpy.array(band['LUT'], numpy.uint8)
    return result


def astrotiff_description_metadata(
    description: str, /, sep: str = ':'
) -> dict[str, Any]:
    """Return metatata from AstroTIFF image description."""
    logmsg = '<tifffile.astrotiff_description_metadata> '
    counts: dict[str, int] = {}
    result: dict[str, Any] = {}
    value: Any
    for line in description.splitlines():
        line = line.strip()
        if not line:
            continue

        key = line[:8].strip()
        value = line[8:]

        if not value.startswith('='):
            # for example, COMMENT or HISTORY
            if key + f'{sep}0' not in result:
                result[key + f'{sep}0'] = value
                counts[key] = 1
            else:
                result[key + f'{sep}{counts[key]}'] = value
                counts[key] += 1
            continue

        value = value[1:]
        if '/' in value:
            value, comment = value.split('/', 1)
            comment = comment.strip()
        else:
            comment = ''
        value = value.strip()

        if not value:
            # undefined
            value = None
        elif value[0] == "'":
            # string
            if len(value) < 2:
                logger().warning(logmsg + f'{key}: invalid string {value!r}')
                continue
            if value[-1] == "'":
                value = value[1:-1]
            else:
                # string containing '/'
                if not ("'" in comment and '/' in comment):
                    logger().warning(
                        logmsg + f'{key}: invalid string {value!r}'
                    )
                    continue
                value, comment = line[9:].strip()[1:].split("'", 1)
                comment = comment.split('/', 1)[-1].strip()
            # TODO: string containing single quote '
        elif value[0] == '(' and value[-1] == ')':
            # complex number
            value = value[1:-1]
            dtype = float if '.' in value else int
            value = tuple(dtype(v.strip()) for v in value.split(','))
        elif value == 'T':
            value = True
        elif value == 'F':
            value = False
        elif '.' in value:
            value = float(value)
        else:
            try:
                value = int(value)
            except Exception:
                logger().warning(logmsg + f'{key}: invalid value {value!r}')
                continue

        if key in result:
            logger().warning(logmsg + f'{key}: duplicate key')

        result[key] = value
        if comment:
            result[key + f'{sep}COMMENT'] = comment
            if comment[0] == '[' and ']' in comment:
                result[key + f'{sep}UNIT'] = comment[1:].split(']', 1)[0]

    return result


def streak_description_metadata(
    description: str, fh: FileHandle, /
) -> dict[str, Any]:
    """Return metatata from Hamamatsu streak image description."""
    section_pattern = re.compile(
        r'\[([a-zA-Z0-9 _\-\.]+)\],([^\[]*)', re.DOTALL
    )
    properties_pattern = re.compile(
        r'([a-zA-Z0-9 _\-\.]+)=(\"[^\"]*\"|[\+\-0-9\.]+|[^,]*)'
    )
    result: dict[str, Any] = {}
    for section, values in section_pattern.findall(description.strip()):
        properties = {}
        for key, value in properties_pattern.findall(values):
            value = value.strip()
            if not value or value == '"':
                value = None
            elif value[0] == '"' and value[-1] == '"':
                value = value[1:-1]
            if ',' in value:
                try:
                    value = tuple(
                        (
                            float(v)
                            if '.' in value
                            else int(v[1:] if v[0] == '#' else v)
                        )
                        for v in value.split(',')
                    )
                except ValueError:
                    pass
            elif '.' in value:
                try:
                    value = float(value)
                except ValueError:
                    pass
            else:
                try:
                    value = int(value)
                except ValueError:
                    pass
            properties[key] = value
        result[section] = properties

    if not fh.closed:
        pos = fh.tell()
        for scaling in ('ScalingXScaling', 'ScalingYScaling'):
            try:
                offset, count = result['Scaling'][scaling + 'File']
                fh.seek(offset)
                result['Scaling'][scaling] = fh.read_array(
                    dtype='<f4', count=count
                )
            except Exception:
                pass
        fh.seek(pos)

    return result


def unpack_rgb(
    data: bytes,
    /,
    dtype: DTypeLike | None = None,
    bitspersample: tuple[int, ...] | None = None,
    rescale: bool = True,
) -> NDArray[Any]:
    """Return array from bytes containing packed samples.

    Use to unpack RGB565 or RGB555 to RGB888 format.
    Works on little-endian platforms only.

    Parameters:
        data:
            Bytes to be decoded.
            Samples in each pixel are stored consecutively.
            Pixels are aligned to 8, 16, or 32 bit boundaries.
        dtype:
            Data type of samples.
            The byte order applies also to the data stream.
        bitspersample:
            Number of bits for each sample in pixel.
        rescale:
            Upscale samples to number of bits in dtype.

    Returns:
        Flattened array of unpacked samples of native dtype.

    Examples:
        >>> data = struct.pack('BBBB', 0x21, 0x08, 0xFF, 0xFF)
        >>> print(unpack_rgb(data, '<B', (5, 6, 5), False))
        [ 1  1  1 31 63 31]
        >>> print(unpack_rgb(data, '<B', (5, 6, 5)))
        [  8   4   8 255 255 255]
        >>> print(unpack_rgb(data, '<B', (5, 5, 5)))
        [ 16   8   8 255 255 255]

    """
    if bitspersample is None:
        bitspersample = (5, 6, 5)
    if dtype is None:
        dtype = '<B'
    dtype = numpy.dtype(dtype)
    bits = int(numpy.sum(bitspersample))
    if not (
        bits <= 32 and all(i <= dtype.itemsize * 8 for i in bitspersample)
    ):
        raise ValueError(f'sample size not supported: {bitspersample}')
    dt = next(i for i in 'BHI' if numpy.dtype(i).itemsize * 8 >= bits)
    data_array = numpy.frombuffer(data, dtype.byteorder + dt)
    result = numpy.empty((data_array.size, len(bitspersample)), dtype.char)
    for i, bps in enumerate(bitspersample):
        t = data_array >> int(numpy.sum(bitspersample[i + 1 :]))
        t &= int('0b' + '1' * bps, 2)
        if rescale:
            o = ((dtype.itemsize * 8) // bps + 1) * bps
            if o > data_array.dtype.itemsize * 8:
                t = t.astype('I')
            t *= (2**o - 1) // (2**bps - 1)
            t //= 2 ** (o - (dtype.itemsize * 8))
        result[:, i] = t
    return result.reshape(-1)


def apply_colormap(
    image: NDArray[Any], colormap: NDArray[Any], /, contig: bool = True
) -> NDArray[Any]:
    """Return palette-colored image.

    The image array values are used to index the colormap on axis 1.
    The returned image array is of shape `image.shape+colormap.shape[0]`
    and dtype `colormap.dtype`.

    Parameters:
        image:
            Array of indices into colormap.
        colormap:
            RGB lookup table aka palette of shape `(3, 2**bitspersample)`.
        contig:
            Return contiguous array.

    Examples:
        >>> import numpy
        >>> im = numpy.arange(256, dtype='uint8')
        >>> colormap = numpy.vstack([im, im, im]).astype('uint16') * 256
        >>> apply_colormap(im, colormap)[-1]
        array([65280, 65280, 65280], dtype=uint16)

    """
    image = numpy.take(colormap, image, axis=1)
    image = numpy.rollaxis(image, 0, image.ndim)
    if contig:
        image = numpy.ascontiguousarray(image)
    return image


def parse_filenames(
    files: Sequence[str],
    /,
    pattern: str | None = None,
    axesorder: Sequence[int] | None = None,
    categories: dict[str, dict[str, int]] | None = None,
    *,
    _shape: Sequence[int] | None = None,
) -> tuple[
    tuple[str, ...], tuple[int, ...], list[tuple[int, ...]], Sequence[str]
]:
    r"""Return shape and axes from sequence of file names matching pattern.

    Parameters:
        files:
            Sequence of file names to parse.
        pattern:
            Regular expression pattern matching axes names and chunk indices
            in file names.
            By default, no pattern matching is performed.
            Axes names can be specified by matching groups preceding the index
            groups in the file name, be provided as group names for the index
            groups, or be omitted.
            The predefined 'axes' pattern matches Olympus OIF and Leica TIFF
            series.
        axesorder:
            Indices of axes in pattern. By default, axes are returned in the
            order they appear in pattern.
        categories:
            Map of index group matches to integer indices.
            `{'axislabel': {'category': index}}`
        _shape:
            Shape of file sequence. The default is
            `maximum - minimum + 1` of the parsed indices for each dimension.

    Returns:
        - Axes names for each dimension.
        - Shape of file series.
        - Index of each file in shape.
        - Filtered sequence of file names.

    Examples:
        >>> parse_filenames(
        ...     ['c1001.ext', 'c2002.ext'], r'([^\d])(\d)(?P<t>\d+)\.ext'
        ... )
        (('c', 't'), (2, 2), [(0, 0), (1, 1)], ['c1001.ext', 'c2002.ext'])

    """
    # TODO: add option to filter files that do not match pattern

    shape = None if _shape is None else tuple(_shape)
    if pattern is None:
        if shape is not None and (len(shape) != 1 or shape[0] < len(files)):
            raise ValueError(
                f'shape {(len(files),)} does not fit provided shape {shape}'
            )
        return (
            ('I',),
            (len(files),),
            [(i,) for i in range(len(files))],
            files,
        )

    pattern = TIFF.FILE_PATTERNS.get(pattern, pattern)
    if not pattern:
        raise ValueError('invalid pattern')
    pattern_compiled: Any
    if isinstance(pattern, str):
        pattern_compiled = re.compile(pattern)
    elif hasattr(pattern, 'groupindex'):
        pattern_compiled = pattern
    else:
        raise ValueError('invalid pattern')

    if categories is None:
        categories = {}

    def parse(fname: str, /) -> tuple[tuple[str, ...], tuple[int, ...]]:
        # return axes names and indices from file name
        assert categories is not None
        dims: list[str] = []
        indices: list[int] = []
        groupindex = {v: k for k, v in pattern_compiled.groupindex.items()}
        match = pattern_compiled.search(fname)
        if match is None:
            raise ValueError(f'pattern does not match file name {fname!r}')
        ax = None
        for i, m in enumerate(match.groups()):
            if m is None:
                continue
            if i + 1 in groupindex:
                ax = groupindex[i + 1]
            elif m[0].isalpha():
                ax = m  # axis label for next index
                continue
            if ax is None:
                ax = 'Q'  # no preceding axis letter
            try:
                if ax in categories:
                    m = categories[ax][m]
                m = int(m)
            except Exception as exc:
                raise ValueError(f'invalid index {m!r}') from exc
            indices.append(m)
            dims.append(ax)
            ax = None
        return tuple(dims), tuple(indices)

    normpaths = [os.path.normpath(f) for f in files]
    if len(normpaths) == 1:
        prefix_str = os.path.dirname(normpaths[0])
    else:
        prefix_str = os.path.commonpath(normpaths)
    prefix = len(prefix_str)

    dims: tuple[str, ...] | None = None
    indices: list[tuple[int, ...]] = []
    for fname in normpaths:
        lbl, idx = parse(fname[prefix:])
        if dims is None:
            dims = lbl
            if axesorder is not None and (
                len(axesorder) != len(dims)
                or any(i not in axesorder for i in range(len(dims)))
            ):
                raise ValueError(
                    f'invalid axesorder {axesorder!r} for {dims!r}'
                )
        elif dims != lbl:
            raise ValueError('dims do not match within image sequence')
        if axesorder is not None:
            idx = tuple(idx[i] for i in axesorder)
        indices.append(idx)

    assert dims is not None
    if axesorder is not None:
        dims = tuple(dims[i] for i in axesorder)

    # determine shape
    indices_array = numpy.array(indices, dtype=numpy.intp)
    parsedshape = numpy.max(indices, axis=0)

    if shape is None:
        startindex = numpy.min(indices_array, axis=0)
        indices_array -= startindex
        parsedshape -= startindex
        parsedshape += 1
        shape = tuple(int(i) for i in parsedshape.tolist())
    elif len(parsedshape) != len(shape) or any(
        i > j for i, j in zip(shape, parsedshape)
    ):
        raise ValueError(
            f'parsed shape {parsedshape} does not fit provided shape {shape}'
        )

    indices_list: list[list[int]]
    indices_list = indices_array.tolist()
    indices = [tuple(index) for index in indices_list]

    return dims, shape, indices, files


def iter_images(data: NDArray[Any], /) -> Iterator[NDArray[Any]]:
    """Return iterator over pages in data array of normalized shape."""
    yield from data


def iter_strips(
    pageiter: Iterator[NDArray[Any] | None],
    shape: tuple[int, ...],
    dtype: numpy.dtype[Any],
    rowsperstrip: int,
    /,
) -> Iterator[NDArray[Any]]:
    """Return iterator over strips in pages."""
    numstrips = (shape[-3] + rowsperstrip - 1) // rowsperstrip

    for iteritem in pageiter:
        if iteritem is None:
            # for _ in range(numstrips):
            #     yield None
            # continue
            pagedata = numpy.zeros(shape, dtype)
        else:
            pagedata = iteritem.reshape(shape)
        for plane in pagedata:
            for depth in plane:
                for i in range(numstrips):
                    yield depth[i * rowsperstrip : (i + 1) * rowsperstrip]


def iter_tiles(
    data: NDArray[Any],
    tile: tuple[int, ...],
    tiles: tuple[int, ...],
    /,
) -> Iterator[NDArray[Any]]:
    """Return iterator over full tiles in data array of normalized shape.

    Tiles are zero-padded if necessary.

    """
    if not 1 < len(tile) < 4 or len(tile) != len(tiles):
        raise ValueError('invalid tile or tiles shape')
    chunkshape = tile + (data.shape[-1],)
    chunksize = product(chunkshape)
    dtype = data.dtype
    sz, sy, sx = data.shape[2:5]
    if len(tile) == 2:
        y, x = tile
        for page in data:
            for plane in page:
                for ty in range(tiles[0]):
                    ty *= y
                    cy = min(y, sy - ty)
                    for tx in range(tiles[1]):
                        tx *= x
                        cx = min(x, sx - tx)
                        chunk = plane[0, ty : ty + cy, tx : tx + cx]
                        if chunk.size != chunksize:
                            chunk_ = numpy.zeros(chunkshape, dtype)
                            chunk_[:cy, :cx] = chunk
                            chunk = chunk_
                        yield chunk
    else:
        z, y, x = tile
        for page in data:
            for plane in page:
                for tz in range(tiles[0]):
                    tz *= z
                    cz = min(z, sz - tz)
                    for ty in range(tiles[1]):
                        ty *= y
                        cy = min(y, sy - ty)
                        for tx in range(tiles[2]):
                            tx *= x
                            cx = min(x, sx - tx)
                            chunk = plane[
                                tz : tz + cz, ty : ty + cy, tx : tx + cx
                            ]
                            if chunk.size != chunksize:
                                chunk_ = numpy.zeros(chunkshape, dtype)
                                chunk_[:cz, :cy, :cx] = chunk
                                chunk = chunk_
                            yield chunk[0] if z == 1 else chunk


def encode_chunks(
    numchunks: int,
    chunkiter: Iterator[NDArray[Any] | None],
    encode: Callable[[NDArray[Any]], bytes],
    shape: Sequence[int],
    dtype: numpy.dtype[Any],
    maxworkers: int | None,
    buffersize: int | None,
    istile: bool,
    /,
) -> Iterator[bytes]:
    """Return iterator over encoded chunks."""
    if numchunks <= 0:
        return

    chunksize = product(shape) * dtype.itemsize

    if istile:
        # pad tiles
        def func(chunk: NDArray[Any] | None, /) -> bytes:
            if chunk is None:
                return b''
            chunk = numpy.ascontiguousarray(chunk, dtype)
            if chunk.nbytes != chunksize:
                # if chunk.dtype != dtype:
                #     raise ValueError('dtype of chunk does not match data')
                pad = tuple((0, i - j) for i, j in zip(shape, chunk.shape))
                chunk = numpy.pad(chunk, pad)
            return encode(chunk)

    else:
        # strips
        def func(chunk: NDArray[Any] | None, /) -> bytes:
            if chunk is None:
                return b''
            chunk = numpy.ascontiguousarray(chunk, dtype)
            return encode(chunk)

    if maxworkers is None or maxworkers < 2 or numchunks < 2:
        for _ in range(numchunks):
            chunk = next(chunkiter)
            # assert chunk is None or isinstance(chunk, numpy.ndarray)
            yield func(chunk)
            del chunk
        return

    # because ThreadPoolExecutor.map is not collecting items lazily, reduce
    # memory overhead by processing chunks iterator maxchunks items at a time
    if buffersize is None:
        buffersize = TIFF.BUFFERSIZE * 2
    maxchunks = max(maxworkers, buffersize // chunksize)

    if numchunks <= maxchunks:

        def chunks() -> Iterator[NDArray[Any] | None]:
            for _ in range(numchunks):
                chunk = next(chunkiter)
                # assert chunk is None or isinstance(chunk, numpy.ndarray)
                yield chunk
                del chunk

        with ThreadPoolExecutor(maxworkers) as executor:
            yield from executor.map(func, chunks())
        return

    with ThreadPoolExecutor(maxworkers) as executor:
        count = 1
        chunk_list = []
        for _ in range(numchunks):
            chunk = next(chunkiter)
            if chunk is not None:
                count += 1
            # assert chunk is None or isinstance(chunk, numpy.ndarray)
            chunk_list.append(chunk)
            if count == maxchunks:
                yield from executor.map(func, chunk_list)
                chunk_list.clear()
                count = 0
        if chunk_list:
            yield from executor.map(func, chunk_list)


def zarr_selection(
    store: ZarrStore,
    selection: Any,
    /,
    *,
    groupindex: int | None = None,
    close: bool = True,
    out: OutputType = None,
) -> NDArray[Any]:
    """Return selection from Zarr 2 store.

    Parameters:
        store:
            ZarrStore instance to read selection from.
        selection:
            Subset of image to be extracted and returned.
            Refer to the Zarr 2 documentation for valid selections.
        groupindex:
            Index of array if store is Zarr 2 group.
        close:
            Close store before returning.
        out:
            Specifies how image array is returned.
            By default, create a new array.
            If a *numpy.ndarray*, a writable array to which the images
            are copied.
            If *'memmap'*, create a memory-mapped array in a temporary
            file.
            If a *string* or *open file*, the file used to create a
            memory-mapped array.

    """
    import zarr

    try:
        import zarr.hierarchy
        import zarr.indexing
    except ImportError as exc:
        raise ValueError(
            f'zarr {zarr.__version__} >= 3 is not supported'
        ) from exc

    z = zarr.open(store, mode='r')
    try:
        if isinstance(z, zarr.hierarchy.Group):
            if groupindex is None:
                groupindex = 0
            z = z[groupindex]
        if out is not None:
            shape = zarr.indexing.BasicIndexer(selection, z).shape
            out = create_output(out, shape, z.dtype)
        result = z.get_basic_selection(selection, out=out)
    finally:
        if close:
            store.close()
    return result


def reorient(
    image: NDArray[Any], orientation: ORIENTATION | int | str, /
) -> NDArray[Any]:
    """Return reoriented view of image array.

    Parameters:
        image:
            Non-squeezed output of `asarray` functions.
            Axes -3 and -2 must be image length and width respectively.
        orientation:
            Value of Orientation tag.

    """
    orientation = cast(ORIENTATION, enumarg(ORIENTATION, orientation))

    if orientation == ORIENTATION.TOPLEFT:
        return image
    if orientation == ORIENTATION.TOPRIGHT:
        return image[..., ::-1, :]
    if orientation == ORIENTATION.BOTLEFT:
        return image[..., ::-1, :, :]
    if orientation == ORIENTATION.BOTRIGHT:
        return image[..., ::-1, ::-1, :]
    if orientation == ORIENTATION.LEFTTOP:
        return numpy.swapaxes(image, -3, -2)
    if orientation == ORIENTATION.RIGHTTOP:
        return numpy.swapaxes(image, -3, -2)[..., ::-1, :]
    if orientation == ORIENTATION.RIGHTBOT:
        return numpy.swapaxes(image, -3, -2)[..., ::-1, :, :]
    if orientation == ORIENTATION.LEFTBOT:
        return numpy.swapaxes(image, -3, -2)[..., ::-1, ::-1, :]
    return image


def repeat_nd(a: ArrayLike, repeats: Sequence[int], /) -> NDArray[Any]:
    """Return read-only view into input array with elements repeated.

    Zoom image array by integer factors using nearest neighbor interpolation
    (box filter).

    Parameters:
        a: Input array.
        repeats: Number of repetitions to apply along each dimension of input.

    Examples:
        >>> repeat_nd([[1, 2], [3, 4]], (2, 2))
        array([[1, 1, 2, 2],
               [1, 1, 2, 2],
               [3, 3, 4, 4],
               [3, 3, 4, 4]])

    """
    reshape: list[int] = []
    shape: list[int] = []
    strides: list[int] = []
    a = numpy.asarray(a)
    for i, j, k in zip(a.strides, a.shape, repeats):
        shape.extend((j, k))
        strides.extend((i, 0))
        reshape.append(j * k)
    return numpy.lib.stride_tricks.as_strided(
        a, shape, strides, writeable=False
    ).reshape(reshape)


@overload
def reshape_nd(
    data_or_shape: tuple[int, ...], ndim: int, /
) -> tuple[int, ...]: ...


@overload
def reshape_nd(data_or_shape: NDArray[Any], ndim: int, /) -> NDArray[Any]: ...


def reshape_nd(
    data_or_shape: tuple[int, ...] | NDArray[Any], ndim: int, /
) -> tuple[int, ...] | NDArray[Any]:
    """Return image array or shape with at least `ndim` dimensions.

    Prepend 1s to image shape as necessary.

    >>> import numpy
    >>> reshape_nd(numpy.empty(0), 1).shape
    (0,)
    >>> reshape_nd(numpy.empty(1), 2).shape
    (1, 1)
    >>> reshape_nd(numpy.empty((2, 3)), 3).shape
    (1, 2, 3)
    >>> reshape_nd(numpy.empty((3, 4, 5)), 3).shape
    (3, 4, 5)
    >>> reshape_nd((2, 3), 3)
    (1, 2, 3)

    """
    if isinstance(data_or_shape, tuple):
        shape = data_or_shape
    else:
        shape = data_or_shape.shape
    if len(shape) >= ndim:
        return data_or_shape
    shape = (1,) * (ndim - len(shape)) + shape
    if isinstance(data_or_shape, tuple):
        return shape
    return data_or_shape.reshape(shape)


@overload
def squeeze_axes(
    shape: Sequence[int],
    axes: str,
    /,
    skip: str | None = None,
) -> tuple[tuple[int, ...], str, tuple[bool, ...]]: ...


@overload
def squeeze_axes(
    shape: Sequence[int],
    axes: Sequence[str],
    /,
    skip: Sequence[str] | None = None,
) -> tuple[tuple[int, ...], Sequence[str], tuple[bool, ...]]: ...


def squeeze_axes(
    shape: Sequence[int],
    axes: str | Sequence[str],
    /,
    skip: str | Sequence[str] | None = None,
) -> tuple[tuple[int, ...], str | Sequence[str], tuple[bool, ...]]:
    """Return shape and axes with length-1 dimensions removed.

    Remove unused dimensions unless their axes are listed in `skip`.

    Parameters:
        shape:
            Sequence of dimension sizes.
        axes:
            Character codes for dimensions in `shape`.
        skip:
            Character codes for dimensions whose length-1 dimensions are
            not removed. The default is 'XY'.

    Returns:
        shape:
            Sequence of dimension sizes with length-1 dimensions removed.
        axes:
            Character codes for dimensions in output `shape`.
        squeezed:
            Dimensions were kept (True) or removed (False).

    Examples:
        >>> squeeze_axes((5, 1, 2, 1, 1), 'TZYXC')
        ((5, 2, 1), 'TYX', (True, False, True, True, False))
        >>> squeeze_axes((1,), 'Q')
        ((1,), 'Q', (True,))

    """
    if len(shape) != len(axes):
        raise ValueError('dimensions of axes and shape do not match')
    if not axes:
        return tuple(shape), axes, ()
    if skip is None:
        skip = 'X', 'Y', 'width', 'height', 'length'
    squeezed: list[bool] = []
    shape_squeezed: list[int] = []
    axes_squeezed: list[str] = []
    for size, ax in zip(shape, axes):
        if size > 1 or ax in skip:
            squeezed.append(True)
            shape_squeezed.append(size)
            axes_squeezed.append(ax)
        else:
            squeezed.append(False)
    if len(shape_squeezed) == 0:
        squeezed[-1] = True
        shape_squeezed.append(shape[-1])
        axes_squeezed.append(axes[-1])
    if isinstance(axes, str):
        axes = ''.join(axes_squeezed)
    else:
        axes = tuple(axes_squeezed)
    return (tuple(shape_squeezed), axes, tuple(squeezed))


def transpose_axes(
    image: NDArray[Any],
    axes: str,
    /,
    asaxes: Sequence[str] | None = None,
) -> NDArray[Any]:
    """Return image array with its axes permuted to match specified axes.

    Parameters:
        image:
            Image array to permute.
        axes:
            Character codes for dimensions in image array.
        asaxes:
            Character codes for dimensions in output image array.
            The default is 'CTZYX'.

    Returns:
        Transposed image array.
        A length-1 dimension is added for added dimensions.
        A view of the input array is returned if possible.

    Examples:
        >>> import numpy
        >>> transpose_axes(
        ...     numpy.zeros((2, 3, 4, 5)), 'TYXC', asaxes='CTZYX'
        ... ).shape
        (5, 2, 1, 3, 4)

    """
    if asaxes is None:
        asaxes = 'CTZYX'
    for ax in axes:
        if ax not in asaxes:
            raise ValueError(f'unknown axis {ax}')
    # add missing axes to image
    shape = image.shape
    for ax in reversed(asaxes):
        if ax not in axes:
            axes = ax + axes
            shape = (1,) + shape
    image = image.reshape(shape)
    # transpose axes
    image = image.transpose([axes.index(ax) for ax in asaxes])
    return image


@overload
def reshape_axes(
    axes: str,
    shape: Sequence[int],
    newshape: Sequence[int],
    /,
    unknown: str | None = None,
) -> str: ...


@overload
def reshape_axes(
    axes: Sequence[str],
    shape: Sequence[int],
    newshape: Sequence[int],
    /,
    unknown: str | None = None,
) -> Sequence[str]: ...


def reshape_axes(
    axes: str | Sequence[str],
    shape: Sequence[int],
    newshape: Sequence[int],
    /,
    unknown: str | None = None,
) -> str | Sequence[str]:
    """Return axes matching new shape.

    Parameters:
        axes:
            Character codes for dimensions in `shape`.
        shape:
            Input shape matching `axes`.
        newshape:
            Output shape matching output axes.
            Size must match size of `shape`.
        unknown:
            Character used for new axes in output. The default is 'Q'.

    Returns:
        Character codes for dimensions in `newshape`.

    Examples:
        >>> reshape_axes('YXS', (219, 301, 1), (219, 301))
        'YX'
        >>> reshape_axes('IYX', (12, 219, 301), (3, 4, 219, 1, 301, 1))
        'QQYQXQ'

    """
    shape = tuple(shape)
    newshape = tuple(newshape)
    if len(axes) != len(shape):
        raise ValueError('axes do not match shape')

    size = product(shape)
    newsize = product(newshape)
    if size != newsize:
        raise ValueError(f'cannot reshape {shape} to {newshape}')
    if not axes or not newshape:
        return '' if isinstance(axes, str) else tuple()

    lendiff = max(0, len(shape) - len(newshape))
    if lendiff:
        newshape = newshape + (1,) * lendiff

    i = len(shape) - 1
    prodns = 1
    prods = 1
    result = []
    for ns in newshape[::-1]:
        prodns *= ns
        while i > 0 and shape[i] == 1 and ns != 1:
            i -= 1
        if ns == shape[i] and prodns == prods * shape[i]:
            prods *= shape[i]
            result.append(axes[i])
            i -= 1
        elif unknown:
            result.append(unknown)
        else:
            unknown = 'Q'
            result.append(unknown)

    if isinstance(axes, str):
        axes = ''.join(reversed(result[lendiff:]))
    else:
        axes = tuple(reversed(result[lendiff:]))
    return axes


def order_axes(
    indices: ArrayLike,
    /,
    squeeze: bool = False,
) -> tuple[int, ...]:
    """Return order of axes sorted by variations in indices.

    Parameters:
        indices:
            Multi-dimensional indices of chunks in array.
        squeeze:
            Remove length-1 dimensions of nonvarying axes.

    Returns:
        Order of axes sorted by variations in indices.
        The axis with the least variations in indices is returned first,
        the axis varying fastest is last.

    Examples:
        First axis varies fastest, second axis is squeezed:
        >>> order_axes([(0, 2, 0), (1, 2, 0), (0, 2, 1), (1, 2, 1)], True)
        (2, 0)

    """
    diff = numpy.sum(numpy.abs(numpy.diff(indices, axis=0)), axis=0).tolist()
    order = tuple(sorted(range(len(diff)), key=diff.__getitem__))
    if squeeze:
        order = tuple(i for i in order if diff[i] != 0)
    return order


def check_shape(
    page_shape: Sequence[int], series_shape: Sequence[int]
) -> bool:
    """Return if page and series shapes are compatible."""
    pi = product(page_shape)
    pj = product(series_shape)
    if pi == 0 and pj == 0:
        return True
    if pi == 0 or pj == 0:
        return False
    if pj % pi:
        return False

    series_shape = tuple(reversed(series_shape))
    a = 0
    pi = pj = 1
    for i in reversed(page_shape):
        pi *= i
        # if a == len(series_shape):
        #     return not pj % pi
        for j in series_shape[a:]:
            a += 1
            pj *= j
            if i == j or pi == pj:
                break
            if j == 1:
                continue
            if pj != pi:
                return False
    return True


@overload
def subresolution(
    a: TiffPage, b: TiffPage, /, p: int = 2, n: int = 16
) -> int | None: ...


@overload
def subresolution(
    a: TiffPageSeries, b: TiffPageSeries, /, p: int = 2, n: int = 16
) -> int | None: ...


def subresolution(
    a: TiffPage | TiffPageSeries,
    b: TiffPage | TiffPageSeries,
    /,
    p: int = 2,
    n: int = 16,
) -> int | None:
    """Return level of subresolution of series or page b vs a."""
    if a.axes != b.axes or a.dtype != b.dtype:
        return None
    level = None
    for ax, i, j in zip(a.axes.lower(), a.shape, b.shape):
        if ax in 'xyz':
            if level is None:
                for r in range(n):
                    d = p**r
                    if d > i:
                        return None
                    if abs((i / d) - j) < 1.0:
                        level = r
                        break
                else:
                    return None
            else:
                d = p**level
                if d > i:
                    return None
                if abs((i / d) - j) >= 1.0:
                    return None
        elif i != j:
            return None
    return level


def pyramidize_series(
    series: list[TiffPageSeries], /, isreduced: bool = False
) -> None:
    """Pyramidize list of TiffPageSeries in-place.

    TiffPageSeries that are a subresolution of another TiffPageSeries are
    appended to the other's TiffPageSeries levels and removed from the list.
    Levels are to be ordered by size using the same downsampling factor.
    TiffPageSeries of subifds cannot be pyramid top levels.

    """
    samplingfactors = (2, 3, 4)
    i = 0
    while i < len(series):
        a = series[i]
        p = None
        j = i + 1
        if a.keyframe.is_subifd:
            # subifds cannot be pyramid top levels
            i += 1
            continue
        while j < len(series):
            b = series[j]
            if isreduced and not b.keyframe.is_reduced:
                # pyramid levels must be reduced
                j += 1
                continue  # not a pyramid level
            if p is None:
                for f in samplingfactors:
                    if subresolution(a.levels[-1], b, p=f) == 1:
                        p = f
                        break  # not a pyramid level
                else:
                    j += 1
                    continue  # not a pyramid level
            elif subresolution(a.levels[-1], b, p=p) != 1:
                j += 1
                continue
            a.levels.append(b)
            del series[j]
        i += 1


def stack_pages(
    pages: Sequence[TiffPage | TiffFrame | None],
    /,
    *,
    tiled: TiledSequence | None = None,
    lock: threading.RLock | NullContext | None = None,
    maxworkers: int | None = None,
    out: OutputType = None,
    **kwargs: Any,
) -> NDArray[Any]:
    """Return vertically stacked image arrays from sequence of TIFF pages.

    Parameters:
        pages:
            TIFF pages or frames to stack.
        tiled:
            Organize pages in non-overlapping grid.
        lock:
            Reentrant lock to synchronize seeks and reads from file.
        maxworkers:
            Maximum number of threads to concurrently decode pages or segments.
            By default, use up to :py:attr:`_TIFF.MAXWORKERS` threads.
        out:
            Specifies how image array is returned.
            By default, a new NumPy array is created.
            If a *numpy.ndarray*, a writable array to which the images
            are copied.
            If a string or open file, the file used to create a memory-mapped
            array.
        **kwargs:
            Additional arguments passed to :py:meth:`TiffPage.asarray`.

    """
    npages = len(pages)
    if npages == 0:
        raise ValueError('no pages')

    if npages == 1:
        kwargs['maxworkers'] = maxworkers
        assert pages[0] is not None
        return pages[0].asarray(out=out, **kwargs)

    page0 = next(p.keyframe for p in pages if p is not None)
    assert page0 is not None
    if tiled is None:
        shape = (npages,) + page0.shape
    else:
        shape = tiled.shape
    dtype = page0.dtype
    assert dtype is not None
    out = create_output(out, shape, dtype)

    # TODO: benchmark and optimize this
    if maxworkers is None or maxworkers < 1:
        # auto-detect
        page_maxworkers = page0.maxworkers
        maxworkers = min(npages, TIFF.MAXWORKERS)
        if maxworkers == 1 or page_maxworkers < 1:
            maxworkers = page_maxworkers = 1
        elif npages < 3:
            maxworkers = 1
        elif (
            page_maxworkers <= 2
            and page0.compression == 1
            and page0.fillorder == 1
            and page0.predictor == 1
        ):
            maxworkers = 1
        else:
            page_maxworkers = 1
    elif maxworkers == 1:
        maxworkers = page_maxworkers = 1
    elif npages > maxworkers or page0.maxworkers < 2:
        page_maxworkers = 1
    else:
        page_maxworkers = maxworkers
        maxworkers = 1

    kwargs['maxworkers'] = page_maxworkers

    fh = page0.parent.filehandle
    if lock is None:
        haslock = fh.has_lock
        if not haslock and maxworkers > 1 or page_maxworkers > 1:
            fh.set_lock(True)
        lock = fh.lock
    else:
        haslock = True
    filecache = FileCache(size=max(4, maxworkers), lock=lock)

    if tiled is None:

        def func(
            page: TiffPage | TiffFrame | None,
            index: int,
            out: Any = out,
            filecache: FileCache = filecache,
            kwargs: dict[str, Any] = kwargs,
            /,
        ) -> None:
            # read, decode, and copy page data
            if page is not None:
                filecache.open(page.parent.filehandle)
                page.asarray(lock=lock, out=out[index], **kwargs)
                filecache.close(page.parent.filehandle)

        if maxworkers < 2:
            for index, page in enumerate(pages):
                func(page, index)
        else:
            page0.decode  # init TiffPage.decode function
            with ThreadPoolExecutor(maxworkers) as executor:
                for _ in executor.map(func, pages, range(npages)):
                    pass

    else:
        # TODO: not used or tested

        def func_tiled(
            page: TiffPage | TiffFrame | None,
            index: tuple[int | slice, ...],
            out: Any = out,
            filecache: FileCache = filecache,
            kwargs: dict[str, Any] = kwargs,
            /,
        ) -> None:
            # read, decode, and copy page data
            if page is not None:
                filecache.open(page.parent.filehandle)
                out[index] = page.asarray(lock=lock, **kwargs)
                filecache.close(page.parent.filehandle)

        if maxworkers < 2:
            for index_tiled, page in zip(tiled.slices(), pages):
                func_tiled(page, index_tiled)
        else:
            page0.decode  # init TiffPage.decode function
            with ThreadPoolExecutor(maxworkers) as executor:
                for _ in executor.map(func_tiled, pages, tiled.slices()):
                    pass

    filecache.clear()
    if not haslock:
        fh.set_lock(False)
    return out


def create_output(
    out: OutputType,
    /,
    shape: Sequence[int],
    dtype: DTypeLike,
    *,
    mode: Literal['r+', 'w+', 'r', 'c'] = 'w+',
    suffix: str | None = None,
    fillvalue: int | float | None = 0,
) -> NDArray[Any] | numpy.memmap[Any, Any]:
    """Return NumPy array where images of shape and dtype can be copied.

    Parameters:
        out:
            Specifies kind of array to return:

                `None`:
                    A new array of shape and dtype is created and returned.
                `numpy.ndarray`:
                    An existing, writable array compatible with `dtype` and
                    `shape`. A view of the array is returned.
                `'memmap'` or `'memmap:tempdir'`:
                    A memory-map to an array stored in a temporary binary file
                    on disk is created and returned.
                `str` or open file:
                    File name or file object used to create a memory-map
                    to an array stored in a binary file on disk.
                    The memory-mapped array is returned.
        shape:
            Shape of NumPy array to return.
        dtype:
            Data type of NumPy array to return.
        suffix:
            Suffix of `NamedTemporaryFile` if `out` is 'memmap'.
            The default suffix is 'memmap'.
        fillvalue:
            Value to initialize newly created arrays.
            If *None*, return an uninitialized array.

    """
    shape = tuple(shape)
    if out is None:
        if fillvalue is None:
            return numpy.empty(shape, dtype)
        if fillvalue:
            out = numpy.empty(shape, dtype)
            out[:] = fillvalue
            return out
        return numpy.zeros(shape, dtype)
    if isinstance(out, numpy.ndarray):
        if product(shape) != product(out.shape):
            raise ValueError('incompatible output shape')
        if not numpy.can_cast(dtype, out.dtype):
            raise ValueError('incompatible output dtype')
        return out.reshape(shape)
    if isinstance(out, str) and out[:6] == 'memmap':
        import tempfile

        tempdir = out[7:] if len(out) > 7 else None
        if suffix is None:
            suffix = '.memmap'
        with tempfile.NamedTemporaryFile(dir=tempdir, suffix=suffix) as fh:
            out = numpy.memmap(fh, shape=shape, dtype=dtype, mode=mode)
            if fillvalue:
                out[:] = fillvalue
            return out
    out = numpy.memmap(out, shape=shape, dtype=dtype, mode=mode)
    if fillvalue:
        out[:] = fillvalue
    return out


def matlabstr2py(matlabstr: str, /) -> Any:
    r"""Return Python object from Matlab string representation.

    Use to access ScanImage metadata.

    Parameters:
        matlabstr: String representation of Matlab objects.

    Returns:
        Matlab structures are returned as `dict`.
        Matlab arrays or cells are returned as `lists`.
        Other Matlab objects are returned as `str`, `bool`, `int`, or `float`.

    Examples:
        >>> matlabstr2py('1')
        1
        >>> matlabstr2py("['x y z' true false; 1 2.0 -3e4; NaN Inf @class]")
        [['x y z', True, False], [1, 2.0, -30000.0], [nan, inf, '@class']]
        >>> d = matlabstr2py(
        ...     "SI.hChannels.channelType = {'stripe' 'stripe'}\n"
        ...     "SI.hChannels.channelsActive = 2"
        ... )
        >>> d['SI.hChannels.channelType']
        ['stripe', 'stripe']

    """
    # TODO: handle invalid input
    # TODO: review unboxing of multidimensional arrays

    def lex(s: str, /) -> list[str]:
        # return sequence of tokens from Matlab string representation
        tokens = ['[']
        while True:
            t, i = next_token(s)
            if t is None:
                break
            if t == ';':
                tokens.extend((']', '['))
            elif t == '[':
                tokens.extend(('[', '['))
            elif t == ']':
                tokens.extend((']', ']'))
            else:
                tokens.append(t)
            s = s[i:]
        tokens.append(']')
        return tokens

    def next_token(s: str, /) -> tuple[str | None, int]:
        # return next token in Matlab string
        length = len(s)
        if length == 0:
            return None, 0
        i = 0
        while i < length and s[i] == ' ':
            i += 1
        if i == length:
            return None, i
        if s[i] in '{[;]}':
            return s[i], i + 1
        if s[i] == "'":
            j = i + 1
            while j < length and s[j] != "'":
                j += 1
            return s[i : j + 1], j + 1
        if s[i] == '<':
            j = i + 1
            while j < length and s[j] != '>':
                j += 1
            return s[i : j + 1], j + 1
        j = i
        while j < length and s[j] not in ' {[;]}':
            j += 1
        return s[i:j], j

    def value(s: str, fail: bool = False, /) -> Any:
        # return Python value of token
        s = s.strip()
        if not s:
            return s
        if len(s) == 1:
            try:
                return int(s)
            except Exception as exc:
                if fail:
                    raise ValueError from exc
                return s
        if s[0] == "'":
            if fail and s[-1] != "'" or "'" in s[1:-1]:
                raise ValueError
            return s[1:-1]
        if s[0] == '<':
            if fail and s[-1] != '>' or '<' in s[1:-1]:
                raise ValueError
            return s
        if fail and any(i in s for i in " ';[]{}"):
            raise ValueError
        if s[0] == '@':
            return s
        if s in {'true', 'True'}:
            return True
        if s in {'false', 'False'}:
            return False
        if s[:6] == 'zeros(':
            return numpy.zeros([int(i) for i in s[6:-1].split(',')]).tolist()
        if s[:5] == 'ones(':
            return numpy.ones([int(i) for i in s[5:-1].split(',')]).tolist()
        if '.' in s or 'e' in s:
            try:
                return float(s)
            except Exception:
                pass
        try:
            return int(s)
        except Exception:
            pass
        try:
            return float(s)  # nan, inf
        except Exception as exc:
            if fail:
                raise ValueError from exc
        return s

    def parse(s: str, /) -> Any:
        # return Python value from string representation of Matlab value
        s = s.strip()
        try:
            return value(s, True)
        except ValueError:
            pass
        result: list[Any]
        addto: list[Any]
        result = addto = []
        levels = [addto]
        for t in lex(s):
            if t in '[{':
                addto = []
                levels.append(addto)
            elif t in ']}':
                x = levels.pop()
                addto = levels[-1]
                if len(x) == 1 and isinstance(x[0], (list, str)):
                    addto.append(x[0])
                else:
                    addto.append(x)
            else:
                addto.append(value(t))
        if len(result) == 1 and isinstance(result[0], (list, str)):
            return result[0]
        return result

    if '\r' in matlabstr or '\n' in matlabstr:
        # structure
        d = {}
        for line in matlabstr.splitlines():
            line = line.strip()
            if not line or line[0] == '%':
                continue
            k, v = line.split('=', 1)
            k = k.strip()
            if any(c in k for c in " ';[]{}<>"):
                continue
            d[k] = parse(v)
        return d
    return parse(matlabstr)


def strptime(datetime_string: str, format: str | None = None, /) -> DateTime:
    """Return datetime corresponding to date string using common formats.

    Parameters:
        datetime_string:
            String representation of date and time.
        format:
            Format of `datetime_string`.
            By default, several datetime formats commonly found in TIFF files
            are parsed.

    Raises:
        ValueError: `datetime_string` does not match any format.

    Examples:
        >>> strptime('2022:08:01 22:23:24')
        datetime.datetime(2022, 8, 1, 22, 23, 24)

    """
    formats = {
        '%Y:%m:%d %H:%M:%S': 1,  # TIFF6 specification
        '%Y%m%d %H:%M:%S.%f': 2,  # MetaSeries
        '%Y-%m-%dT%H %M %S.%f': 3,  # Pilatus
        '%Y-%m-%dT%H:%M:%S.%f': 4,  # ISO
        '%Y-%m-%dT%H:%M:%S': 5,  # ISO, microsecond is 0
        '%Y:%m:%d %H:%M:%S.%f': 6,
        '%d/%m/%Y %H:%M:%S': 7,
        '%d/%m/%Y %H:%M:%S.%f': 8,
        '%m/%d/%Y %I:%M:%S %p': 9,
        '%m/%d/%Y %I:%M:%S.%f %p': 10,
        '%Y%m%d %H:%M:%S': 11,
        '%Y/%m/%d %H:%M:%S': 12,
        '%Y/%m/%d %H:%M:%S.%f': 13,
        '%Y-%m-%dT%H:%M:%S%z': 14,
        '%Y-%m-%dT%H:%M:%S.%f%z': 15,
    }
    if format is not None:
        formats[format] = 0  # highest priority; replaces existing key if any
    for format, _ in sorted(formats.items(), key=lambda item: item[1]):
        try:
            return DateTime.strptime(datetime_string, format)
        except ValueError:
            pass
    raise ValueError(
        f'time data {datetime_string!r} does not match any format'
    )


@overload
def stripnull(
    string: bytes, /, null: bytes | None = None, *, first: bool = True
) -> bytes: ...


@overload
def stripnull(
    string: str, /, null: str | None = None, *, first: bool = True
) -> str: ...


def stripnull(
    string: str | bytes,
    /,
    null: str | bytes | None = None,
    *,
    first: bool = True,
) -> str | bytes:
    r"""Return string truncated at first null character.

    Use to clean NULL terminated C strings.

    >>> stripnull(b'bytes\x00\x00')
    b'bytes'
    >>> stripnull(b'bytes\x00bytes\x00\x00', first=False)
    b'bytes\x00bytes'
    >>> stripnull('string\x00')
    'string'

    """
    # TODO: enable deprecation warning
    # warnings.warn(
    #     '<tifffile.stripnull is deprecated since 2025.3.18',
    #     DeprecationWarning,
    #     stacklevel=2,
    # )
    if null is None:
        if isinstance(string, bytes):
            null = b'\x00'
        else:
            null = '\0'
    if first:
        i = string.find(null)  # type: ignore[arg-type]
        return string if i < 0 else string[:i]
    return string.rstrip(null)  # type: ignore[arg-type]


def stripascii(string: bytes, /) -> bytes:
    r"""Return string truncated at last byte that is 7-bit ASCII.

    Use to clean NULL separated and terminated TIFF strings.

    >>> stripascii(b'string\x00string\n\x01\x00')
    b'string\x00string\n'
    >>> stripascii(b'\x00')
    b''

    """
    # TODO: pythonize this
    i = len(string)
    while i:
        i -= 1
        if 8 < string[i] < 127:
            break
    else:
        i = -1
    return string[: i + 1]


@overload
def asbool(
    value: str,
    /,
    true: Sequence[str] | None = None,
    false: Sequence[str] | None = None,
) -> bool: ...


@overload
def asbool(
    value: bytes,
    /,
    true: Sequence[bytes] | None = None,
    false: Sequence[bytes] | None = None,
) -> bool: ...


def asbool(
    value: str | bytes,
    /,
    true: Sequence[str | bytes] | None = None,
    false: Sequence[str | bytes] | None = None,
) -> bool | bytes:
    """Return string as bool if possible, else raise TypeError.

    >>> asbool(b' False ')
    False
    >>> asbool('ON', ['on'], ['off'])
    True

    """
    value = value.strip().lower()
    isbytes = False
    if true is None:
        if isinstance(value, bytes):
            if value == b'true':
                return True
            isbytes = True
        elif value == 'true':
            return True
    elif value in true:
        return True
    if false is None:
        if isbytes or isinstance(value, bytes):
            if value == b'false':
                return False
        elif value == 'false':
            return False
    elif value in false:
        return False
    raise TypeError


def astype(value: Any, /, types: Sequence[Any] | None = None) -> Any:
    """Return argument as one of types if possible.

    >>> astype('42')
    42
    >>> astype('3.14')
    3.14
    >>> astype('True')
    True
    >>> astype(b'Neee-Wom')
    'Neee-Wom'

    """
    if types is None:
        types = int, float, asbool, bytes2str
    for typ in types:
        try:
            return typ(value)
        except (ValueError, AttributeError, TypeError, UnicodeEncodeError):
            pass
    return value


def rational(arg: float | tuple[int, int], /) -> tuple[int, int]:
    """Return rational numerator and denominator from float or two integers."""
    from fractions import Fraction

    if isinstance(arg, Sequence):
        f = Fraction(arg[0], arg[1])
    else:
        f = Fraction.from_float(arg)

    numerator, denominator = f.as_integer_ratio()
    if numerator > 4294967295 or denominator > 4294967295:
        s = 4294967295 / max(numerator, denominator)
        numerator = round(numerator * s)
        denominator = round(denominator * s)
    return numerator, denominator


def unique_strings(strings: Iterator[str], /) -> Iterator[str]:
    """Return iterator over unique strings.

    >>> list(unique_strings(iter(('a', 'b', 'a'))))
    ['a', 'b', 'a2']

    """
    known = set()
    for i, string in enumerate(strings):
        if string in known:
            string += str(i)
        known.add(string)
        yield string


def format_size(size: int | float, /, threshold: int | float = 1536) -> str:
    """Return file size as string from byte size.

    >>> format_size(1234)
    '1234 B'
    >>> format_size(12345678901)
    '11.50 GiB'

    """
    if size < threshold:
        return f'{size} B'
    for unit in ('KiB', 'MiB', 'GiB', 'TiB', 'PiB'):
        size /= 1024.0
        if size < threshold:
            return f'{size:.2f} {unit}'
    return 'ginormous'


def identityfunc(arg: Any, /, *args: Any, **kwargs: Any) -> Any:
    """Single argument identity function.

    >>> identityfunc('arg')
    'arg'

    """
    return arg


def nullfunc(*args: Any, **kwargs: Any) -> None:
    """Null function.

    >>> nullfunc('arg', kwarg='kwarg')

    """
    return


def sequence(value: Any, /) -> Sequence[Any]:
    """Return tuple containing value if value is not tuple or list.

    >>> sequence(1)
    (1,)
    >>> sequence([1])
    [1]
    >>> sequence('ab')
    ('ab',)

    """
    return value if isinstance(value, (tuple, list)) else (value,)


def product(iterable: Iterable[int], /) -> int:
    """Return product of integers.

    Equivalent of ``math.prod(iterable)``, but multiplying NumPy integers
    does not overflow.

    >>> product([2**8, 2**30])
    274877906944
    >>> product([])
    1

    """
    prod = 1
    for i in iterable:
        prod *= int(i)
    return prod


def peek_iterator(iterator: Iterator[Any], /) -> tuple[Any, Iterator[Any]]:
    """Return first item of iterator and iterator.

    >>> first, it = peek_iterator(iter((0, 1, 2)))
    >>> first
    0
    >>> list(it)
    [0, 1, 2]

    """
    first = next(iterator)

    def newiter(
        first: Any = first, iterator: Iterator[Any] = iterator
    ) -> Iterator[Any]:
        yield first
        yield from iterator

    return first, newiter()


def natural_sorted(iterable: Iterable[str], /) -> list[str]:
    """Return human-sorted list of strings.

    Use to sort file names.

    >>> natural_sorted(['f1', 'f2', 'f10'])
    ['f1', 'f2', 'f10']

    """

    def sortkey(x: str, /) -> list[int | str]:
        return [(int(c) if c.isdigit() else c) for c in re.split(numbers, x)]

    numbers = re.compile(r'(\d+)')
    return sorted(iterable, key=sortkey)


def epics_datetime(sec: int, nsec: int, /) -> DateTime:
    """Return datetime object from epicsTSSec and epicsTSNsec tag values.

    >>> epics_datetime(802117916, 103746502)
    datetime.datetime(2015, 6, 2, 11, 31, 56, 103746)

    """
    return DateTime.fromtimestamp(sec + 631152000 + nsec / 1e9)


def excel_datetime(timestamp: float, epoch: int | None = None, /) -> DateTime:
    """Return datetime object from timestamp in Excel serial format.

    Use to convert LSM time stamps.

    >>> excel_datetime(40237.029999999795)
    datetime.datetime(2010, 2, 28, 0, 43, 11, 999982)

    """
    if epoch is None:
        epoch = 693594
    return DateTime.fromordinal(epoch) + TimeDelta(timestamp)


def julian_datetime(julianday: int, millisecond: int = 0, /) -> DateTime:
    """Return datetime from days since 1/1/4713 BC and ms since midnight.

    Convert Julian dates according to MetaMorph.

    >>> julian_datetime(2451576, 54362783)
    datetime.datetime(2000, 2, 2, 15, 6, 2, 783000)

    """
    if julianday <= 1721423:
        # return DateTime.min  # ?
        raise ValueError(f'no datetime before year 1 ({julianday=})')

    a = julianday + 1
    if a > 2299160:
        alpha = math.trunc((a - 1867216.25) / 36524.25)
        a += 1 + alpha - alpha // 4
    b = a + (1524 if a > 1721423 else 1158)
    c = math.trunc((b - 122.1) / 365.25)
    d = math.trunc(365.25 * c)
    e = math.trunc((b - d) / 30.6001)

    day = b - d - math.trunc(30.6001 * e)
    month = e - (1 if e < 13.5 else 13)
    year = c - (4716 if month > 2.5 else 4715)

    hour, millisecond = divmod(millisecond, 1000 * 60 * 60)
    minute, millisecond = divmod(millisecond, 1000 * 60)
    second, millisecond = divmod(millisecond, 1000)

    return DateTime(year, month, day, hour, minute, second, millisecond * 1000)


def byteorder_isnative(byteorder: str, /) -> bool:
    """Return if byteorder matches system's byteorder.

    >>> byteorder_isnative('=')
    True

    """
    if byteorder in {'=', sys.byteorder}:
        return True
    keys = {'big': '>', 'little': '<'}
    return keys.get(byteorder, byteorder) == keys[sys.byteorder]


def byteorder_compare(byteorder: str, other: str, /) -> bool:
    """Return if byteorders match.

    >>> byteorder_compare('<', '<')
    True
    >>> byteorder_compare('>', '<')
    False

    """
    if byteorder in {other, '|'} or other == '|':
        return True
    if byteorder == '=':
        byteorder = {'big': '>', 'little': '<'}[sys.byteorder]
    elif other == '=':
        other = {'big': '>', 'little': '<'}[sys.byteorder]
    return byteorder == other


def recarray2dict(recarray: numpy.recarray[Any, Any], /) -> dict[str, Any]:
    """Return numpy.recarray as dictionary.

    >>> r = numpy.array(
    ...     [(1.0, 2, 'a'), (3.0, 4, 'bc')],
    ...     dtype=[('x', '<f4'), ('y', '<i4'), ('s', 'S2')],
    ... )
    >>> recarray2dict(r)
    {'x': [1.0, 3.0], 'y': [2, 4], 's': ['a', 'bc']}
    >>> recarray2dict(r[1])
    {'x': 3.0, 'y': 4, 's': 'bc'}

    """
    # TODO: subarrays
    value: Any
    result = {}
    for descr in recarray.dtype.descr:
        name, dtype = descr[:2]
        value = recarray[name]
        if value.ndim == 0:
            value = value.tolist()
            if dtype[1] == 'S':
                value = bytes2str(value)
        elif value.ndim == 1:
            value = value.tolist()
            if dtype[1] == 'S':
                value = [bytes2str(v) for v in value]
        result[name] = value
    return result


def xml2dict(
    xml: str,
    /,
    *,
    sanitize: bool = True,
    prefix: tuple[str, str] | None = None,
    sep: str = ',',
) -> dict[str, Any]:
    """Return XML as dictionary.

    Parameters:
        xml: XML data to convert.
        sanitize: Remove prefix from from etree Element.
        prefix: Prefixes for dictionary keys.
        sep: Sequence separator.

    Examples:
        >>> xml2dict(
        ...     '<?xml version="1.0" ?><root attr="name"><key>1</key></root>'
        ... )
        {'root': {'key': 1, 'attr': 'name'}}
        >>> xml2dict('<level1><level2>3.5322,-3.14</level2></level1>')
        {'level1': {'level2': (3.5322, -3.14)}}

    """
    try:
        from defusedxml import ElementTree as etree
    except ImportError:
        from xml.etree import ElementTree as etree

    at, tx = prefix if prefix else ('', '')

    def astype(value: Any, /) -> Any:
        # return string value as int, float, bool, tuple, or unchanged
        if not isinstance(value, str):
            return value
        if sep and sep in value:
            # sequence of numbers?
            values = []
            for val in value.split(sep):
                v = astype(val)
                if isinstance(v, str):
                    return value
                values.append(v)
            return tuple(values)
        for t in (int, float, asbool):
            try:
                return t(value)
            except (TypeError, ValueError):
                pass
        return value

    def etree2dict(t: Any, /) -> dict[str, Any]:
        # adapted from https://stackoverflow.com/a/10077069/453463
        key = t.tag
        if sanitize:
            key = key.rsplit('}', 1)[-1]
        d: dict[str, Any] = {key: {} if t.attrib else None}
        children = list(t)
        if children:
            dd = collections.defaultdict(list)
            for dc in map(etree2dict, children):
                for k, v in dc.items():
                    dd[k].append(astype(v))
            d = {
                key: {
                    k: astype(v[0]) if len(v) == 1 else astype(v)
                    for k, v in dd.items()
                }
            }
        if t.attrib:
            d[key].update((at + k, astype(v)) for k, v in t.attrib.items())
        if t.text:
            text = t.text.strip()
            if children or t.attrib:
                if text:
                    d[key][tx + 'value'] = astype(text)
            else:
                d[key] = astype(text)
        return d

    return etree2dict(etree.fromstring(xml))


def hexdump(
    data: bytes,
    /,
    *,
    width: int = 75,
    height: int = 24,
    snipat: int | float | None = 0.75,
    modulo: int = 2,
    ellipsis: str | None = None,
) -> str:
    """Return hexdump representation of bytes.

    Parameters:
        data:
            Bytes to represent as hexdump.
        width:
            Maximum width of hexdump.
        height:
            Maximum number of lines of hexdump.
        snipat:
            Approximate position at which to split long hexdump.
        modulo:
            Number of bytes represented in line of hexdump are modulus
            of this value.
        ellipsis:
            Characters to insert for snipped content of long hexdump.
            The default is '...'.

    Examples:
        >>> import binascii
        >>> hexdump(binascii.unhexlify('49492a00080000000e00fe0004000100'))
        '49 49 2a 00 08 00 00 00 0e 00 fe 00 04 00 01 00 II*.............'

    """
    size = len(data)
    if size < 1 or width < 2 or height < 1:
        return ''
    if height == 1:
        addr = b''
        bytesperline = min(
            modulo * (((width - len(addr)) // 4) // modulo), size
        )
        if bytesperline < 1:
            return ''
        nlines = 1
    else:
        addr = b'%%0%ix: ' % len(b'%x' % size)
        bytesperline = min(
            modulo * (((width - len(addr % 1)) // 4) // modulo), size
        )
        if bytesperline < 1:
            return ''
        width = 3 * bytesperline + len(addr % 1)
        nlines = (size - 1) // bytesperline + 1

    if snipat is None or snipat == 1:
        snipat = height
    elif 0 < abs(snipat) < 1:
        snipat = int(math.floor(height * snipat))
    if snipat < 0:
        snipat += height
    assert isinstance(snipat, int)

    blocks: list[tuple[int, bytes | None]]

    if height == 1 or nlines == 1:
        blocks = [(0, data[:bytesperline])]
        addr = b''
        height = 1
        width = 3 * bytesperline
    elif not height or nlines <= height:
        blocks = [(0, data)]
    elif snipat <= 0:
        start = bytesperline * (nlines - height)
        blocks = [(start, data[start:])]  # (start, None)
    elif snipat >= height or height < 3:
        end = bytesperline * height
        blocks = [(0, data[:end])]  # (end, None)
    else:
        end1 = bytesperline * snipat
        end2 = bytesperline * (height - snipat - 2)
        if size % bytesperline:
            end2 += size % bytesperline
        else:
            end2 += bytesperline
        blocks = [
            (0, data[:end1]),
            (size - end1 - end2, None),
            (size - end2, data[size - end2 :]),
        ]

    if ellipsis is None:
        if addr and bytesperline > 3:
            elps = b' ' * (len(addr % 1) + bytesperline // 2 * 3 - 2)
            elps += b'...'
        else:
            elps = b'...'
    else:
        elps = ellipsis.encode('cp1252')

    result = []
    for start, bstr in blocks:
        if bstr is None:
            result.append(elps)  # 'skip %i bytes' % start)
            continue
        hexstr = binascii.hexlify(bstr)
        strstr = re.sub(br'[^\x20-\x7f]', b'.', bstr)
        for i in range(0, len(bstr), bytesperline):
            h = hexstr[2 * i : 2 * i + bytesperline * 2]
            r = (addr % (i + start)) if height > 1 else addr
            r += b' '.join(h[i : i + 2] for i in range(0, 2 * bytesperline, 2))
            r += b' ' * (width - len(r))
            r += strstr[i : i + bytesperline]
            result.append(r)
    return b'\n'.join(result).decode('ascii')


def isprintable(string: str | bytes, /) -> bool:
    r"""Return if all characters in string are printable.

    >>> isprintable('abc')
    True
    >>> isprintable(b'\01')
    False

    """
    string = string.strip()
    if not string:
        return True
    try:
        return string.isprintable()  # type: ignore[union-attr]
    except Exception:
        pass
    try:
        return string.decode().isprintable()  # type: ignore[union-attr]
    except Exception:
        pass
    return False


def clean_whitespace(string: str, /, compact: bool = False) -> str:
    r"""Return string with compressed whitespace.

    >>> clean_whitespace('  a  \n\n  b ')
    'a\n b'

    """
    string = (
        string.replace('\r\n', '\n')
        .replace('\r', '\n')
        .replace('\n\n', '\n')
        .replace('\t', ' ')
        .replace('  ', ' ')
        .replace('  ', ' ')
        .replace(' \n', '\n')
    )
    if compact:
        string = (
            string.replace('\n', ' ')
            .replace('[ ', '[')
            .replace('  ', ' ')
            .replace('  ', ' ')
            .replace('  ', ' ')
        )
    return string.strip()


def indent(*args: Any) -> str:
    """Return joined string representations of objects with indented lines.

    >>> print(indent('Title:', 'Text'))
    Title:
      Text

    """
    text = '\n'.join(str(arg) for arg in args)
    return '\n'.join(
        ('  ' + line if line else line) for line in text.splitlines() if line
    )[2:]


def pformat_xml(xml: str | bytes, /) -> str:
    """Return pretty formatted XML."""
    try:
        from lxml import etree

        if not isinstance(xml, bytes):
            xml = xml.encode()
        tree = etree.parse(io.BytesIO(xml))
        xml = etree.tostring(
            tree,
            pretty_print=True,
            xml_declaration=True,
            encoding=tree.docinfo.encoding,
        )
        assert isinstance(xml, bytes)
        xml = bytes2str(xml)
    except Exception:
        if isinstance(xml, bytes):
            xml = bytes2str(xml)
        xml = xml.replace('><', '>\n<')
    return xml.replace('  ', ' ').replace('\t', ' ')


def pformat(
    arg: Any,
    /,
    *,
    height: int | None = 24,
    width: int | None = 79,
    linewidth: int | None = 288,
    compact: bool = True,
) -> str:
    """Return pretty formatted representation of object as string.

    Whitespace might be altered. Long lines are cut off.

    """
    if height is None or height < 1:
        height = 1024
    if width is None or width < 1:
        width = 256
    if linewidth is None or linewidth < 1:
        linewidth = width

    npopt = numpy.get_printoptions()
    numpy.set_printoptions(threshold=100, linewidth=width)

    if isinstance(arg, bytes):
        if arg[:5].lower() == b'<?xml' or arg[-4:] == b'OME>':
            arg = bytes2str(arg)

    if isinstance(arg, bytes):
        if isprintable(arg):
            arg = bytes2str(arg)
            arg = clean_whitespace(arg)
        else:
            numpy.set_printoptions(**npopt)
            return hexdump(arg, width=width, height=height, modulo=1)
        arg = arg.rstrip()
    elif isinstance(arg, str):
        if arg[:5].lower() == '<?xml' or arg[-4:] == 'OME>':
            arg = arg[: 4 * width] if height == 1 else pformat_xml(arg)
        # too slow
        # else:
        #    import textwrap
        #    return '\n'.join(
        #        textwrap.wrap(arg, width=width, max_lines=height, tabsize=2)
        #    )
        arg = arg.rstrip()
    elif isinstance(arg, numpy.record):
        arg = arg.pprint()
    else:
        import pprint

        arg = pprint.pformat(arg, width=width, compact=compact)

    numpy.set_printoptions(**npopt)

    if height == 1:
        arg = arg[: width * width]
        arg = clean_whitespace(arg, compact=True)
        return arg[:linewidth]

    argl = list(arg.splitlines())
    if len(argl) > height:
        arg = '\n'.join(
            line[:linewidth]
            for line in argl[: height // 2] + ['...'] + argl[-height // 2 :]
        )
    else:
        arg = '\n'.join(line[:linewidth] for line in argl[:height])
    return arg


def snipstr(
    string: str,
    /,
    width: int = 79,
    *,
    snipat: int | float | None = None,
    ellipsis: str | None = None,
) -> str:
    """Return string cut to specified length.

    Parameters:
        string:
            String to snip.
        width:
            Maximum length of returned string.
        snipat:
            Approximate position at which to split long strings.
            The default is 0.5.
        ellipsis:
            Characters to insert between splits of long strings.
            The default is '...'.

    Examples:
        >>> snipstr('abcdefghijklmnop', 8)
        'abc...op'

    """
    if snipat is None:
        snipat = 0.5
    if ellipsis is None:
        if isinstance(string, bytes):  # type: ignore[unreachable]
            ellipsis = b'...'
        else:
            ellipsis = '\u2026'
    esize = len(ellipsis)

    splitlines = string.splitlines()
    # TODO: finish and test multiline snip

    result = []
    for line in splitlines:
        if line is None:
            result.append(ellipsis)
            continue
        linelen = len(line)
        if linelen <= width:
            result.append(string)
            continue

        if snipat is None or snipat == 1:
            split = linelen
        elif 0 < abs(snipat) < 1:
            split = int(math.floor(linelen * snipat))
        else:
            split = int(snipat)

        if split < 0:
            split += linelen
            split = max(split, 0)

        if esize == 0 or width < esize + 1:
            if split <= 0:
                result.append(string[-width:])
            else:
                result.append(string[:width])
        elif split <= 0:
            result.append(ellipsis + string[esize - width :])
        elif split >= linelen or width < esize + 4:
            result.append(string[: width - esize] + ellipsis)
        else:
            splitlen = linelen - width + esize
            end1 = split - splitlen // 2
            end2 = end1 + splitlen
            result.append(string[:end1] + ellipsis + string[end2:])

    if isinstance(string, bytes):  # type: ignore[unreachable]
        return b'\n'.join(result)
    return '\n'.join(result)


def enumstr(enum: Any, /) -> str:
    """Return short string representation of Enum member.

    >>> enumstr(PHOTOMETRIC.RGB)
    'RGB'

    """
    name = enum.name
    if name is None:
        name = str(enum)
    return name


def enumarg(enum: type[enum.IntEnum], arg: Any, /) -> enum.IntEnum:
    """Return enum member from its name or value.

    Parameters:
        enum: Type of IntEnum.
        arg: Name or value of enum member.

    Returns:
        Enum member matching name or value.

    Raises:
        ValueError: No enum member matches name or value.

    Examples:
        >>> enumarg(PHOTOMETRIC, 2)
        <PHOTOMETRIC.RGB: 2>
        >>> enumarg(PHOTOMETRIC, 'RGB')
        <PHOTOMETRIC.RGB: 2>

    """
    try:
        return enum(arg)
    except Exception:
        try:
            return enum[arg.upper()]
        except Exception as exc:
            raise ValueError(f'invalid argument {arg!r}') from exc


def parse_kwargs(
    kwargs: dict[str, Any], /, *keys: str, **keyvalues: Any
) -> dict[str, Any]:
    """Return dict with keys from keys|keyvals and values from kwargs|keyvals.

    Existing keys are deleted from `kwargs`.

    >>> kwargs = {'one': 1, 'two': 2, 'four': 4}
    >>> kwargs2 = parse_kwargs(kwargs, 'two', 'three', four=None, five=5)
    >>> kwargs == {'one': 1}
    True
    >>> kwargs2 == {'two': 2, 'four': 4, 'five': 5}
    True

    """
    result = {}
    for key in keys:
        if key in kwargs:
            result[key] = kwargs[key]
            del kwargs[key]
    for key, value in keyvalues.items():
        if key in kwargs:
            result[key] = kwargs[key]
            del kwargs[key]
        else:
            result[key] = value
    return result


def update_kwargs(kwargs: dict[str, Any], /, **keyvalues: Any) -> None:
    """Update dict with keys and values if keys do not already exist.

    >>> kwargs = {'one': 1}
    >>> update_kwargs(kwargs, one=None, two=2)
    >>> kwargs == {'one': 1, 'two': 2}
    True

    """
    for key, value in keyvalues.items():
        if key not in kwargs:
            kwargs[key] = value


def kwargs_notnone(**kwargs: Any) -> dict[str, Any]:
    """Return dict of kwargs which values are not None.

    >>> kwargs_notnone(one=1, none=None)
    {'one': 1}

    """
    return dict(item for item in kwargs.items() if item[1] is not None)


def logger() -> logging.Logger:
    """Return logger for tifffile module."""
    return logging.getLogger('tifffile')


def validate_jhove(
    filename: str,
    /,
    jhove: str | None = None,
    ignore: Collection[str] | None = None,
) -> None:
    """Validate TIFF file with ``jhove -m TIFF-hul``.

    JHOVE does not support the BigTIFF format, more than 50 IFDs, and
    many TIFF extensions.

    Parameters:
        filename:
            Name of TIFF file to validate.
        jhove:
            Path of jhove app. The default is 'jhove'.
        ignore:
            Jhove error message to ignore.

    Raises:
        ValueError:
            Jhove printed error message and did not contain one of strings
            in `ignore`.

    References:
        - `JHOVE TIFF-hul Module <http://jhove.sourceforge.net/tiff-hul.html>`_

    """
    import subprocess

    if ignore is None:
        ignore = {'More than 50 IFDs', 'Predictor value out of range'}
    if jhove is None:
        jhove = 'jhove'
    out = subprocess.check_output([jhove, filename, '-m', 'TIFF-hul'])
    if b'ErrorMessage: ' in out:
        for line in out.splitlines():
            line = line.strip()
            if line.startswith(b'ErrorMessage: '):
                error = line[14:].decode()
                for i in ignore:
                    if i in error:
                        break
                else:
                    raise ValueError(error)
                break


def tiffcomment(
    arg: str | os.PathLike[Any] | FileHandle | IO[bytes],
    /,
    comment: str | bytes | None = None,
    pageindex: int | None = None,
    tagcode: int | str | None = None,
) -> str | None:
    """Return or replace ImageDescription value in first page of TIFF file.

    Parameters:
        arg:
            Specifies TIFF file to open.
        comment:
            7-bit ASCII string or bytes to replace existing tag value.
            The existing value is zeroed.
        pageindex:
            Index of page which ImageDescription tag value to
            read or replace. The default is 0.
        tagcode:
            Code of tag which value to read or replace.
            The default is 270 (ImageDescription).

    Returns:
        None, if `comment` is specified. Else, the current value of the
        specified tag in the specified page.


    """
    if pageindex is None:
        pageindex = 0
    if tagcode is None:
        tagcode = 270
    mode: Any = None if comment is None else 'r+'
    with TiffFile(arg, mode=mode) as tif:
        page = tif.pages[pageindex]
        if not isinstance(page, TiffPage):
            raise IndexError(f'TiffPage {pageindex} not found')
        tag = page.tags.get(tagcode, None)
        if tag is None:
            raise ValueError(f'no {TIFF.TAGS[tagcode]} tag found')
        if comment is None:
            return tag.value
        tag.overwrite(comment)
        return None


def tiff2fsspec(
    filename: str | os.PathLike[Any],
    /,
    url: str,
    *,
    out: str | None = None,
    key: int | None = None,
    series: int | None = None,
    level: int | None = None,
    chunkmode: CHUNKMODE | int | str | None = None,
    fillvalue: int | float | None = None,
    zattrs: dict[str, Any] | None = None,
    squeeze: bool | None = None,
    groupname: str | None = None,
    version: int | None = None,
) -> None:
    """Write fsspec ReferenceFileSystem in JSON format for data in TIFF file.

    By default, the first series, including all levels, is exported.

    Parameters:
        filename:
            Name of TIFF file to reference.
        url:
            Remote location of TIFF file without file name(s).
        out:
            Name of output JSON file.
            The default is the `filename` with a '.json' extension.
        key, series, level, chunkmode, fillvalue, zattrs, squeeze:
            Passed to :py:meth:`TiffFile.aszarr`.
        groupname, version:
            Passed to :py:meth:`ZarrTiffStore.write_fsspec`.

    """
    if out is None:
        out = os.fspath(filename) + '.json'
    with TiffFile(filename) as tif:
        store: ZarrTiffStore
        with tif.aszarr(
            key=key,
            series=series,
            level=level,
            chunkmode=chunkmode,
            fillvalue=fillvalue,
            zattrs=zattrs,
            squeeze=squeeze,
        ) as store:
            store.write_fsspec(out, url, groupname=groupname, version=version)


def lsm2bin(
    lsmfile: str,
    /,
    binfile: str | None = None,
    *,
    tile: tuple[int, int] | None = None,
    verbose: bool = True,
) -> None:
    """Convert [MP]TZCYX LSM file to series of BIN files.

    One BIN file containing 'ZCYX' data is created for each position, time,
    and tile. The position, time, and tile indices are encoded at the end
    of the filenames.

    Parameters:
        lsmfile:
            Name of LSM file to convert.
        binfile:
            Common name of output BIN files.
            The default is the name of the LSM file without extension.
        tile:
            Y and X dimension sizes of BIN files.
            The default is (256, 256).
        verbose:
            Print status of conversion.

    """
    prints: Any = print if verbose else nullfunc

    if tile is None:
        tile = (256, 256)

    if binfile is None:
        binfile = lsmfile
    elif binfile.lower() == 'none':
        binfile = None
    if binfile:
        binfile += '_(z%ic%iy%ix%i)_m%%ip%%it%%03iy%%ix%%i.bin'

    prints('\nOpening LSM file... ', end='', flush=True)
    timer = Timer()

    with TiffFile(lsmfile) as lsm:
        if not lsm.is_lsm:
            prints('\n', lsm, flush=True)
            raise ValueError('not a LSM file')
        series = lsm.series[0]  # first series contains the image
        shape = series.get_shape(False)
        axes = series.get_axes(False)
        dtype = series.dtype
        size = product(shape) * dtype.itemsize

        prints(timer)
        # verbose(lsm, flush=True)
        prints(
            indent(
                'Image',
                f'axes:  {axes}',
                f'shape: {shape}',
                f'dtype: {dtype}',
                f'size:  {size}',
            ),
            flush=True,
        )
        if axes == 'CYX':
            shape = (1, 1) + shape
        elif axes == 'ZCYX':
            shape = (1,) + shape
        elif axes == 'MPCYX':
            shape = shape[:2] + (1, 1) + shape[2:]
        elif axes == 'MPZCYX':
            shape = shape[:2] + (1,) + shape[2:]
        elif not axes.endswith('TZCYX'):
            raise ValueError('not a *TZCYX LSM file')

        prints('Copying image from LSM to BIN files', end='', flush=True)
        timer.start()
        tiles = shape[-2] // tile[-2], shape[-1] // tile[-1]
        if binfile:
            binfile = binfile % (shape[-4], shape[-3], tile[0], tile[1])
        shape = (1,) * (7 - len(shape)) + shape
        # cache for ZCYX stacks and output files
        data = numpy.empty(shape[3:], dtype=dtype)
        out = numpy.empty(
            (shape[-4], shape[-3], tile[0], tile[1]), dtype=dtype
        )
        # iterate over Tiff pages containing data
        pages = iter(series.pages)
        for m in range(shape[0]):  # mosaic axis
            for p in range(shape[1]):  # position axis
                for t in range(shape[2]):  # time axis
                    for z in range(shape[3]):  # z slices
                        page = next(pages)
                        assert page is not None
                        data[z] = page.asarray()
                    for y in range(tiles[0]):  # tile y
                        for x in range(tiles[1]):  # tile x
                            out[:] = data[
                                ...,
                                y * tile[0] : (y + 1) * tile[0],
                                x * tile[1] : (x + 1) * tile[1],
                            ]
                            if binfile:
                                out.tofile(binfile % (m, p, t, y, x))
                            prints('.', end='', flush=True)
        prints(timer, flush=True)


def imshow(
    data: NDArray[Any],
    /,
    *,
    photometric: PHOTOMETRIC | int | str | None = None,
    planarconfig: PLANARCONFIG | int | str | None = None,
    bitspersample: int | None = None,
    nodata: int | float = 0,
    interpolation: str | None = None,
    cmap: Any | None = None,
    vmin: int | float | None = None,
    vmax: int | float | None = None,
    figure: Any = None,
    subplot: Any = None,
    title: str | None = None,
    window_title: str | None = None,
    dpi: int = 96,
    maxdim: int | None = None,
    background: tuple[float, float, float] | str | None = None,
    show: bool = False,
    **kwargs: Any,
) -> tuple[Any, Any, Any]:
    """Plot n-dimensional images with `matplotlib.pyplot`.

    Parameters:
        data:
            Image array to display.
        photometric:
            Color space of image.
        planarconfig:
            How components of each pixel are stored.
        bitspersample:
            Number of bits per channel in integer RGB images.
        interpolation:
            Image interpolation method used in `matplotlib.imshow`.
           The default is 'nearest' for image dimensions > 512,
           else 'bilinear'.
        cmap:
            Colormap mapping non-RGBA scalar data to colors.
            See `matplotlib.colors.Colormap`.
        vmin:
            Minimum of data range covered by colormap.
            By default, the complete range of the data is covered.
        vmax:
            Maximum of data range covered by colormap.
            By default, the complete range of the data is covered.
        figure:
            Matplotlib figure to use for plotting.
            See `matplotlib.figure.Figure`.
        subplot:
            A `matplotlib.pyplot.subplot` axis.
        title:
            Subplot title.
        window_title:
            Window title.
        dpi:
            Resolution of figure.
        maxdim:
            Maximum image width and length.
        background:
            Background color.
        show:
            Display figure.
        **kwargs:
            Additional arguments passed to `matplotlib.pyplot.imshow`.

    Returns:
        Matplotlib figure, subplot, and plot axis.

    """
    # TODO: rewrite detection of isrgb, iscontig
    # TODO: use planarconfig
    if photometric is None:
        photometric = 'RGB'
    if maxdim is None:
        maxdim = 2**16
    isrgb = photometric in {'RGB', 'YCBCR'}  # 'PALETTE', 'YCBCR'

    if data.dtype == 'float16':
        data = data.astype(numpy.float32)

    if data.dtype.kind == 'b':
        isrgb = False

    if isrgb and not (
        data.shape[-1] in {3, 4}
        or (data.ndim > 2 and data.shape[-3] in {3, 4})
    ):
        isrgb = False
        photometric = 'MINISBLACK'

    data = data.squeeze()
    if photometric in {
        None,
        'MINISWHITE',
        'MINISBLACK',
        'CFA',
        'MASK',
        'PALETTE',
        'LOGL',
        'LOGLUV',
        'DEPTH_MAP',
        'SEMANTIC_MASK',
    }:
        data = reshape_nd(data, 2)
    else:
        data = reshape_nd(data, 3)

    dims = data.ndim
    if dims < 2:
        raise ValueError('not an image')
    if dims == 2:
        dims = 0
        isrgb = False
    else:
        if isrgb and data.shape[-3] in {3, 4} and data.shape[-1] not in {3, 4}:
            data = numpy.swapaxes(data, -3, -2)
            data = numpy.swapaxes(data, -2, -1)
        elif not isrgb and (
            data.shape[-1] < data.shape[-2] // 8
            and data.shape[-1] < data.shape[-3] // 8
        ):
            data = numpy.swapaxes(data, -3, -1)
            data = numpy.swapaxes(data, -2, -1)
        isrgb = isrgb and data.shape[-1] in {3, 4}
        dims -= 3 if isrgb else 2

    if interpolation is None:
        threshold = 512
    elif isinstance(interpolation, int):
        threshold = interpolation  # type: ignore[unreachable]
    else:
        threshold = 0

    if isrgb:
        data = data[..., :maxdim, :maxdim, :maxdim]
        if threshold:
            if data.shape[-2] > threshold or data.shape[-3] > threshold:
                interpolation = 'bilinear'
            else:
                interpolation = 'nearest'
    else:
        data = data[..., :maxdim, :maxdim]
        if threshold:
            if data.shape[-1] > threshold or data.shape[-2] > threshold:
                interpolation = 'bilinear'
            else:
                interpolation = 'nearest'

    if photometric == 'PALETTE' and isrgb:
        try:
            datamax = numpy.max(data)
        except ValueError:
            datamax = 1
        if datamax > 255:
            data = data >> 8  # possible precision loss
        data = data.astype('B', copy=False)
    elif data.dtype.kind in 'ui':
        if not (isrgb and data.dtype.itemsize <= 1) or bitspersample is None:
            try:
                bitspersample = int(math.ceil(math.log(data.max(), 2)))
            except Exception:
                bitspersample = data.dtype.itemsize * 8
        elif not isinstance(bitspersample, (int, numpy.integer)):
            # bitspersample can be tuple, such as (5, 6, 5)
            bitspersample = data.dtype.itemsize * 8
        assert bitspersample is not None
        datamax = 2**bitspersample
        if isrgb:
            if bitspersample < 8:
                data = data << (8 - bitspersample)
            elif bitspersample > 8:
                data = data >> (bitspersample - 8)  # precision loss
            data = data.astype('B', copy=False)
    elif data.dtype.kind == 'f':
        if nodata:
            data = data.copy()
            data[data == nodata] = numpy.nan
        try:
            datamax = numpy.nanmax(data)
        except ValueError:
            datamax = 1
        if isrgb and datamax > 1.0:
            if data.dtype.char == 'd':
                data = data.astype('f')
                data /= datamax
            else:
                data = data / datamax
    elif data.dtype.kind == 'b':
        datamax = 1
    elif data.dtype.kind == 'c':
        data = numpy.absolute(data)
        try:
            datamax = numpy.nanmax(data)
        except ValueError:
            datamax = 1

    if isrgb:
        vmin = 0
    else:
        if vmax is None:
            vmax = datamax
        if vmin is None:
            if data.dtype.kind == 'i':
                imin = numpy.iinfo(data.dtype).min
                try:
                    vmin = numpy.min(data)
                except ValueError:
                    vmin = -1
                if vmin == imin:
                    vmin = numpy.min(data[data > imin])
            elif data.dtype.kind == 'f':
                fmin = float(numpy.finfo(data.dtype).min)
                try:
                    vmin = numpy.nanmin(data)
                except ValueError:
                    vmin = 0.0
                if vmin == fmin:
                    vmin = numpy.nanmin(data[data > fmin])
            else:
                vmin = 0

    from matplotlib import pyplot
    from matplotlib.widgets import Slider

    if figure is None:
        pyplot.rc('font', family='sans-serif', weight='normal', size=8)
        figure = pyplot.figure(
            dpi=dpi,
            figsize=(10.3, 6.3),
            frameon=True,
            facecolor='1.0',
            edgecolor='w',
        )
        if window_title is not None:
            try:
                figure.canvas.manager.window.title(window_title)
            except Exception:
                pass
        size = len(title.splitlines()) if title else 1
        pyplot.subplots_adjust(
            bottom=0.03 * (dims + 2),
            top=0.98 - size * 0.03,
            left=0.1,
            right=0.95,
            hspace=0.05,
            wspace=0.0,
        )
    if subplot is None:
        subplot = 111
    subplot = pyplot.subplot(subplot)
    if background is None:
        background = (0.382, 0.382, 0.382)
    subplot.set_facecolor(background)

    if title:
        if isinstance(title, bytes):
            title = title.decode('Windows-1252')
        pyplot.title(title, size=11)

    if cmap is None:
        if data.dtype.char == '?':
            cmap = 'gray'
        elif data.dtype.kind in 'buf' or vmin == 0:
            cmap = 'viridis'
        else:
            cmap = 'coolwarm'
        if photometric == 'MINISWHITE':
            cmap += '_r'

    image = pyplot.imshow(
        numpy.atleast_2d(data[(0,) * dims].squeeze()),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        interpolation=interpolation,
        **kwargs,
    )

    if not isrgb:
        pyplot.colorbar()  # panchor=(0.55, 0.5), fraction=0.05

    def format_coord(x: float, y: float, /) -> str:
        # callback function to format coordinate display in toolbar
        x = int(x + 0.5)
        y = int(y + 0.5)
        try:
            if dims:
                return f'{curaxdat[1][y, x]} @ {current} [{y:4}, {x:4}]'
            return f'{data[y, x]} @ [{y:4}, {x:4}]'
        except IndexError:
            return ''

    def none(event: Any) -> str:
        return ''

    subplot.format_coord = format_coord
    image.get_cursor_data = none  # type: ignore[assignment, method-assign]
    image.format_cursor_data = none  # type: ignore[assignment, method-assign]

    if dims:
        current = list((0,) * dims)
        curaxdat = [0, data[tuple(current)].squeeze()]
        sliders = [
            Slider(
                ax=pyplot.axes((0.125, 0.03 * (axis + 1), 0.725, 0.025)),
                label=f'Dimension {axis}',
                valmin=0,
                valmax=data.shape[axis] - 1,
                valinit=0,
                valfmt=f'%.0f [{data.shape[axis]}]',
            )
            for axis in range(dims)
        ]
        for slider in sliders:
            slider.drawon = False

        def set_image(current, sliders=sliders, data=data):
            # change image and redraw canvas
            curaxdat[1] = data[tuple(current)].squeeze()
            image.set_data(curaxdat[1])
            for ctrl, index in zip(sliders, current):
                ctrl.eventson = False
                ctrl.set_val(index)
                ctrl.eventson = True
            figure.canvas.draw()

        def on_changed(index, axis, data=data, current=current):
            # callback function for slider change event
            index = int(round(index))
            curaxdat[0] = axis
            if index == current[axis]:
                return
            if index >= data.shape[axis]:
                index = 0
            elif index < 0:
                index = data.shape[axis] - 1
            current[axis] = index
            set_image(current)

        def on_keypressed(event, data=data, current=current):
            # callback function for key press event
            key = event.key
            axis = curaxdat[0]
            if str(key) in '0123456789':
                on_changed(key, axis)
            elif key == 'right':
                on_changed(current[axis] + 1, axis)
            elif key == 'left':
                on_changed(current[axis] - 1, axis)
            elif key == 'up':
                curaxdat[0] = 0 if axis == len(data.shape) - 1 else axis + 1
            elif key == 'down':
                curaxdat[0] = len(data.shape) - 1 if axis == 0 else axis - 1
            elif key == 'end':
                on_changed(data.shape[axis] - 1, axis)
            elif key == 'home':
                on_changed(0, axis)

        figure.canvas.mpl_connect('key_press_event', on_keypressed)
        for axis, ctrl in enumerate(sliders):
            ctrl.on_changed(
                lambda k, a=axis: on_changed(k, a)  # type: ignore[misc]
            )

    if show:
        pyplot.show()

    return figure, subplot, image


def askopenfilename(**kwargs: Any) -> str:
    """Return file name(s) from Tkinter's file open dialog."""
    from tkinter import Tk, filedialog

    root = Tk()
    root.withdraw()
    root.update()
    filenames = filedialog.askopenfilename(**kwargs)
    root.destroy()
    return filenames


def main() -> int:
    """Tifffile command line usage main function."""
    import optparse  # TODO: use argparse

    logger().setLevel(logging.INFO)

    parser = optparse.OptionParser(
        usage='usage: %prog [options] path',
        description='Display image and metadata in TIFF file.',
        version=f'%prog {__version__}',
        prog='tifffile',
    )
    opt = parser.add_option
    opt(
        '-p',
        '--page',
        dest='page',
        type='int',
        default=-1,
        help='display single page',
    )
    opt(
        '-s',
        '--series',
        dest='series',
        type='int',
        default=-1,
        help='display select series',
    )
    opt(
        '-l',
        '--level',
        dest='level',
        type='int',
        default=-1,
        help='display pyramid level of series',
    )
    opt(
        '--nomultifile',
        dest='nomultifile',
        action='store_true',
        default=False,
        help='do not read OME series from multiple files',
    )
    opt(
        '--maxplots',
        dest='maxplots',
        type='int',
        default=10,
        help='maximum number of plot windows',
    )
    opt(
        '--interpol',
        dest='interpol',
        metavar='INTERPOL',
        default=None,
        help='image interpolation method',
    )
    opt('--dpi', dest='dpi', type='int', default=96, help='plot resolution')
    opt(
        '--vmin',
        dest='vmin',
        type='int',
        default=None,
        help='minimum value for colormapping',
    )
    opt(
        '--vmax',
        dest='vmax',
        type='int',
        default=None,
        help='maximum value for colormapping',
    )
    opt(
        '--cmap',
        dest='cmap',
        type='str',
        default=None,
        help='colormap name used to map data to colors',
    )
    opt(
        '--maxworkers',
        dest='maxworkers',
        type='int',
        default=0,
        help='maximum number of threads',
    )
    opt(
        '--debug',
        dest='debug',
        action='store_true',
        default=False,
        help='raise exception on failures',
    )
    opt('-v', '--detail', dest='detail', type='int', default=2)
    opt('-q', '--quiet', dest='quiet', action='store_true')

    settings, path_list = parser.parse_args()
    path = ' '.join(path_list)

    if not path:
        path = askopenfilename(
            title='Select a TIFF file', filetypes=TIFF.FILEOPEN_FILTER
        )
        if not path:
            parser.error('No file specified')

    if any(i in path for i in '?*'):
        path_list = glob.glob(path)
        if not path_list:
            print('No files match the pattern')
            return 0
        # TODO: handle image sequences
        path = path_list[0]

    if not settings.quiet:
        print('\nReading TIFF header:', end=' ', flush=True)
    timer = Timer()
    try:
        tif = TiffFile(path, _multifile=not settings.nomultifile)
    except Exception as exc:
        if settings.debug:
            raise
        print(f'\n\n{exc.__class__.__name__}: {exc}')
        return 0

    if not settings.quiet:
        print(timer)

    if tif.is_ome:
        settings.norgb = True

    images: list[tuple[Any, Any, Any]] = []
    if settings.maxplots > 0:
        if not settings.quiet:
            print('Reading image data:', end=' ', flush=True)

        def notnone(x: Any, /) -> Any:
            return next(i for i in x if i is not None)

        timer.start()
        try:
            if settings.page >= 0:
                images = [
                    (
                        tif.asarray(
                            key=settings.page, maxworkers=settings.maxworkers
                        ),
                        tif.pages[settings.page],
                        None,
                    )
                ]
            elif settings.series >= 0:
                series = tif.series[settings.series]
                if settings.level >= 0:
                    level = settings.level
                elif series.is_pyramidal and product(series.shape) > 2**32:
                    level = -1
                    for r in series.levels:
                        level += 1
                        if product(r.shape) < 2**32:
                            break
                else:
                    level = 0
                images = [
                    (
                        tif.asarray(
                            series=settings.series,
                            level=level,
                            maxworkers=settings.maxworkers,
                        ),
                        notnone(tif.series[settings.series]._pages),
                        tif.series[settings.series],
                    )
                ]
            else:
                for i, s in enumerate(tif.series[: settings.maxplots]):
                    if settings.level < 0:
                        level = -1
                        for r in s.levels:
                            level += 1
                            if product(r.shape) < 2**31:
                                break
                    else:
                        level = settings.level
                    try:
                        images.append(
                            (
                                tif.asarray(
                                    series=i,
                                    level=level,
                                    maxworkers=settings.maxworkers,
                                ),
                                notnone(s._pages),
                                tif.series[i],
                            )
                        )
                    except Exception as exc:
                        images.append((None, notnone(s.pages), None))
                        if settings.debug:
                            raise
                        print(f'\nSeries {i} raised {exc!r:.128}... ', end='')
        except Exception as exc:
            if settings.debug:
                raise
            print(f'{exc.__class__.__name__}: {exc}')

        if not settings.quiet:
            print(timer)

    if not settings.quiet:
        print('Generating report:', end='   ', flush=True)
        timer.start()
        try:
            width = os.get_terminal_size()[0]
        except Exception:
            width = 80
        info = tif._str(detail=int(settings.detail), width=width - 1)
        print(timer)
        print()
        print(info)
        print()

    if images and settings.maxplots > 0:
        try:
            from matplotlib import pyplot
        except ImportError as exc:
            logger().warning(f'<tifffile.main> raised {exc!r:.128}')
        else:
            for img, page, series in images:
                if img is None:
                    continue
                keyframe = page.keyframe
                vmin, vmax = settings.vmin, settings.vmax
                if keyframe.nodata:
                    try:
                        if img.dtype.kind == 'f':
                            img[img == keyframe.nodata] = numpy.nan
                            vmin = numpy.nanmin(img)
                        else:
                            vmin = numpy.min(img[img > keyframe.nodata])
                    except ValueError:
                        pass
                if tif.is_stk:
                    try:
                        vmin = tif.stk_metadata[
                            'MinScale'  # type: ignore[index]
                        ]
                        vmax = tif.stk_metadata[
                            'MaxScale'  # type: ignore[index]
                        ]
                    except KeyError:
                        pass
                    else:
                        if vmax <= vmin:
                            vmin, vmax = settings.vmin, settings.vmax
                if series:
                    title = f'{tif}\n{page}\n{series}'
                    window_title = f'{tif.filename} series {series.index}'
                else:
                    title = f'{tif}\n{page}'
                    window_title = f'{tif.filename} page {page.index}'
                photometric = 'MINISBLACK'
                if keyframe.photometric != 3:
                    photometric = PHOTOMETRIC(keyframe.photometric).name
                imshow(
                    img,
                    title=title,
                    window_title=window_title,
                    vmin=vmin,
                    vmax=vmax,
                    cmap=settings.cmap,
                    bitspersample=keyframe.bitspersample,
                    nodata=keyframe.nodata,
                    photometric=photometric,
                    interpolation=settings.interpol,
                    dpi=settings.dpi,
                    show=False,
                )
            pyplot.show()

    tif.close()
    return 0


def bytes2str(
    b: bytes, /, encoding: str | None = None, errors: str = 'strict'
) -> str:
    """Return Unicode string from encoded bytes up to first NULL character."""
    if encoding is None or '16' not in encoding:
        i = b.find(b'\x00')
        if i >= 0:
            b = b[:i]
    else:
        # utf-16
        i = b.find(b'\x00\x00')
        if i >= 0:
            b = b[: i + i % 2]

    try:
        return b.decode('utf-8' if encoding is None else encoding, errors)
    except UnicodeDecodeError:
        if encoding is not None:
            raise
        return b.decode('cp1252', errors)


def bytestr(s: str | bytes, /, encoding: str = 'cp1252') -> bytes:
    """Return bytes from Unicode string, else pass through."""
    return s.encode(encoding) if isinstance(s, str) else s


# aliases and deprecated
TiffReader = TiffFile


if __name__ == '__main__':
    sys.exit(main())

# mypy: allow-untyped-defs, allow-untyped-calls
# mypy: disable-error-code="no-any-return, unreachable, redundant-expr"

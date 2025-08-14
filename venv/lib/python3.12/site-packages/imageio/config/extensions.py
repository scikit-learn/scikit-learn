"""
A set of objects representing each file extension recognized by ImageIO. If an
extension is not listed here it is still supported, as long as there exists a
supporting backend.

"""


class FileExtension:
    """File Extension Metadata

    This class holds information about a image file format associated with a
    given extension. This information is used to track plugins that are known to
    be able to handle a particular format. It also contains additional
    information about a format, which is used when creating the supported format
    docs.

    Plugins known to be able to handle this format are ordered by a ``priority``
    list. This list is used to determine the ideal plugin to use when choosing a
    plugin based on file extension.

    Parameters
    ----------
    extension : str
        The name of the extension including the initial dot, e.g. ".png".
    priority : List
        A list of plugin names (entries in config.known_plugins) that can handle
        this format. The position of a plugin expresses a preference, e.g.
        ["plugin1", "plugin2"] indicates that, if available, plugin1 should be
        preferred over plugin2 when handling a request related to this format.
    name : str
        The full name of the format.
    description : str
        A description of the format.
    external_link : str
        A link to further information about the format. Typically, the format's
        specification.
    volume_support : str
        If True, the format/extension supports volumetric image data.

    Examples
    --------
    >>> FileExtension(
            name="Bitmap",
            extension=".bmp",
            priority=["pillow", "BMP-PIL", "BMP-FI", "ITK"],
            external_link="https://en.wikipedia.org/wiki/BMP_file_format",
        )

    """

    def __init__(
        self,
        *,
        extension,
        priority,
        name=None,
        description=None,
        external_link=None,
        volume_support=False,
    ):
        self.extension = extension
        self.priority = priority
        self.name = name
        self.description = description
        self.external_link = external_link
        self.default_priority = priority.copy()
        self.volume_support = volume_support

    def reset(self):
        self.priority = self.default_priority.copy()


extension_list = [
    FileExtension(
        name="Hasselblad raw",
        extension=".3fr",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Sony alpha",
        extension=".arw",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Animated Portable Network Graphics",
        external_link="https://en.wikipedia.org/wiki/APNG",
        extension=".apng",
        priority=["pillow", "pyav"],
    ),
    FileExtension(
        name="Audio Video Interleave",
        extension=".avi",
        priority=["FFMPEG"],
    ),
    FileExtension(
        name="Casio raw format",
        extension=".bay",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".blp",
        priority=["pillow"],
    ),
    FileExtension(
        name="Bitmap",
        extension=".bmp",
        priority=["pillow", "BMP-PIL", "BMP-FI", "ITK", "pyav", "opencv"],
        external_link="https://en.wikipedia.org/wiki/BMP_file_format",
    ),
    FileExtension(
        name="Device-Independent Bitmap",
        extension=".dip",
        priority=["opencv"],
        external_link="https://en.wikipedia.org/wiki/BMP_file_format",
    ),
    FileExtension(
        name="Re-Volt mipmap",
        extension=".bmq",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Binary Structured Data Format",
        extension=".bsdf",
        priority=["BSDF"],
        external_link="http://bsdf.io/",
    ),
    FileExtension(
        name="Binary Universal Form for the Representation of meteorological data",
        extension=".bufr",
        priority=["pillow", "BUFR-PIL"],
    ),
    FileExtension(
        name="Silicon Graphics Image",
        extension=".bw",
        priority=["pillow", "SGI-PIL", "SGI-FI"],
    ),
    FileExtension(
        name="Scirra Construct",
        extension=".cap",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="AMETEK High Speed Camera Format",
        extension=".cine",
        priority=["RAW-FI"],
        external_link="https://phantomhighspeed-knowledge.secure.force.com/servlet/fileField?id=0BE1N000000kD2i#:~:text=Cine%20is%20a%20video%20file,camera%20model%20and%20image%20resolution",
    ),
    FileExtension(extension=".cr2", priority=["RAW-FI"]),
    FileExtension(
        extension=".crw",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".cs1",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Computerized Tomography",
        extension=".ct",
        priority=["DICOM"],
    ),
    FileExtension(
        name="Windows Cursor Icons",
        extension=".cur",
        priority=["pillow", "CUR-PIL"],
    ),
    FileExtension(
        name="Dr. Halo",
        extension=".cut",
        priority=["CUT-FI"],
    ),
    FileExtension(
        extension=".dc2",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="DICOM file format",
        extension=".dcm",
        priority=["DICOM", "ITK"],
    ),
    FileExtension(
        extension=".dcr",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Intel DCX",
        extension=".dcx",
        priority=["pillow", "DCX-PIL"],
    ),
    FileExtension(
        name="DirectX Texture Container",
        extension=".dds",
        priority=["pillow", "DDS-FI", "DDS-PIL"],
    ),
    FileExtension(
        name="Windows Bitmap",
        extension=".dib",
        priority=["pillow", "DIB-PIL"],
    ),
    FileExtension(
        name="DICOM file format",
        extension=".dicom",
        priority=["ITK"],
    ),
    FileExtension(
        extension=".dng",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".drf",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".dsc",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Enhanced Compression Wavelet",
        extension=".ecw",
        priority=["GDAL"],
    ),
    FileExtension(
        name="Windows Metafile",
        extension=".emf",
        priority=["pillow", "WMF-PIL"],
    ),
    FileExtension(
        name="Encapsulated Postscript",
        extension=".eps",
        priority=["pillow", "EPS-PIL"],
    ),
    FileExtension(
        extension=".erf",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="OpenEXR",
        extension=".exr",
        external_link="https://openexr.readthedocs.io/en/latest/",
        priority=["EXR-FI", "pyav", "opencv"],
    ),
    FileExtension(
        extension=".fff",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Flexible Image Transport System File",
        extension=".fit",
        priority=["pillow", "FITS-PIL", "FITS"],
    ),
    FileExtension(
        name="Flexible Image Transport System File",
        extension=".fits",
        priority=["pillow", "FITS-PIL", "FITS", "pyav"],
    ),
    FileExtension(
        name="Autodesk FLC Animation",
        extension=".flc",
        priority=["pillow", "FLI-PIL"],
    ),
    FileExtension(
        name="Autodesk FLI Animation",
        extension=".fli",
        priority=["pillow", "FLI-PIL"],
    ),
    FileExtension(
        name="Kodak FlashPix",
        extension=".fpx",
        priority=["pillow", "FPX-PIL"],
    ),
    FileExtension(
        name="Independence War 2: Edge Of Chaos Texture Format",
        extension=".ftc",
        priority=["pillow", "FTEX-PIL"],
    ),
    FileExtension(
        name="Flexible Image Transport System File",
        extension=".fts",
        priority=["FITS"],
    ),
    FileExtension(
        name="Independence War 2: Edge Of Chaos Texture Format",
        extension=".ftu",
        priority=["pillow", "FTEX-PIL"],
    ),
    FileExtension(
        name="Flexible Image Transport System File",
        extension=".fz",
        priority=["FITS"],
    ),
    FileExtension(
        name="Raw fax format CCITT G.3",
        extension=".g3",
        priority=["G3-FI"],
    ),
    FileExtension(
        name="GIMP brush file",
        extension=".gbr",
        priority=["pillow", "GBR-PIL"],
    ),
    FileExtension(
        name="Grassroots DICOM",
        extension=".gdcm",
        priority=["ITK"],
    ),
    FileExtension(
        name="Graphics Interchange Format",
        extension=".gif",
        priority=["pillow", "GIF-PIL", "pyav"],
    ),
    FileExtension(
        name="UMDS GIPL",
        extension=".gipl",
        priority=["ITK"],
    ),
    FileExtension(
        name="gridded meteorological data",
        extension=".grib",
        priority=["pillow", "GRIB-PIL"],
    ),
    FileExtension(
        name="Hierarchical Data Format 5",
        extension=".h5",
        priority=["pillow", "HDF5-PIL"],
    ),
    FileExtension(
        name="Hierarchical Data Format 5",
        extension=".hdf",
        priority=["pillow", "HDF5-PIL"],
    ),
    FileExtension(
        name="Hierarchical Data Format 5",
        extension=".hdf5",
        priority=["ITK"],
    ),
    FileExtension(
        name="JPEG Extended Range",
        extension=".hdp",
        priority=["JPEG-XR-FI"],
    ),
    FileExtension(
        name="High Dynamic Range Image",
        extension=".hdr",
        priority=["HDR-FI", "ITK", "opencv"],
    ),
    FileExtension(
        extension=".ia",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".icb",
        priority=["pillow"],
    ),
    FileExtension(
        name="Mac OS Icon File",
        extension=".icns",
        priority=["pillow", "ICNS-PIL"],
    ),
    FileExtension(
        name="Windows Icon File",
        extension=".ico",
        priority=["pillow", "ICO-FI", "ICO-PIL", "pyav"],
    ),
    FileExtension(
        name="ILBM Interleaved Bitmap",
        extension=".iff",
        priority=["IFF-FI"],
    ),
    FileExtension(
        name="IPTC/NAA",
        extension=".iim",
        priority=["pillow", "IPTC-PIL"],
    ),
    FileExtension(
        extension=".iiq",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="IFUNC Image Memory",
        extension=".im",
        priority=["pillow", "IM-PIL"],
    ),
    FileExtension(
        extension=".img",
        priority=["ITK", "GDAL"],
    ),
    FileExtension(
        extension=".img.gz",
        priority=["ITK"],
    ),
    FileExtension(
        name="IM Tools",
        extension=".IMT",
        priority=["pillow", "IMT-PIL"],
    ),
    FileExtension(
        name="Image Processing Lab",
        extension=".ipl",
        priority=["ITK"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".j2c",
        priority=["pillow", "J2K-FI", "JPEG2000-PIL", "pyav"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".j2k",
        priority=["pillow", "J2K-FI", "JPEG2000-PIL", "pyav"],
    ),
    FileExtension(
        name="JPEG",
        extension=".jfif",
        priority=["pillow", "JPEG-PIL"],
    ),
    FileExtension(
        name="JPEG",
        extension=".jif",
        priority=["JPEG-FI"],
    ),
    FileExtension(
        name="JPEG Network Graphics",
        extension=".jng",
        priority=["JNG-FI"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".jp2",
        priority=["pillow", "JP2-FI", "JPEG2000-PIL", "pyav", "opencv"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".jpc",
        priority=["pillow", "JPEG2000-PIL"],
    ),
    FileExtension(
        name="JPEG",
        extension=".jpe",
        priority=["pillow", "JPEG-FI", "JPEG-PIL", "opencv"],
    ),
    FileExtension(
        name="Joint Photographic Experts Group",
        extension=".jpeg",
        priority=["pillow", "JPEG-PIL", "JPEG-FI", "ITK", "GDAL", "pyav", "opencv"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".jpf",
        priority=["pillow", "JPEG2000-PIL"],
    ),
    FileExtension(
        name="Joint Photographic Experts Group",
        extension=".jpg",
        priority=["pillow", "JPEG-PIL", "JPEG-FI", "ITK", "GDAL", "pyav", "opencv"],
    ),
    FileExtension(
        name="JPEG 2000",
        extension=".jpx",
        priority=["pillow", "JPEG2000-PIL"],
    ),
    FileExtension(
        name="JPEG Extended Range",
        extension=".jxr",
        priority=["JPEG-XR-FI"],
    ),
    FileExtension(
        extension=".k25",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".kc2",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".kdc",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="C64 Koala Graphics",
        extension=".koa",
        priority=["KOALA-FI"],
    ),
    FileExtension(
        name="ILBM Interleaved Bitmap",
        extension=".lbm",
        priority=["IFF-FI"],
    ),
    FileExtension(
        name="Lytro F01",
        extension=".lfp",
        priority=["LYTRO-LFP"],
    ),
    FileExtension(
        name="Lytro Illum",
        extension=".lfr",
        priority=["LYTRO-LFR"],
    ),
    FileExtension(
        name="ZEISS LSM",
        extension=".lsm",
        priority=["tifffile", "ITK", "TIFF"],
    ),
    FileExtension(
        name="McIdas area file",
        extension=".MCIDAS",
        priority=["pillow", "MCIDAS-PIL"],
        external_link="https://www.ssec.wisc.edu/mcidas/doc/prog_man/2003print/progman2003-formats.html",
    ),
    FileExtension(
        extension=".mdc",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".mef",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="FreeSurfer File Format",
        extension=".mgh",
        priority=["ITK"],
    ),
    FileExtension(
        name="ITK MetaImage",
        extension=".mha",
        priority=["ITK"],
    ),
    FileExtension(
        name="ITK MetaImage Header",
        extension=".mhd",
        priority=["ITK"],
    ),
    FileExtension(
        name="Microsoft Image Composer",
        extension=".mic",
        priority=["pillow", "MIC-PIL"],
    ),
    FileExtension(
        name="Matroska Multimedia Container",
        extension=".mkv",
        priority=["FFMPEG", "pyav"],
    ),
    FileExtension(
        name="Medical Imaging NetCDF",
        extension=".mnc",
        priority=["ITK"],
    ),
    FileExtension(
        name="Medical Imaging NetCDF 2",
        extension=".mnc2",
        priority=["ITK"],
    ),
    FileExtension(
        name="Leaf Raw Image Format",
        extension=".mos",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="QuickTime File Format",
        extension=".mov",
        priority=["FFMPEG", "pyav"],
    ),
    FileExtension(
        name="MPEG-4 Part 14",
        extension=".mp4",
        priority=["FFMPEG", "pyav"],
    ),
    FileExtension(
        name="MPEG-1 Moving Picture Experts Group",
        extension=".mpeg",
        priority=["FFMPEG", "pyav"],
    ),
    FileExtension(
        name="Moving Picture Experts Group",
        extension=".mpg",
        priority=["pillow", "FFMPEG", "pyav"],
    ),
    FileExtension(
        name="JPEG Multi-Picture Format",
        extension=".mpo",
        priority=["pillow", "MPO-PIL"],
    ),
    FileExtension(
        name="Magnetic resonance imaging",
        extension=".mri",
        priority=["DICOM"],
    ),
    FileExtension(
        extension=".mrw",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Windows Paint",
        extension=".msp",
        priority=["pillow", "MSP-PIL"],
    ),
    FileExtension(
        extension=".nef",
        priority=["RAW-FI", "rawpy"],
    ),
    FileExtension(
        extension=".nhdr",
        priority=["ITK"],
    ),
    FileExtension(
        extension=".nia",
        priority=["ITK"],
    ),
    FileExtension(
        extension=".nii",
        priority=["ITK"],
    ),
    FileExtension(
        name="nii.gz",
        extension=".nii.gz",
        priority=["ITK"],
    ),
    FileExtension(
        name="Numpy Array",
        extension=".npz",
        priority=["NPZ"],
        volume_support=True,
    ),
    FileExtension(
        extension=".nrrd",
        priority=["ITK"],
    ),
    FileExtension(
        extension=".nrw",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".orf",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".palm",
        priority=["pillow"],
    ),
    FileExtension(
        name="Portable Bitmap",
        extension=".pbm",
        priority=["PGM-FI", "PGMRAW-FI", "pyav", "opencv"],
    ),
    FileExtension(
        name="Kodak PhotoCD",
        extension=".pcd",
        priority=["pillow", "PCD-FI", "PCD-PIL"],
    ),
    FileExtension(
        name="Macintosh PICT",
        extension=".pct",
        priority=["PICT-FI"],
    ),
    FileExtension(
        name="Zsoft Paintbrush",
        extension=".PCX",
        priority=["pillow", "PCX-FI", "PCX-PIL"],
    ),
    FileExtension(
        extension=".pdf",
        priority=["pillow"],
    ),
    FileExtension(
        extension=".pef",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".pfm",
        priority=["PFM-FI", "pyav", "opencv"],
    ),
    FileExtension(
        name="Portable Greymap",
        extension=".pgm",
        priority=["pillow", "PGM-FI", "PGMRAW-FI", "pyav", "opencv"],
    ),
    FileExtension(
        name="Macintosh PICT",
        extension=".pic",
        priority=["PICT-FI", "ITK", "opencv"],
    ),
    FileExtension(
        name="Macintosh PICT",
        extension=".pict",
        priority=["PICT-FI"],
    ),
    FileExtension(
        name="Portable Network Graphics",
        extension=".png",
        priority=["pillow", "PNG-PIL", "PNG-FI", "ITK", "pyav", "opencv"],
    ),
    FileExtension(
        name="Portable Image Format",
        extension=".pnm",
        priority=["pillow", "opencv"],
    ),
    FileExtension(
        name="Pbmplus image",
        extension=".ppm",
        priority=["pillow", "PPM-PIL", "pyav"],
    ),
    FileExtension(
        name="Pbmplus image",
        extension=".pbm",
        priority=["pillow", "PPM-PIL", "PPM-FI"],
    ),
    FileExtension(
        name="Portable image format",
        extension=".pxm",
        priority=["opencv"],
    ),
    FileExtension(
        name="Portable Pixelmap (ASCII)",
        extension=".ppm",
        priority=["PPM-FI", "opencv"],
    ),
    FileExtension(
        name="Portable Pixelmap (Raw)",
        extension=".ppm",
        priority=["PPMRAW-FI"],
    ),
    FileExtension(
        name="Ghostscript",
        extension=".ps",
        priority=["pillow", "EPS-PIL"],
    ),
    FileExtension(
        name="Adope Photoshop 2.5 and 3.0",
        extension=".psd",
        priority=["pillow", "PSD-PIL", "PSD-FI"],
    ),
    FileExtension(
        extension=".ptx",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".pxn",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="PIXAR raster image",
        extension=".pxr",
        priority=["pillow", "PIXAR-PIL"],
    ),
    FileExtension(
        extension=".qtk",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".raf",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Sun Raster File",
        extension=".ras",
        priority=["pillow", "SUN-PIL", "RAS-FI", "pyav", "opencv"],
    ),
    FileExtension(
        name="Sun Raster File",
        extension=".sr",
        priority=["opencv"],
    ),
    FileExtension(
        extension=".raw",
        priority=["RAW-FI", "LYTRO-ILLUM-RAW", "LYTRO-F01-RAW", "rawpy"],
    ),
    FileExtension(
        extension=".rdc",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Silicon Graphics Image",
        extension=".rgb",
        priority=["pillow", "SGI-PIL"],
    ),
    FileExtension(
        name="Silicon Graphics Image",
        extension=".rgba",
        priority=["pillow", "SGI-PIL"],
    ),
    FileExtension(
        extension=".rw2",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".rwl",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".rwz",
        priority=["RAW-FI"],
    ),
    FileExtension(
        name="Silicon Graphics Image",
        extension=".sgi",
        priority=["pillow", "SGI-PIL", "pyav"],
    ),
    FileExtension(
        name="SPE File Format",
        extension=".spe",
        priority=["SPE"],
    ),
    FileExtension(
        extension=".SPIDER",
        priority=["pillow", "SPIDER-PIL"],
    ),
    FileExtension(
        extension=".sr2",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".srf",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".srw",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".sti",
        priority=["RAW-FI"],
    ),
    FileExtension(
        extension=".stk",
        priority=["tifffile", "TIFF"],
    ),
    FileExtension(
        name="ShockWave Flash",
        extension=".swf",
        priority=["SWF", "pyav"],
    ),
    FileExtension(
        name="Truevision TGA",
        extension=".targa",
        priority=["pillow", "TARGA-FI"],
    ),
    FileExtension(
        name="Truevision TGA",
        extension=".tga",
        priority=["pillow", "TGA-PIL", "TARGA-FI", "pyav"],
    ),
    FileExtension(
        name="Tagged Image File",
        extension=".tif",
        priority=[
            "tifffile",
            "TIFF",
            "pillow",
            "TIFF-PIL",
            "TIFF-FI",
            "FEI",
            "ITK",
            "GDAL",
            "pyav",
            "opencv",
        ],
        volume_support=True,
    ),
    FileExtension(
        name="Tagged Image File Format",
        extension=".tiff",
        priority=[
            "tifffile",
            "TIFF",
            "pillow",
            "TIFF-PIL",
            "TIFF-FI",
            "FEI",
            "ITK",
            "GDAL",
            "pyav",
            "opencv",
        ],
        volume_support=True,
    ),
    FileExtension(
        extension=".vda",
        priority=["pillow"],
    ),
    FileExtension(
        extension=".vst",
        priority=["pillow"],
    ),
    FileExtension(
        extension=".vtk",
        priority=["ITK"],
    ),
    FileExtension(
        name="Wireless Bitmap",
        extension=".wap",
        priority=["WBMP-FI"],
    ),
    FileExtension(
        name="Wireless Bitmap",
        extension=".wbm",
        priority=["WBMP-FI"],
    ),
    FileExtension(
        name="Wireless Bitmap",
        extension=".wbmp",
        priority=["WBMP-FI"],
    ),
    FileExtension(
        name="JPEG Extended Range",
        extension=".wdp",
        priority=["JPEG-XR-FI"],
    ),
    FileExtension(
        name="Matroska",
        extension=".webm",
        priority=["FFMPEG", "pyav"],
    ),
    FileExtension(
        name="Google WebP",
        extension=".webp",
        priority=["pillow", "WEBP-FI", "pyav", "opencv"],
    ),
    FileExtension(
        name="Windows Meta File",
        extension=".wmf",
        priority=["pillow", "WMF-PIL"],
    ),
    FileExtension(
        name="Windows Media Video",
        extension=".wmv",
        priority=["FFMPEG"],
    ),
    FileExtension(
        name="X11 Bitmap",
        extension=".xbm",
        priority=["pillow", "XBM-PIL", "XBM-FI", "pyav"],
    ),
    FileExtension(
        name="X11 Pixel Map",
        extension=".xpm",
        priority=["pillow", "XPM-PIL", "XPM-FI"],
    ),
    FileExtension(
        name="Thumbnail Image",
        extension=".XVTHUMB",
        priority=["pillow", "XVTHUMB-PIL"],
    ),
    FileExtension(
        extension=".dpx",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".im1",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".im24",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".im8",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".jls",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".ljpg",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".pam",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".pcx",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".pgmyuv",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".pix",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".ppm",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".rs",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".sun",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".sunras",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".xface",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".xwd",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".y",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".3g2",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".3gp",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".f4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".ism",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".isma",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".ismv",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".m4a",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".m4b",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".mj2",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP (3GPP file format)",
        extension=".psp",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".3g2",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".3gp",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".f4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".ism",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".isma",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".ismv",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".m4a",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".m4b",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".mj2",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GP2 (3GPP2 file format)",
        extension=".psp",
        priority=["pyav"],
    ),
    FileExtension(
        name="3GPP AMR",
        extension=".amr",
        priority=["pyav"],
    ),
    FileExtension(
        name="a64 - video for Commodore 64",
        extension=".A64",
        priority=["pyav"],
    ),
    FileExtension(
        name="a64 - video for Commodore 64",
        extension=".a64",
        priority=["pyav"],
    ),
    FileExtension(
        name="Adobe Filmstrip",
        extension=".flm",
        priority=["pyav"],
    ),
    FileExtension(
        name="AMV",
        extension=".amv",
        priority=["pyav"],
    ),
    FileExtension(
        name="ASF (Advanced / Active Streaming Format)",
        extension=".asf",
        priority=["pyav"],
    ),
    FileExtension(
        name="ASF (Advanced / Active Streaming Format)",
        extension=".asf",
        priority=["pyav"],
    ),
    FileExtension(
        name="ASF (Advanced / Active Streaming Format)",
        extension=".wmv",
        priority=["pyav"],
    ),
    FileExtension(
        name="ASF (Advanced / Active Streaming Format)",
        extension=".wmv",
        priority=["pyav"],
    ),
    FileExtension(
        name="AV1 Annex B",
        extension=".obu",
        priority=["pyav"],
    ),
    FileExtension(
        name="AV1 low overhead OBU",
        extension=".obu",
        priority=["pyav"],
    ),
    FileExtension(
        name="AVI (Audio Video Interleaved)",
        extension=".avi",
        priority=["pyav"],
    ),
    FileExtension(
        name="AVR (Audio Visual Research)",
        extension=".avr",
        priority=["pyav"],
    ),
    FileExtension(
        name="Beam Software SIFF",
        extension=".vb",
        priority=["pyav"],
    ),
    FileExtension(
        name="CD Graphics",
        extension=".cdg",
        priority=["pyav"],
    ),
    FileExtension(
        name="Commodore CDXL video",
        extension=".cdxl",
        priority=["pyav"],
    ),
    FileExtension(
        name="Commodore CDXL video",
        extension=".xl",
        priority=["pyav"],
    ),
    FileExtension(
        name="DASH Muxer",
        extension=".mpd",
        priority=["pyav"],
    ),
    FileExtension(
        name="Digital Pictures SGA",
        extension=".sga",
        priority=["pyav"],
    ),
    FileExtension(
        name="Discworld II BMV",
        extension=".bmv",
        priority=["pyav"],
    ),
    FileExtension(
        name="DV (Digital Video)",
        extension=".dif",
        priority=["pyav"],
    ),
    FileExtension(
        name="DV (Digital Video)",
        extension=".dv",
        priority=["pyav"],
    ),
    FileExtension(
        name="F4V Adobe Flash Video",
        extension=".f4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="FLV (Flash Video)",
        extension=".flv",
        priority=["pyav"],
    ),
    FileExtension(
        name="GXF (General eXchange Format)",
        extension=".gxf",
        priority=["pyav"],
    ),
    FileExtension(
        name="iCE Draw File",
        extension=".idf",
        priority=["pyav"],
    ),
    FileExtension(
        name="IFV CCTV DVR",
        extension=".ifv",
        priority=["pyav"],
    ),
    FileExtension(
        name="iPod H.264 MP4 (MPEG-4 Part 14)",
        extension=".m4a",
        priority=["pyav"],
    ),
    FileExtension(
        name="iPod H.264 MP4 (MPEG-4 Part 14)",
        extension=".m4b",
        priority=["pyav"],
    ),
    FileExtension(
        name="iPod H.264 MP4 (MPEG-4 Part 14)",
        extension=".m4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="IVR (Internet Video Recording)",
        extension=".ivr",
        priority=["pyav"],
    ),
    FileExtension(
        name="Konami PS2 SVAG",
        extension=".svag",
        priority=["pyav"],
    ),
    FileExtension(
        name="KUX (YouKu)",
        extension=".kux",
        priority=["pyav"],
    ),
    FileExtension(
        name="live RTMP FLV (Flash Video)",
        extension=".flv",
        priority=["pyav"],
    ),
    FileExtension(
        name="Loki SDL MJPEG",
        extension=".mjpg",
        priority=["pyav"],
    ),
    FileExtension(
        name="LVF",
        extension=".lvf",
        priority=["pyav"],
    ),
    FileExtension(
        name="Matroska / WebM",
        extension=".mk3d",
        priority=["pyav"],
    ),
    FileExtension(
        name="Matroska / WebM",
        extension=".mka",
        priority=["pyav"],
    ),
    FileExtension(
        name="Matroska / WebM",
        extension=".mks",
        priority=["pyav"],
    ),
    FileExtension(
        name="Microsoft XMV",
        extension=".xmv",
        priority=["pyav"],
    ),
    FileExtension(
        name="MIME multipart JPEG",
        extension=".mjpg",
        priority=["pyav"],
    ),
    FileExtension(
        name="MobiClip MODS",
        extension=".mods",
        priority=["pyav"],
    ),
    FileExtension(
        name="MobiClip MOFLEX",
        extension=".moflex",
        priority=["pyav"],
    ),
    FileExtension(
        name="Motion Pixels MVI",
        extension=".mvi",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".3g2",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".3gp",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".f4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".ism",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".isma",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".ismv",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".m4a",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".m4b",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".mj2",
        priority=["pyav"],
    ),
    FileExtension(
        name="MP4 (MPEG-4 Part 14)",
        extension=".psp",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-2 PS (DVD VOB)",
        extension=".dvd",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-2 PS (SVCD)",
        extension=".vob",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-2 PS (VOB)",
        extension=".vob",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-TS (MPEG-2 Transport Stream)",
        extension=".m2t",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-TS (MPEG-2 Transport Stream)",
        extension=".m2ts",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-TS (MPEG-2 Transport Stream)",
        extension=".mts",
        priority=["pyav"],
    ),
    FileExtension(
        name="MPEG-TS (MPEG-2 Transport Stream)",
        extension=".ts",
        priority=["pyav"],
    ),
    FileExtension(
        name="Musepack",
        extension=".mpc",
        priority=["pyav"],
    ),
    FileExtension(
        name="MXF (Material eXchange Format) Operational Pattern Atom",
        extension=".mxf",
        priority=["pyav"],
    ),
    FileExtension(
        name="MXF (Material eXchange Format)",
        extension=".mxf",
        priority=["pyav"],
    ),
    FileExtension(
        name="MxPEG clip",
        extension=".mxg",
        priority=["pyav"],
    ),
    FileExtension(
        name="NC camera feed",
        extension=".v",
        priority=["pyav"],
    ),
    FileExtension(
        name="NUT",
        extension=".nut",
        priority=["pyav"],
    ),
    FileExtension(
        name="Ogg Video",
        extension=".ogv",
        priority=["pyav"],
    ),
    FileExtension(
        name="Ogg",
        extension=".ogg",
        priority=["pyav"],
    ),
    FileExtension(
        name="On2 IVF",
        extension=".ivf",
        priority=["pyav"],
    ),
    FileExtension(
        name="PSP MP4 (MPEG-4 Part 14)",
        extension=".psp",
        priority=["pyav"],
    ),
    FileExtension(
        name="Psygnosis YOP",
        extension=".yop",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".3g2",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".3gp",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".f4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".ism",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".isma",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".ismv",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".m4a",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".m4b",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".mj2",
        priority=["pyav"],
    ),
    FileExtension(
        name="QuickTime / MOV",
        extension=".psp",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw AVS2-P2/IEEE1857.4 video",
        extension=".avs",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw AVS2-P2/IEEE1857.4 video",
        extension=".avs2",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw AVS3-P2/IEEE1857.10",
        extension=".avs3",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw Chinese AVS (Audio Video Standard) video",
        extension=".cavs",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw Dirac",
        extension=".drc",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw Dirac",
        extension=".vc2",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw DNxHD (SMPTE VC-3)",
        extension=".dnxhd",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw DNxHD (SMPTE VC-3)",
        extension=".dnxhr",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw GSM",
        extension=".gsm",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw H.261",
        extension=".h261",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw H.263",
        extension=".h263",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw H.264 video",
        extension=".264",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw H.264 video",
        extension=".avc",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw H.264 video",
        extension=".h264",
        priority=["pyav", "FFMPEG"],
    ),
    FileExtension(
        name="raw H.264 video",
        extension=".h26l",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw HEVC video",
        extension=".265",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw HEVC video",
        extension=".h265",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw HEVC video",
        extension=".hevc",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw id RoQ",
        extension=".roq",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw Ingenient MJPEG",
        extension=".cgi",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw IPU Video",
        extension=".ipu",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MJPEG 2000 video",
        extension=".j2k",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MJPEG video",
        extension=".mjpeg",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MJPEG video",
        extension=".mjpg",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MJPEG video",
        extension=".mpo",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MPEG-1 video",
        extension=".m1v",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MPEG-1 video",
        extension=".mpeg",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MPEG-1 video",
        extension=".mpg",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MPEG-2 video",
        extension=".m2v",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw MPEG-4 video",
        extension=".m4v",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw VC-1 video",
        extension=".vc1",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw video",
        extension=".cif",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw video",
        extension=".qcif",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw video",
        extension=".rgb",
        priority=["pyav"],
    ),
    FileExtension(
        name="raw video",
        extension=".yuv",
        priority=["pyav"],
    ),
    FileExtension(
        name="RealMedia",
        extension=".rm",
        priority=["pyav"],
    ),
    FileExtension(
        name="SDR2",
        extension=".sdr2",
        priority=["pyav"],
    ),
    FileExtension(
        name="Sega FILM / CPK",
        extension=".cpk",
        priority=["pyav"],
    ),
    FileExtension(
        name="SER (Simple uncompressed video format for astronomical capturing)",
        extension=".ser",
        priority=["pyav"],
    ),
    FileExtension(
        name="Simbiosis Interactive IMX",
        extension=".imx",
        priority=["pyav"],
    ),
    FileExtension(
        name="Square SVS",
        extension=".svs",
        priority=["tifffile", "pyav"],
    ),
    FileExtension(
        name="TiVo TY Stream",
        extension=".ty",
        priority=["pyav"],
    ),
    FileExtension(
        name="TiVo TY Stream",
        extension=".ty+",
        priority=["pyav"],
    ),
    FileExtension(
        name="Uncompressed 4:2:2 10-bit",
        extension=".v210",
        priority=["pyav"],
    ),
    FileExtension(
        name="Uncompressed 4:2:2 10-bit",
        extension=".yuv10",
        priority=["pyav"],
    ),
    FileExtension(
        name="VC-1 test bitstream",
        extension=".rcv",
        priority=["pyav"],
    ),
    FileExtension(
        name="Video CCTV DAT",
        extension=".dat",
        priority=["pyav"],
    ),
    FileExtension(
        name="Video DAV",
        extension=".dav",
        priority=["pyav"],
    ),
    FileExtension(
        name="Vivo",
        extension=".viv",
        priority=["pyav"],
    ),
    FileExtension(
        name="WebM Chunk Muxer",
        extension=".chk",
        priority=["pyav"],
    ),
    FileExtension(
        name="WebM",
        extension=".mk3d",
        priority=["pyav"],
    ),
    FileExtension(
        name="WebM",
        extension=".mka",
        priority=["pyav"],
    ),
    FileExtension(
        name="WebM",
        extension=".mks",
        priority=["pyav"],
    ),
    FileExtension(
        name="Windows Television (WTV)",
        extension=".wtv",
        priority=["pyav"],
    ),
    FileExtension(
        name="Xilam DERF",
        extension=".adp",
        priority=["pyav"],
    ),
    FileExtension(
        name="YUV4MPEG pipe",
        extension=".y4m",
        priority=["pyav"],
    ),
    FileExtension(
        extension=".qpi",
        priority=["tifffile"],
    ),
    FileExtension(
        name="PCO Camera",
        extension=".pcoraw",
        priority=["tifffile"],
    ),
    FileExtension(
        name="PCO Camera",
        extension=".rec",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Perkin Elmer Vectra",
        extension=".qptiff",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Pyramid Encoded TIFF",
        extension=".ptiff",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Pyramid Encoded TIFF",
        extension=".ptif",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Opticks Gel",
        extension=".gel",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Zoomify Image Format",
        extension=".zif",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Hamamatsu Slide Scanner",
        extension=".ndpi",
        priority=["tifffile"],
    ),
    FileExtension(
        name="Roche Digital Pathology",
        extension=".bif",
        priority=["tifffile"],
    ),
    FileExtension(
        extension=".tf8",
        priority=["tifffile"],
    ),
    FileExtension(
        extension=".btf",
        priority=["tifffile"],
    ),
    FileExtension(
        name="High Efficiency Image File Format",
        extension=".heic",
        priority=["pillow"],
    ),
    FileExtension(
        name="AV1 Image File Format",
        extension=".avif",
        priority=["pillow"],
    ),
]
extension_list.sort(key=lambda x: x.extension)


known_extensions = dict()
for ext in extension_list:
    if ext.extension not in known_extensions:
        known_extensions[ext.extension] = list()
    known_extensions[ext.extension].append(ext)

extension_list = [ext for ext_list in known_extensions.values() for ext in ext_list]

_video_extension_strings = [
    ".264",
    ".265",
    ".3g2",
    ".3gp",
    ".a64",
    ".A64",
    ".adp",
    ".amr",
    ".amv",
    ".asf",
    ".avc",
    ".avi",
    ".avr",
    ".avs",
    ".avs2",
    ".avs3",
    ".bmv",
    ".cavs",
    ".cdg",
    ".cdxl",
    ".cgi",
    ".chk",
    ".cif",
    ".cpk",
    ".dat",
    ".dav",
    ".dif",
    ".dnxhd",
    ".dnxhr",
    ".drc",
    ".dv",
    ".dvd",
    ".f4v",
    ".flm",
    ".flv",
    ".gsm",
    ".gxf",
    ".h261",
    ".h263",
    ".h264",
    ".h265",
    ".h26l",
    ".hevc",
    ".idf",
    ".ifv",
    ".imx",
    ".ipu",
    ".ism",
    ".isma",
    ".ismv",
    ".ivf",
    ".ivr",
    ".j2k",
    ".kux",
    ".lvf",
    ".m1v",
    ".m2t",
    ".m2ts",
    ".m2v",
    ".m4a",
    ".m4b",
    ".m4v",
    ".mj2",
    ".mjpeg",
    ".mjpg",
    ".mk3d",
    ".mka",
    ".mks",
    ".mkv",
    ".mods",
    ".moflex",
    ".mov",
    ".mp4",
    ".mpc",
    ".mpd",
    ".mpeg",
    ".mpg",
    ".mpo",
    ".mts",
    ".mvi",
    ".mxf",
    ".mxg",
    ".nut",
    ".obu",
    ".ogg",
    ".ogv",
    ".psp",
    ".qcif",
    ".rcv",
    ".rgb",
    ".rm",
    ".roq",
    ".sdr2",
    ".ser",
    ".sga",
    ".svag",
    ".svs",
    ".ts",
    ".ty",
    ".ty+",
    ".v",
    ".v210",
    ".vb",
    ".vc1",
    ".vc2",
    ".viv",
    ".vob",
    ".webm",
    ".wmv",
    ".wtv",
    ".xl",
    ".xmv",
    ".y4m",
    ".yop",
    ".yuv",
    ".yuv10",
]
video_extensions = list()
for ext_string in _video_extension_strings:
    formats = known_extensions[ext_string]
    video_extensions.append(formats[0])
video_extensions.sort(key=lambda x: x.extension)

"""skimage viewer"""


def main():
    import skimage.io as io
    import sys

    if len(sys.argv) != 2:
        print("Usage: skivi <image-file>")
        sys.exit(-1)

    io.use_plugin('qt')
    io.imshow(io.imread(sys.argv[1]), fancy=True)
    io.show()

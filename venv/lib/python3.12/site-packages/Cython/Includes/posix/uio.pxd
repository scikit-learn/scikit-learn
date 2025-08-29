# https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/sys_uio.h.html

from posix.types cimport off_t


cdef extern from "<sys/uio.h>" nogil:

    cdef struct iovec:
        void  *iov_base
        size_t iov_len

    ssize_t readv (int fd, const iovec *iov, int iovcnt)
    ssize_t writev(int fd, const iovec *iov, int iovcnt)

    # Linux-specific, https://man7.org/linux/man-pages/man2/readv.2.html
    ssize_t preadv (int fd, const iovec *iov, int iovcnt, off_t offset)
    ssize_t pwritev(int fd, const iovec *iov, int iovcnt, off_t offset)

    enum: RWF_DSYNC
    enum: RWF_HIPRI
    enum: RWF_SYNC
    enum: RWF_NOWAIT
    enum: RWF_APPEND

    ssize_t preadv2 (int fd, const iovec *iov, int iovcnt, off_t offset, int flags)
    ssize_t pwritev2(int fd, const iovec *iov, int iovcnt, off_t offset, int flags)

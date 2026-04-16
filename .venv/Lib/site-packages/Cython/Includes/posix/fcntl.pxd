# https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/fcntl.h.html

from posix.types cimport mode_t, off_t, pid_t

cdef extern from "<fcntl.h>" nogil:

    enum: F_DUPFD
    enum: F_DUPFD_CLOEXEC
    enum: F_GETFD
    enum: F_SETFD
    enum: F_GETFL
    enum: F_SETFL
    enum: F_GETLK
    enum: F_SETLK
    enum: F_SETLKW
    enum: F_GETOWN
    enum: F_SETOWN

    enum: FD_CLOEXEC

    enum: F_RDLCK
    enum: F_UNLCK
    enum: F_WRLCK

    enum: SEEK_SET
    enum: SEEK_CUR
    enum: SEEK_END

    enum: O_CLOEXEC
    enum: O_CREAT
    enum: O_DIRECT
    enum: O_DIRECTORY
    enum: O_EXCL
    enum: O_NOCTTY
    enum: O_TRUNC
    enum: O_TTY_INIT

    enum: O_APPEND
    enum: O_DSYNC
    enum: O_NONBLOCK
    enum: O_RSYNC
    enum: O_SYNC

    enum: O_ACCMODE # O_RDONLY|O_WRONLY|O_RDWR

    enum: O_EXEC
    enum: O_RDONLY
    enum: O_WRONLY
    enum: O_RDWR
    enum: O_SEARCH

    enum: AT_FDCWD
    enum: AT_EACCESS
    enum: AT_SYMLINK_NOFOLLOW
    enum: AT_SYMLINK_FOLLOW
    enum: AT_REMOVEDIR

    enum: S_IFMT
    enum: S_IFBLK
    enum: S_IFCHR
    enum: S_IFIFO
    enum: S_IFREG
    enum: S_IFDIR
    enum: S_IFLNK
    enum: S_IFSOCK

    enum: POSIX_FADV_DONTNEED
    enum: POSIX_FADV_NOREUSE
    enum: POSIX_FADV_NORMAL
    enum: POSIX_FADV_RANDOM
    enum: POSIX_FADV_SEQUENTIAL
    enum: POSIX_FADV_WILLNEED

    struct flock:
        short l_type
        short l_whence
        off_t l_start
        off_t l_len
        pid_t l_pid

    int creat(const char *, mode_t)
    int fcntl(int, int, ...)
    int open(const char *, int, ...)
    int openat(int, const char *, int, ...)
    int posix_fadvise(int, off_t, off_t, int)
    int posix_fallocate(int, off_t, off_t)

# https://pubs.opengroup.org/onlinepubs/009695399/basedefs/sys/stat.h.html
# https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/sys_stat.h.html

from posix.types cimport (blkcnt_t, blksize_t, dev_t, gid_t, ino_t, mode_t,
                          nlink_t, off_t, time_t, uid_t)
from posix.time cimport timespec


cdef extern from "<sys/stat.h>" nogil:
    cdef struct struct_stat "stat":
        dev_t   st_dev
        ino_t   st_ino
        mode_t  st_mode
        nlink_t st_nlink
        uid_t   st_uid
        gid_t   st_gid
        dev_t   st_rdev
        off_t   st_size
        blksize_t st_blksize
        blkcnt_t st_blocks
        # POSIX.1-2001
        time_t  st_atime
        time_t  st_mtime
        time_t  st_ctime
        # POSIX.1-2008
        timespec st_atim
        timespec st_mtim
        timespec st_ctim

        # st_birthtime exists on *BSD and OS X.
        # Under Linux, defining it here does not hurt. Compilation under Linux
        # will only (and rightfully) fail when attempting to use the field.
        time_t  st_birthtime

# POSIX prescribes including both <sys/stat.h> and <unistd.h> for these
cdef extern from "<unistd.h>" nogil:
    int chmod(const char *, mode_t)
    int fchmod(int, mode_t)
    int fchmodat(int, const char *, mode_t, int flags)

    int stat(const char *, struct_stat *)
    int lstat(const char *, struct_stat *)
    int fstat(int, struct_stat *)
    int fstatat(int, const char *, struct_stat *, int flags)

    int mkdir(const char *, mode_t)
    int mkdirat(int, const char *, mode_t)
    int mkfifo(const char *, mode_t)
    int mkfifoat(int, const char *, mode_t)
    int mknod(const char *, mode_t, dev_t)
    int mknodat(int, const char *, mode_t, dev_t)

    int futimens(int, const timespec *)
    int utimensat(int, const char *, const timespec *, int flags)

    # Macros for st_mode
    mode_t S_ISREG(mode_t)
    mode_t S_ISDIR(mode_t)
    mode_t S_ISCHR(mode_t)
    mode_t S_ISBLK(mode_t)
    mode_t S_ISFIFO(mode_t)
    mode_t S_ISLNK(mode_t)
    mode_t S_ISSOCK(mode_t)

    mode_t S_IFMT
    mode_t S_IFREG
    mode_t S_IFDIR
    mode_t S_IFCHR
    mode_t S_IFBLK
    mode_t S_IFIFO
    mode_t S_IFLNK
    mode_t S_IFSOCK

    # Permissions
    mode_t S_ISUID
    mode_t S_ISGID
    mode_t S_ISVTX

    mode_t S_IRWXU
    mode_t S_IRUSR
    mode_t S_IWUSR
    mode_t S_IXUSR

    mode_t S_IRWXG
    mode_t S_IRGRP
    mode_t S_IWGRP
    mode_t S_IXGRP

    mode_t S_IRWXO
    mode_t S_IROTH
    mode_t S_IWOTH
    mode_t S_IXOTH

    # test file types
    bint S_TYPEISMQ(struct_stat *buf)
    bint S_TYPEISSEM(struct_stat *buf)
    bint S_TYPEISSHM(struct_stat *buf)
    bint S_TYPEISTMO(struct_stat *buf)

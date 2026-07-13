-Pinned :class:`cython` to `<3.3` in the build requirements to avoid a
compiler crash on fused types with const/pointer qualifier, which were introduced in cython 3.3.0a2.
By :user:`imann128`

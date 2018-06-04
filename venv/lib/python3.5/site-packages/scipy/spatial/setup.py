from __future__ import division, print_function, absolute_import

from os.path import join, dirname
import glob


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from numpy.distutils.misc_util import get_info as get_misc_info
    from numpy.distutils.system_info import get_info as get_sys_info
    from distutils.sysconfig import get_python_inc

    config = Configuration('spatial', parent_package, top_path)

    config.add_data_dir('tests')

    # qhull
    qhull_src = list(glob.glob(join(dirname(__file__), 'qhull',
                                    'src', '*.c')))

    inc_dirs = [get_python_inc()]
    if inc_dirs[0] != get_python_inc(plat_specific=1):
        inc_dirs.append(get_python_inc(plat_specific=1))
    inc_dirs.append(get_numpy_include_dirs())
    inc_dirs.append(join(dirname(dirname(__file__)), '_lib'))

    cfg = dict(get_sys_info('lapack_opt'))
    cfg.setdefault('include_dirs', []).extend(inc_dirs)
    config.add_extension('qhull',
                         sources=['qhull.c'] + qhull_src,
                         **cfg)

    # cKDTree
    ckdtree_src = ['query.cxx',
                   'build.cxx',
                   'globals.cxx',
                   'cpp_exc.cxx',
                   'query_pairs.cxx',
                   'count_neighbors.cxx',
                   'query_ball_point.cxx',
                   'query_ball_tree.cxx',
                   'sparse_distances.cxx']

    ckdtree_src = [join('ckdtree', 'src', x) for x in ckdtree_src]

    ckdtree_headers = ['ckdtree_decl.h',
                       'ckdtree_methods.h',
                       'coo_entries.h',
                       'cpp_exc.h',
                       'cpp_utils.h',
                       'distance_base.h',
                       'distance.h',
                       'ordered_pair.h',
                       'partial_sort.h',
                       'rectangle.h']

    ckdtree_headers = [join('ckdtree', 'src', x) for x in ckdtree_headers]

    ckdtree_dep = ['ckdtree.cxx'] + ckdtree_headers + ckdtree_src
    config.add_extension('ckdtree',
                         sources=['ckdtree.cxx'] + ckdtree_src,
                         depends=ckdtree_dep,
                         include_dirs=inc_dirs + [join('ckdtree', 'src')])
    # _distance_wrap
    config.add_extension('_distance_wrap',
                         sources=[join('src', 'distance_wrap.c')],
                         depends=[join('src', 'distance_impl.h')],
                         include_dirs=[get_numpy_include_dirs()],
                         extra_info=get_misc_info("npymath"))

    config.add_extension('_voronoi',
                         sources=['_voronoi.c'])

    config.add_extension('_hausdorff',
                         sources=['_hausdorff.c'])

    # Add license files
    config.add_data_files('qhull/COPYING.txt')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

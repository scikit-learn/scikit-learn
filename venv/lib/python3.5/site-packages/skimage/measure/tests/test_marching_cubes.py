import numpy as np
from skimage.draw import ellipsoid, ellipsoid_stats
from skimage.measure import (marching_cubes_classic, marching_cubes_lewiner,
                             mesh_surface_area, correct_mesh_orientation)
from skimage._shared import testing
from skimage._shared.testing import assert_array_equal


def test_marching_cubes_isotropic():
    ellipsoid_isotropic = ellipsoid(6, 10, 16, levelset=True)
    _, surf = ellipsoid_stats(6, 10, 16)
    
    # Classic
    verts, faces = marching_cubes_classic(ellipsoid_isotropic, 0.)
    surf_calc = mesh_surface_area(verts, faces)
    # Test within 1% tolerance for isotropic. Will always underestimate.
    assert surf > surf_calc and surf_calc > surf * 0.99
    
    # Lewiner
    verts, faces = marching_cubes_lewiner(ellipsoid_isotropic, 0.)[:2]
    surf_calc = mesh_surface_area(verts, faces)
    # Test within 1% tolerance for isotropic. Will always underestimate.
    assert surf > surf_calc and surf_calc > surf * 0.99


def test_marching_cubes_anisotropic():
    # test spacing as numpy array (and not just tuple)
    spacing = np.array([1., 10 / 6., 16 / 6.])
    ellipsoid_anisotropic = ellipsoid(6, 10, 16, spacing=spacing,
                                      levelset=True)
    _, surf = ellipsoid_stats(6, 10, 16)
    
    # Classic
    verts, faces = marching_cubes_classic(ellipsoid_anisotropic, 0.,
                                          spacing=spacing)
    surf_calc = mesh_surface_area(verts, faces)
    # Test within 1.5% tolerance for anisotropic. Will always underestimate.
    assert surf > surf_calc and surf_calc > surf * 0.985
    
    # Lewiner
    verts, faces = marching_cubes_lewiner(
        ellipsoid_anisotropic, 0., spacing=spacing)[:2]
    surf_calc = mesh_surface_area(verts, faces)
    # Test within 1.5% tolerance for anisotropic. Will always underestimate.
    assert surf > surf_calc and surf_calc > surf * 0.985

    # Test spacing together with allow_degenerate=False
    marching_cubes_lewiner(ellipsoid_anisotropic, 0, spacing=spacing,
                           allow_degenerate=False)


def test_invalid_input():
    # Classic
    with testing.raises(ValueError):
        marching_cubes_classic(np.zeros((2, 2, 1)), 0)
    with testing.raises(ValueError):
        marching_cubes_classic(np.zeros((2, 2, 1)), 1)
    with testing.raises(ValueError):
        marching_cubes_classic(np.ones((3, 3, 3)), 1, spacing=(1, 2))
    with testing.raises(ValueError):
        marching_cubes_classic(np.zeros((20, 20)), 0)
    
    # Lewiner
    with testing.raises(ValueError):
        marching_cubes_lewiner(np.zeros((2, 2, 1)), 0)
    with testing.raises(ValueError):
        marching_cubes_lewiner(np.zeros((2, 2, 1)), 1)
    with testing.raises(ValueError):
        marching_cubes_lewiner(np.ones((3, 3, 3)), 1,
                               spacing=(1, 2))
    with testing.raises(ValueError):
        marching_cubes_lewiner(np.zeros((20, 20)), 0)


def test_correct_mesh_orientation():
    sphere_small = ellipsoid(1, 1, 1, levelset=True)

    # Mesh with incorrectly oriented faces which was previously returned from
    # `marching_cubes`, before it guaranteed correct mesh orientation
    verts = np.array([[1., 2., 2.],
                      [2., 2., 1.],
                      [2., 1., 2.],
                      [2., 2., 3.],
                      [2., 3., 2.],
                      [3., 2., 2.]])

    faces = np.array([[0, 1, 2],
                      [2, 0, 3],
                      [1, 0, 4],
                      [4, 0, 3],
                      [1, 2, 5],
                      [2, 3, 5],
                      [1, 4, 5],
                      [5, 4, 3]])

    # Correct mesh orientation - descent
    corrected_faces1 = correct_mesh_orientation(sphere_small, verts, faces,
                                                gradient_direction='descent')
    corrected_faces2 = correct_mesh_orientation(sphere_small, verts, faces,
                                                gradient_direction='ascent')

    # Ensure ascent is opposite of descent for all faces
    assert_array_equal(corrected_faces1, corrected_faces2[:, ::-1])

    # Ensure correct faces have been reversed: 1, 4, and 5
    idx = [1, 4, 5]
    expected = faces.copy()
    expected[idx] = expected[idx, ::-1]
    assert_array_equal(expected, corrected_faces1)


def test_both_algs_same_result_ellipse():
    # Performing this test on data that does not have ambiguities
    
    sphere_small = ellipsoid(1, 1, 1, levelset=True)
    
    vertices1, faces1 = marching_cubes_classic(sphere_small, 0)[:2]
    vertices2, faces2 = marching_cubes_lewiner(
        sphere_small, 0, allow_degenerate=False)[:2]
    vertices3, faces3 = marching_cubes_lewiner(
        sphere_small, 0, allow_degenerate=False, use_classic=True)[:2]
    
    # Order is different, best we can do is test equal shape and same vertices present
    assert _same_mesh(vertices1, faces1, vertices2, faces2)
    assert _same_mesh(vertices1, faces1, vertices3, faces3)


def _same_mesh(vertices1, faces1, vertices2, faces2, tol=1e-10):
    """ Compare two meshes, using a certain tolerance and invariant to
    the order of the faces.
    """
    # Unwind vertices
    triangles1 = vertices1[np.array(faces1)]
    triangles2 = vertices2[np.array(faces2)]
    # Sort vertices within each triangle
    triang1 = [np.concatenate(sorted(t, key=lambda x:tuple(x)))
               for t in triangles1]
    triang2 = [np.concatenate(sorted(t, key=lambda x:tuple(x)))
               for t in triangles2]
    # Sort the resulting 9-element "tuples"
    triang1 = np.array(sorted([tuple(x) for x in triang1]))
    triang2 = np.array(sorted([tuple(x) for x in triang2]))
    return (triang1.shape == triang2.shape and
            np.allclose(triang1, triang2, 0, tol))


def test_both_algs_same_result_donut():
    # Performing this test on data that does not have ambiguities
    n = 48
    a, b = 2.5/n, -1.25

    vol = np.empty((n, n, n), 'float32')
    for iz in range(vol.shape[0]):
        for iy in range(vol.shape[1]):
            for ix in range(vol.shape[2]):
                # Double-torii formula by Thomas Lewiner
                z, y, x = float(iz)*a+b, float(iy)*a+b, float(ix)*a+b
                vol[iz,iy,ix] = ( ( 
                    (8*x)**2 + (8*y-2)**2 + (8*z)**2 + 16 - 1.85*1.85 ) * ( (8*x)**2 +
                    (8*y-2)**2 + (8*z)**2 + 16 - 1.85*1.85 ) - 64 * ( (8*x)**2 + (8*y-2)**2 )
                    ) * ( ( (8*x)**2 + ((8*y-2)+4)*((8*y-2)+4) + (8*z)**2 + 16 - 1.85*1.85 )
                    * ( (8*x)**2 + ((8*y-2)+4)*((8*y-2)+4) + (8*z)**2 + 16 - 1.85*1.85 ) -
                    64 * ( ((8*y-2)+4)*((8*y-2)+4) + (8*z)**2 
                    ) ) + 1025
    
    vertices1, faces1 = marching_cubes_classic(vol, 0)[:2]
    vertices2, faces2 = marching_cubes_lewiner(vol, 0)[:2]
    vertices3, faces3 = marching_cubes_lewiner(vol, 0, use_classic=True)[:2]
    
    # Old and new alg are different
    assert not _same_mesh(vertices1, faces1, vertices2, faces2)
    # New classic and new Lewiner are different
    assert not _same_mesh(vertices2, faces2, vertices3, faces3)
    # Would have been nice if old and new classic would have been the same
    # assert _same_mesh(vertices1, faces1, vertices3, faces3, 5)

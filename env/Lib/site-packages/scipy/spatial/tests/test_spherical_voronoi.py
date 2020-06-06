from __future__ import print_function
import numpy as np
import itertools
from numpy.testing import (assert_equal,
                           assert_almost_equal,
                           assert_array_equal,
                           assert_array_almost_equal)
import pytest
from pytest import raises as assert_raises
from pytest import warns as assert_warns
from scipy.spatial import SphericalVoronoi, distance
from scipy.spatial import _spherical_voronoi as spherical_voronoi
from scipy._lib._numpy_compat import suppress_warnings
from scipy.spatial.transform import Rotation
from scipy.optimize import linear_sum_assignment


TOL = 1E-10


class TestSphericalVoronoi(object):

    def setup_method(self):
        self.points = np.array([
            [-0.78928481, -0.16341094, 0.59188373],
            [-0.66839141, 0.73309634, 0.12578818],
            [0.32535778, -0.92476944, -0.19734181],
            [-0.90177102, -0.03785291, -0.43055335],
            [0.71781344, 0.68428936, 0.12842096],
            [-0.96064876, 0.23492353, -0.14820556],
            [0.73181537, -0.22025898, -0.6449281],
            [0.79979205, 0.54555747, 0.25039913]]
        )

        # Issue #9386
        self.hemisphere_points = np.array([
            [0.88610999, -0.42383021, 0.18755541],
            [0.51980039, -0.72622668, 0.4498915],
            [0.56540011, -0.81629197, -0.11827989],
            [0.69659682, -0.69972598, 0.15854467]])

        # Issue #8859
        phi = np.linspace(0, 2 * np.pi, 10, endpoint=False)    # azimuth angle
        theta = np.linspace(0.001, np.pi * 0.4, 5)    # polar angle
        theta = theta[np.newaxis, :].T

        phiv, thetav = np.meshgrid(phi, theta)
        phiv = np.reshape(phiv, (50, 1))
        thetav = np.reshape(thetav, (50, 1))

        x = np.cos(phiv) * np.sin(thetav)
        y = np.sin(phiv) * np.sin(thetav)
        z = np.cos(thetav)
        self.hemisphere_points2 = np.concatenate([x, y, z], axis=1)

    def test_constructor(self):
        center = np.array([1, 2, 3])
        radius = 2
        s1 = SphericalVoronoi(self.points)
        # user input checks in SphericalVoronoi now require
        # the radius / center to match the generators so adjust
        # accordingly here
        s2 = SphericalVoronoi(self.points * radius, radius)
        s3 = SphericalVoronoi(self.points + center, center=center)
        s4 = SphericalVoronoi(self.points * radius + center, radius, center)
        assert_array_equal(s1.center, np.array([0, 0, 0]))
        assert_equal(s1.radius, 1)
        assert_array_equal(s2.center, np.array([0, 0, 0]))
        assert_equal(s2.radius, 2)
        assert_array_equal(s3.center, center)
        assert_equal(s3.radius, 1)
        assert_array_equal(s4.center, center)
        assert_equal(s4.radius, radius)

    def test_vertices_regions_translation_invariance(self):
        sv_origin = SphericalVoronoi(self.points)
        center = np.array([1, 1, 1])
        sv_translated = SphericalVoronoi(self.points + center, center=center)
        assert_array_equal(sv_origin.regions, sv_translated.regions)
        assert_array_almost_equal(sv_origin.vertices + center,
                                  sv_translated.vertices)

    def test_vertices_regions_scaling_invariance(self):
        sv_unit = SphericalVoronoi(self.points)
        sv_scaled = SphericalVoronoi(self.points * 2, 2)
        assert_array_equal(sv_unit.regions, sv_scaled.regions)
        assert_array_almost_equal(sv_unit.vertices * 2,
                                  sv_scaled.vertices)

    def test_old_radius_api(self):
        sv_unit = SphericalVoronoi(self.points, radius=1)
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning, "`radius` is `None`")
            sv = SphericalVoronoi(self.points, None)
            assert_array_almost_equal(sv_unit.vertices, sv.vertices)

    def test_old_radius_api_warning(self):
        with assert_warns(DeprecationWarning):
            sv = SphericalVoronoi(self.points, None)

    def test_sort_vertices_of_regions(self):
        sv = SphericalVoronoi(self.points)
        unsorted_regions = sv.regions
        sv.sort_vertices_of_regions()
        assert_array_equal(sorted(sv.regions), sorted(unsorted_regions))

    def test_sort_vertices_of_regions_flattened(self):
        expected = sorted([[0, 6, 5, 2, 3], [2, 3, 10, 11, 8, 7], [0, 6, 4, 1],
                           [4, 8, 7, 5, 6], [9, 11, 10], [2, 7, 5],
                           [1, 4, 8, 11, 9], [0, 3, 10, 9, 1]])
        expected = list(itertools.chain(*sorted(expected)))
        sv = SphericalVoronoi(self.points)
        sv.sort_vertices_of_regions()
        actual = list(itertools.chain(*sorted(sv.regions)))
        assert_array_equal(actual, expected)

    def test_sort_vertices_of_regions_dimensionality(self):
        points = np.array([[1, 0, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1],
                           [0.5, 0.5, 0.5, 0.5]])
        with pytest.raises(TypeError, match="three-dimensional"):
            sv = spherical_voronoi.SphericalVoronoi(points)
            sv.sort_vertices_of_regions()

    def test_num_vertices(self):
        # for any n >= 3, a spherical Voronoi diagram has 2n - 4
        # vertices; this is a direct consequence of Euler's formula
        # as explained by Dinis and Mamede (2010) Proceedings of the
        # 2010 International Symposium on Voronoi Diagrams in Science
        # and Engineering
        sv = SphericalVoronoi(self.points)
        expected = self.points.shape[0] * 2 - 4
        actual = sv.vertices.shape[0]
        assert_equal(actual, expected)

    def test_voronoi_circles(self):
        sv = spherical_voronoi.SphericalVoronoi(self.points)
        for vertex in sv.vertices:
            distances = distance.cdist(sv.points, np.array([vertex]))
            closest = np.array(sorted(distances)[0:3])
            assert_almost_equal(closest[0], closest[1], 7, str(vertex))
            assert_almost_equal(closest[0], closest[2], 7, str(vertex))

    def test_duplicate_point_handling(self):
        # an exception should be raised for degenerate generators
        # related to Issue# 7046
        self.degenerate = np.concatenate((self.points, self.points))
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.degenerate)

    def test_incorrect_radius_handling(self):
        # an exception should be raised if the radius provided
        # cannot possibly match the input generators
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.points,
                                                    radius=0.98)

    def test_incorrect_center_handling(self):
        # an exception should be raised if the center provided
        # cannot possibly match the input generators
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.points,
                                                    center=[0.1, 0, 0])

    def test_single_hemisphere_handling(self):
        # Test solution of Issues #9386, #8859

        for points in [self.hemisphere_points, self.hemisphere_points2]:
            sv = SphericalVoronoi(points)
            triangles = sv._tri.points[sv._tri.simplices]
            dots = np.einsum('ij,ij->i', sv.vertices, triangles[:, 0])
            circumradii = np.arccos(np.clip(dots, -1, 1))
            assert np.max(circumradii) > np.pi / 2

    def test_rank_deficient(self):
        # rank-1 input cannot be triangulated
        points = np.array([[-1, 0, 0], [1, 0, 0]])
        with pytest.raises(ValueError, match="Rank of input points"):
            sv = spherical_voronoi.SphericalVoronoi(points)

    @pytest.mark.parametrize("n", [8, 15, 21])
    @pytest.mark.parametrize("radius", [0.5, 1, 2])
    @pytest.mark.parametrize("center", [(0, 0, 0), (1, 2, 3)])
    def test_geodesic_input(self, n, radius, center):
        U = Rotation.random(random_state=0).as_matrix()
        thetas = np.linspace(0, 2 * np.pi, n, endpoint=False)
        points = np.vstack([np.sin(thetas), np.cos(thetas), np.zeros(n)]).T
        points = radius * points @ U
        sv = SphericalVoronoi(points + center, radius=radius, center=center)

        # each region must have 4 vertices
        region_sizes = np.array([len(region) for region in sv.regions])
        assert (region_sizes == 4).all()
        regions = np.array(sv.regions)

        # vertices are those between each pair of input points + north and
        # south poles
        vertices = sv.vertices - center
        assert len(vertices) == n + 2

        # verify that north and south poles are orthogonal to geodesic on which
        # input points lie
        poles = vertices[n:]
        assert np.abs(np.dot(points, poles.T)).max() < 1E-10

        for point, region in zip(points, sv.regions):
            cosine = np.dot(vertices[region], point)
            sine = np.linalg.norm(np.cross(vertices[region], point), axis=1)
            arclengths = radius * np.arctan2(sine, cosine)
            # test arc lengths to poles
            assert_almost_equal(arclengths[[1, 3]], radius * np.pi / 2)
            # test arc lengths to forward and backward neighbors
            assert_almost_equal(arclengths[[0, 2]], radius * np.pi / n)

        regions = sv.regions.copy()
        sv.sort_vertices_of_regions()
        assert regions == sv.regions

    @pytest.mark.parametrize("dim", range(2, 7))
    def test_higher_dimensions(self, dim):
        n = 100
        rng = np.random.RandomState(seed=0)
        points = rng.randn(n, dim)
        points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
        sv = SphericalVoronoi(points)
        assert sv.vertices.shape[1] == dim
        assert len(sv.regions) == n

        # verify Euler characteristic
        cell_counts = []
        simplices = np.sort(sv._tri.simplices)
        for i in range(1, dim + 1):
            cells = []
            for indices in itertools.combinations(range(dim), i):
                cells.append(simplices[:, list(indices)])
            cells = np.unique(np.concatenate(cells), axis=0)
            cell_counts.append(len(cells))
        expected_euler = 1 + (-1)**(dim-1)
        actual_euler = sum([(-1)**i * e for i, e in enumerate(cell_counts)])
        assert expected_euler == actual_euler

    @pytest.mark.parametrize("dim", range(2, 7))
    def test_cross_polytope_regions(self, dim):
        # The hypercube is the dual of the cross-polytope, so the voronoi
        # vertices of the cross-polytope lie on the points of the hypercube.

        # generate points of the cross-polytope
        points = np.concatenate((-np.eye(dim), np.eye(dim)))
        sv = SphericalVoronoi(points)
        assert all([len(e) == 2**(dim - 1) for e in sv.regions])

        # generate points of the hypercube
        expected = np.vstack(list(itertools.product([-1, 1], repeat=dim)))
        expected = expected.astype(np.float) / np.sqrt(dim)

        # test that Voronoi vertices are correctly placed
        dist = distance.cdist(sv.vertices, expected)
        res = linear_sum_assignment(dist)
        assert dist[res].sum() < TOL

    @pytest.mark.parametrize("dim", range(2, 4))
    def test_hypercube_regions(self, dim):
        # The cross-polytope is the dual of the hypercube, so the voronoi
        # vertices of the hypercube lie on the points of the cross-polytope.

        # generate points of the hypercube
        points = np.vstack(list(itertools.product([-1, 1], repeat=dim)))
        points = points.astype(np.float) / np.sqrt(dim)
        sv = SphericalVoronoi(points)

        # generate points of the cross-polytope
        expected = np.concatenate((-np.eye(dim), np.eye(dim)))

        # test that Voronoi vertices are correctly placed
        dist = distance.cdist(sv.vertices, expected)
        res = linear_sum_assignment(dist)
        assert dist[res].sum() < TOL

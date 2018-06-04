from __future__ import print_function
import numpy as np
import itertools
from numpy.testing import (assert_equal,
                           assert_almost_equal,
                           assert_array_equal,
                           assert_array_almost_equal)
from pytest import raises as assert_raises
from scipy.spatial import SphericalVoronoi, distance
from scipy.spatial import _spherical_voronoi as spherical_voronoi


class TestCircumcenters(object):

    def test_circumcenters(self):
        tetrahedrons = np.array([
            [[1, 2, 3],
             [-1.1, -2.1, -3.1],
             [-1.2, 2.2, 3.2],
             [-1.3, -2.3, 3.3]],
            [[10, 20, 30],
             [-10.1, -20.1, -30.1],
             [-10.2, 20.2, 30.2],
             [-10.3, -20.3, 30.3]]
        ])

        result = spherical_voronoi.calc_circumcenters(tetrahedrons)

        expected = [
            [-0.5680861153262529, -0.133279590288315, 0.1843323216995444],
            [-0.5965330784014926, -0.1480377040397778, 0.1981967854886021]
        ]

        assert_array_almost_equal(result, expected)


class TestProjectToSphere(object):

    def test_unit_sphere(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([0, 0, 0])
        radius = 1
        projected = spherical_voronoi.project_to_sphere(points, center, radius)
        assert_array_almost_equal(points, projected)

    def test_scaled_points(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([0, 0, 0])
        radius = 1
        scaled = points * 2
        projected = spherical_voronoi.project_to_sphere(scaled, center, radius)
        assert_array_almost_equal(points, projected)

    def test_translated_sphere(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([1, 2, 3])
        translated = points + center
        radius = 1
        projected = spherical_voronoi.project_to_sphere(translated, center,
                                                        radius)
        assert_array_almost_equal(translated, projected)


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

    def test_constructor(self):
        center = np.array([1, 2, 3])
        radius = 2
        s1 = SphericalVoronoi(self.points)
        # user input checks in SphericalVoronoi now require
        # the radius / center to match the generators so adjust
        # accordingly here
        s2 = SphericalVoronoi(self.points * radius, radius)
        s3 = SphericalVoronoi(self.points + center, None, center)
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
        sv_translated = SphericalVoronoi(self.points + center, None, center)
        assert_array_equal(sv_origin.regions, sv_translated.regions)
        assert_array_almost_equal(sv_origin.vertices + center,
                                  sv_translated.vertices)

    def test_vertices_regions_scaling_invariance(self):
        sv_unit = SphericalVoronoi(self.points)
        sv_scaled = SphericalVoronoi(self.points * 2, 2)
        assert_array_equal(sv_unit.regions, sv_scaled.regions)
        assert_array_almost_equal(sv_unit.vertices * 2,
                                  sv_scaled.vertices)

    def test_sort_vertices_of_regions(self):
        sv = SphericalVoronoi(self.points)
        unsorted_regions = sv.regions
        sv.sort_vertices_of_regions()
        assert_array_equal(sorted(sv.regions), sorted(unsorted_regions))

    def test_sort_vertices_of_regions_flattened(self):
        expected = sorted([[0, 6, 5, 2, 3], [2, 3, 10, 11, 8, 7], [0, 6, 4, 1], [4, 8,
            7, 5, 6], [9, 11, 10], [2, 7, 5], [1, 4, 8, 11, 9], [0, 3, 10, 9,
                1]])
        expected = list(itertools.chain(*sorted(expected)))
        sv = SphericalVoronoi(self.points)
        sv.sort_vertices_of_regions()
        actual = list(itertools.chain(*sorted(sv.regions)))
        assert_array_equal(actual, expected)

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
            distances = distance.cdist(sv.points,np.array([vertex]))
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
                                                    center=[0.1,0,0])

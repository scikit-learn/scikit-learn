#!/usr/bin/env python

import unittest
import numpy

from unittest import TestCase
from ..barycenter import Barycenter

class TestBarycenter(TestCase):
  def test_create(self):
    test = Barycenter()

if __name__ == "__main__":
  unittest.main()
  
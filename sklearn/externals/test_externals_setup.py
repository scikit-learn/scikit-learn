"""
Fixtures to get the external bundled dependencies tested.

This module gets loaded by test discovery scanners (such as nose) in
their collection scan.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

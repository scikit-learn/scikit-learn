"""
Auto Cython Optimization Analyzer for scikit-learn

This script scans all .pyx Cython files in the repository,
analyzes them for static typing usage, presence of annotations,
NumPy optimizations, OpenMP directives, and generates a basic
report as HTML for performance optimization insights.
"""

import os
import re
from pathlib import Path
from collections import defaultdict

CYTHON_KEYWORDS = {
    'static_typing': ['cdef int', 'cdef double', 'cdef float', 'cdef long', 'cdef unsigned', 'cpdef'],
    'numpy_usage': ['numpy.ndarray', 'np.ndarray', 'cimport numpy'],
    'openmp': ['nogil', 'prange'],
    'annotations': ['# cython: profile=True', '# cython: boundscheck=False', '# cython: wraparound=False']
}

REPORT_HTML = "cython_optimization_report.html"

def find_pyx_files(base_dir):
    return list(Path(base_dir).rglob("*.pyx"))

def analyze_file(filepath):
    results = defaultdict(int)
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    for line in lines:
        for key, patterns in CYTHON_KEYWORDS.items():
            for pattern in patterns:
                if pattern in line:
                    results[key] += 1
    return results

def generate_report(results_dict):
    html = ["<html><head><title>Cython Optimization Report</title></head><body>"]
    html.append("<h1>Cython Optimization Analysis Report</h1>")
    html.append("<table border='1'><tr><th>File</th><th>Static Typing</th><th>NumPy</th><th>OpenMP</th><th>Annotations</th></tr>")
    for file, metrics in results_dict.items():
        html.append(f"<tr><td>{file}</td><td>{metrics['static_typing']}</td><td>{metrics['numpy_usage']}</td><td>{metrics['openmp']}</td><td>{metrics['annotations']}</td></tr>")
    html.append("</table></body></html>")
    with open(REPORT_HTML, "w") as f:
        f.write("
".join(html))

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    pyx_files = find_pyx_files(base_dir)
    if not pyx_files:
        print("No .pyx files found.")
        return

    all_results = {}
    for filepath in pyx_files:
        metrics = analyze_file(filepath)
        all_results[str(filepath)] = metrics

    generate_report(all_results)
    print(f" Report generated: {REPORT_HTML}")

if __name__ == "__main__":
    main()

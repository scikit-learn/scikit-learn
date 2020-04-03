"""
Benchmarks of with vs without CSR conversion

The performances are evaluated in various input dimensions. The input array X is a n x n square matrix. 

Time is calculated by runing the transform 10 times respectively for each method. 

Performance Comparison
======================
Input_Dimensions Time_With_CSR_Conversion Time_Without_Conversion
-----------------------------------------------------------------
      1000               0.3052s                  0.2708s        
      2000               1.6780s                  1.7055s        
      3000               3.9639s                  4.0052s
      4000               7.5915s                  7.2595s
      5000               12.2266s                 12.1459s
      6000               18.1505s                 17.7042s
      7000               25.2273s                 25.0626s
      8000               33.1733s                 33.2724s
      9000               42.9401s                 42.7734s
"""
import scipy.sparse as ss
from sklearn.decomposition import TruncatedSVD
from time import time
from sklearn.utils import check_array
def benchmark_converted_transform(dim, trials=10):
    start_time = time()
    for i in range(trials):
        X = ss.random(dim, dim, format='csc')
        X = check_array(X, accept_sparse='csr')
        tsvd = TruncatedSVD().fit(X)
        tsvd.transform(X)
    return time() - start_time
def benchmark_direct_transform(dim, trials=10):
    start_time = time()
    for i in range(trials):
        X = ss.random(dim, dim, format='csc')
        tsvd = TruncatedSVD().fit(X)
        tsvd.transform(X)
    return time() - start_time
if __name__ == '__main__':
    print("Performance Comparison")
    print("======================")
    headers = ("Input_Dimensions", "Time_With_CSR_Conversion", "Time_Without_Conversion")
    header_len = [len(header) for header in headers]
    print("%s %s %s"
          % headers)
    print("-" * (sum(header_len) + 2))
    for dim in range(1000, 10000, 1000):
        print("%s %s %s" % ( ("%d" % dim).center(header_len[0]),
                                ("%.4fs" % benchmark_converted_transform(dim)).center(header_len[1]),
                                ("%.4fs" % benchmark_direct_transform(dim)).center(header_len[2])))
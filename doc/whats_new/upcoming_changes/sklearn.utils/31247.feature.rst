Feature
------------

- Added :func:`randomized_eigen_decomposition`, an approximate eigen decomposition method based on the approach from 
  **Halko, Martinsson, and Tropp (2011)**: *Finding Structure with Randomness: Probabilistic Algorithms for Constructing 
  Approximate Matrix Decompositions*. This method provides a faster alternative to existing eigen decomposition techniques.

- Integrated :func:`eigen_decomposition_one_pass` into :class:`sklearn.manifold.Isomap` and 
  :class:`sklearn.decomposition.KernelPCA` as an additional option for eigen decomposition.

- Added a test suite comparing the new method to existing solvers (:obj:`arpack`, :obj:`dense`, etc.), ensuring numerical 
  accuracy and stability.

- Included benchmarks to analyze the **performance improvements** over traditional eigen decomposition approaches.

by :user: `Sylvain Marié<@smarie>`, `Mohamed yaich<@yaichm>`, `Oussama Er-rabie<@eroussama>`, `Mohamed Dlimi<@Dlimim>`, 
`Hamza Zeroual<@HamzaLuffy>` and `Amine Hannoun<@AmineHannoun>`.
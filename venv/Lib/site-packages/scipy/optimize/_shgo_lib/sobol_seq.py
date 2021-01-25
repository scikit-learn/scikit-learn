"""
  Licensing:
    This code is distributed under the MIT license.


  Authors:
    Original FORTRAN77 version of i4_sobol by Bennett Fox.
    MATLAB version by John Burkardt.
    PYTHON version by Corrado Chisari

    Original Python version of is_prime by Corrado Chisari

    Original MATLAB versions of other functions by John Burkardt.
    PYTHON versions by Corrado Chisari

    Original code is available from
    http://people.sc.fsu.edu/~jburkardt/py_src/sobol/sobol.html

  Modifications:
    Wrapped into Python class [30.10.2017]
"""
import numpy as np

__all__ = ['Sobol']


class Sobol:
    def __init__(self):
        # Init class variables
        self.atmost = None
        self.dim_max = None
        self.dim_num_save = None
        self.initialized = None
        self.lastq = None
        self.log_max = None
        self.maxcol = None
        self.poly = None
        self.recipd = None
        self.seed_save = None
        self.v = None

    def i4_sobol_generate(self, dim_num, n, skip=1):
        """
        i4_sobol_generate generates a Sobol dataset.

        Parameters:
          Input, integer dim_num, the spatial dimension.
          Input, integer N, the number of points to generate.
          Input, integer SKIP, the number of initial points to skip.

          Output, real R(M,N), the points.
        """
        r = np.full((n, dim_num), np.nan)
        for j in range(n):
            seed = j + skip
            r[j, 0:dim_num], next_seed = self.i4_sobol(dim_num, seed)

        return r

    def i4_bit_hi1(self, n):
        """
        i4_bit_hi1 returns the position of the high 1 bit base 2 in an integer.

        Example:
          +------+-------------+-----
          |    N |      Binary | BIT
          +------|-------------+-----
          |    0 |           0 |   0
          |    1 |           1 |   1
          |    2 |          10 |   2
          |    3 |          11 |   2
          |    4 |         100 |   3
          |    5 |         101 |   3
          |    6 |         110 |   3
          |    7 |         111 |   3
          |    8 |        1000 |   4
          |    9 |        1001 |   4
          |   10 |        1010 |   4
          |   11 |        1011 |   4
          |   12 |        1100 |   4
          |   13 |        1101 |   4
          |   14 |        1110 |   4
          |   15 |        1111 |   4
          |   16 |       10000 |   5
          |   17 |       10001 |   5
          | 1023 |  1111111111 |  10
          | 1024 | 10000000000 |  11
          | 1025 | 10000000001 |  11

        Parameters:
          Input, integer N, the integer to be measured.
          N should be nonnegative. If N is nonpositive,
          the value will always be 0.

          Output, integer BIT, the number of bits base 2.
        """
        i = np.floor(n)
        bit = 0
        while i > 0:
            bit += 1
            i //= 2
        return bit

    def i4_bit_lo0(self, n):
        """
        I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.

        Example:
          +------+------------+----
          |    N |     Binary | BIT
          +------+------------+----
          |    0 |          0 |   1
          |    1 |          1 |   2
          |    2 |         10 |   1
          |    3 |         11 |   3
          |    4 |        100 |   1
          |    5 |        101 |   2
          |    6 |        110 |   1
          |    7 |        111 |   4
          |    8 |       1000 |   1
          |    9 |       1001 |   2
          |   10 |       1010 |   1
          |   11 |       1011 |   3
          |   12 |       1100 |   1
          |   13 |       1101 |   2
          |   14 |       1110 |   1
          |   15 |       1111 |   5
          |   16 |      10000 |   1
          |   17 |      10001 |   2
          | 1023 | 1111111111 |   1
          | 1024 | 0000000000 |   1
          | 1025 | 0000000001 |   1

        Parameters:
          Input, integer N, the integer to be measured.
          N should be nonnegative.

          Output, integer BIT, the position of the low 1 bit.
        """
        bit = 1
        i = np.floor(n)
        while i != 2 * (i // 2):
            bit += 1
            i //= 2
        return bit

    def i4_sobol(self, dim_num, seed):
        """
        i4_sobol generates a new quasirandom Sobol vector with each call.

        Discussion:
          The routine adapts the ideas of Antonov and Saleev.

        Reference:
          Antonov, Saleev,
          USSR Computational Mathematics and Mathematical Physics,
          Volume 19, 1980, pages 252 - 256.

          Paul Bratley, Bennett Fox,
          Algorithm 659:
          Implementing Sobol's Quasirandom Sequence Generator,
          ACM Transactions on Mathematical Software,
          Volume 14, Number 1, pp. 88-100, 1988.

          Bennett Fox,
          Algorithm 647:
          Implementation and Relative Efficiency of Quasirandom
          Sequence Generators,
          ACM Transactions on Mathematical Software,
          Volume 12, Number 4, pp. 362-376, 1986.

          Ilya Sobol,
          USSR Computational Mathematics and Mathematical Physics,
          Volume 16, pp. 236-242, 1977.

          Ilya Sobol, Levitan,
          The Production of Points Uniformly Distributed in a Multidimensional
          Cube (in Russian),
          Preprint IPM Akad. Nauk SSSR,
          Number 40, Moscow 1976.

        Parameters:
          Input, integer DIM_NUM, the number of spatial dimensions.
          DIM_NUM must satisfy 1 <= DIM_NUM <= 40.

          Input/output, integer SEED, the "seed" for the sequence.
          This is essentially the index in the sequence of the quasirandom
          value to be generated. On output, SEED has been set to the
          appropriate next value, usually simply SEED+1.
          If SEED is less than 0 on input, it is treated as though it were 0.
          An input value of 0 requests the first (0th) element of the sequence.

          Output, real QUASI(DIM_NUM), the next quasirandom vector.
        """

        # if 'self.initialized' not in list(globals().keys()):
        if self.initialized is None:
            self.initialized = 0
            self.dim_num_save = -1

        if not self.initialized or dim_num != self.dim_num_save:
            self.initialized = 1
            self.dim_max = 40
            self.dim_num_save = -1
            self.log_max = 30
            self.seed_save = -1

            #  Initialize (part of) V.
            self.v = np.zeros((self.dim_max, self.log_max))
            self.v[0:40, 0] = np.transpose([
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

            self.v[2:40, 1] = np.transpose([
                1, 3, 1, 3, 1, 3, 3, 1,
                3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
                1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
                3, 1, 1, 3, 1, 3, 1, 3, 1, 3])

            self.v[3:40, 2] = np.transpose([
                7, 5, 1, 3, 3, 7, 5,
                5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
                5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
                5, 1, 1, 5, 7, 7, 5, 1, 3, 3])

            self.v[5:40, 3] = np.transpose([
                1, 7, 9, 13, 11,
                1, 3, 7, 9, 5, 13, 13, 11, 3, 15,
                5, 3, 15, 7, 9, 13, 9, 1, 11, 7,
                5, 15, 1, 15, 11, 5, 3, 1, 7, 9])

            self.v[7:40, 4] = np.transpose([
                9, 3, 27,
                15, 29, 21, 23, 19, 11, 25, 7, 13, 17,
                1, 25, 29, 3, 31, 11, 5, 23, 27, 19,
                21, 5, 1, 17, 13, 7, 15, 9, 31, 9])

            self.v[13:40, 5] = np.transpose([
                37, 33, 7, 5, 11, 39, 63,
                27, 17, 15, 23, 29, 3, 21, 13, 31, 25,
                9, 49, 33, 19, 29, 11, 19, 27, 15, 25])

            self.v[19:40, 6] = np.transpose([
                13,
                33, 115, 41, 79, 17, 29, 119, 75, 73, 105,
                7, 59, 65, 21, 3, 113, 61, 89, 45, 107])

            self.v[37:40, 7] = np.transpose([
                7, 23, 39])

            #  Set POLY.
            self.poly = [
                1, 3, 7, 11, 13, 19, 25, 37, 59, 47,
                61, 55, 41, 67, 97, 91, 109, 103, 115, 131,
                193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
                213, 191, 253, 203, 211, 239, 247, 285, 369, 299]

            self.atmost = 2 ** self.log_max - 1

            #  Find the number of bits in ATMOST.
            self.maxcol = self.i4_bit_hi1(self.atmost)

            #  Initialize row 1 of V.
            self.v[0, 0:self.maxcol] = 1

        # Things to do only if the dimension changed.
        if dim_num != self.dim_num_save:
            self.dim_num_save = dim_num

            #  Initialize the remaining rows of V.
            for i in range(2, dim_num + 1):

                # The bits of the integer POLY(I) gives the form of
                # self.polynomial I.
                # Find the degree of self.polynomial I from binary encoding.
                j = self.poly[i - 1]
                m = 0
                j //= 2
                while j > 0:
                    j //= 2
                    m += 1

                # Expand this bit pattern to separate
                # components of the logical array INCLUD.
                j = self.poly[i - 1]
                includ = np.empty(m)
                for k in range(m, 0, -1):
                    j2 = j // 2
                    includ[k - 1] = (j != 2 * j2)
                    j = j2

                # Calculate the remaining elements of row I as explained
                #  in Bratley and Fox, section 2.
                for j in range(m + 1, self.maxcol + 1):
                    newv = self.v[i - 1, j - m - 1]
                    lseed = 1
                    for k in range(1, m + 1):
                        lseed *= 2
                        if includ[k - 1]:
                            newv = np.bitwise_xor(
                                int(newv),
                                int(lseed * self.v[i - 1, j - k - 1]))
                    self.v[i - 1, j - 1] = newv

            # Multiply columns of V by appropriate power of 2.
            lseed = 1
            for j in range(self.maxcol - 1, 0, -1):
                lseed *= 2
                self.v[0:dim_num, j - 1] = self.v[0:dim_num, j - 1] * lseed

            # RECIPD is 1/(common denominator of the elements in V).
            self.recipd = 1.0 / (2 * lseed)
            self.lastq = np.zeros(dim_num)

        seed = int(np.floor(seed))

        if seed < 0:
            seed = 0

        lseed = 1
        if seed == 0:
            self.lastq = np.zeros(dim_num)

        elif seed == self.seed_save + 1:

            #  Find the position of the right-hand zero in SEED.
            lseed = self.i4_bit_lo0(seed)

        elif seed <= self.seed_save:

            self.seed_save = 0
            self.lastq = np.zeros(dim_num)

            for seed_temp in range(int(self.seed_save), int(seed)):
                lseed = self.i4_bit_lo0(seed_temp)
                for i in range(1, dim_num + 1):
                    self.lastq[i - 1] = np.bitwise_xor(
                        int(self.lastq[i - 1]), int(self.v[i - 1, lseed - 1]))

            lseed = self.i4_bit_lo0(seed)

        elif self.seed_save + 1 < seed:

            for seed_temp in range(int(self.seed_save + 1), int(seed)):
                lseed = self.i4_bit_lo0(seed_temp)
                for i in range(1, dim_num + 1):
                    self.lastq[i - 1] = np.bitwise_xor(
                        int(self.lastq[i - 1]), int(self.v[i - 1, lseed - 1]))

            lseed = self.i4_bit_lo0(seed)

        # Check that the user is not calling too many times!
        if self.maxcol < lseed:
            print('I4_SOBOL - Fatal error!')
            print('  Too many calls!')
            print('  MAXCOL = %d\n' % self.maxcol)
            print('  L =      %d\n' % lseed)
            return

        # Calculate the new components of QUASI.
        quasi = np.zeros(dim_num)
        for i in range(1, dim_num + 1):
            quasi[i - 1] = self.lastq[i - 1] * self.recipd
            self.lastq[i - 1] = np.bitwise_xor(
                int(self.lastq[i - 1]), int(self.v[i - 1, lseed - 1]))

        self.seed_save = seed
        seed += 1

        return [quasi, seed]

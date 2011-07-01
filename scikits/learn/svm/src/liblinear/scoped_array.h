#ifndef SCOPED_ARRAY_H_
#define SCOPED_ARRAY_H_

/**
 * RAII-style, exception-safe array wrapper.
 *
 * Lightweight alternative to std::vector, similar to boost::scoped_array.
 * Does not do bounds checking.
 *
 * Author: Lars Buitinck, ILPS, University of Amsterdam.
 */
template <typename T>
class scoped_array
{
    T *a;

  public:
    scoped_array(T *a) : a(a) {}
    ~scoped_array() { delete[] a; }

    T &operator[](size_t i) { return a[i]; }
    T const &operator[](size_t i) const { return a[i]; }

    T *get() const { return a; }

  private:
    scoped_array(scoped_array const &);      // non-copyable
    void operator=(scoped_array const &);    // non-assignable
};

#endif

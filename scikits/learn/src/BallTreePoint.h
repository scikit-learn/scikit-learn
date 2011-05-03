#ifndef BALLTREE_POINT_H
#define BALLTREE_POINT_H

#include <vector>

typedef std::vector<double> Point;

// Needed because cython doesn't handle operator[] as an lvalue
// http://osdir.com/ml/python.cython.devel/2008-03/msg00007.html
inline void SET(Point *p, int idx, double val) {
  (*p)[idx] = val;
}

#endif //BALLTREE_POINT_H

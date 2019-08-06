
/* this is a hack to generate libsvm with both sparse and dense
   methods in the same binary*/

#define _DENSE_REP
#include "svm.cpp"
#undef _DENSE_REP
#include "svm.cpp"

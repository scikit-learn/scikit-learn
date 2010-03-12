#ifndef BALLTREE_POINT_H
#define BALLTREE_POINT_H

#include <Python.h>
#include <numpy/arrayobject.h>
#include <exception>

/************************************************************
 * BallTree_Point
 *  This is a point class for use with the templated Ball
 *  Tree code, which is implemented with an internal reference
 *  to a PyArray object. 
 *
 * It can be initialized either with a PyObject pointer to a
 *  numpy array, or with 
 ************************************************************/

class BallTree_Point{
 public:
  typedef double value_type;

  BallTree_Point(npy_intp size)
    : base_arr_(PyArray_SimpleNew(1,&size,NPY_DOUBLE) ), 
      data_ptr_(0), inc_(1), size_(size)
  {
    data_ptr_ = (double*)PyArray_DATA(base_arr_);
  }

  BallTree_Point(PyObject* base_arr,
		 double* data_ptr,
		 int inc,
		 int size)
    : base_arr_(base_arr), data_ptr_(data_ptr),
      inc_(inc), size_(size){
    if(base_arr != NULL)
      Py_INCREF(base_arr);
  }

  BallTree_Point(PyObject* base_arr) 
    : base_arr_(base_arr), data_ptr_(0), inc_(0), size_(0)
  {
    if(base_arr_==0)
      return;
    Py_INCREF(base_arr_);
    int nd = PyArray_NDIM(base_arr_);
    
    if( PyArray_TYPE(base_arr_)!=NPY_DOUBLE ){
      Py_DECREF(base_arr_);
      throw std::exception();
    }
    
    if(nd==1)
      {
	size_ = PyArray_DIMS(base_arr_)[0];
	inc_ = PyArray_STRIDES(base_arr_)[0] / PyArray_DESCR(base_arr_)->elsize;
	data_ptr_ = (double*)PyArray_GETPTR1(base_arr_,0);
      }
    else
      {
	Py_DECREF(base_arr_);
	throw std::exception();
      }
  }
  
  BallTree_Point(const BallTree_Point& BTP)
    : base_arr_(BTP.base_arr_), data_ptr_(BTP.data_ptr_), 
      inc_(BTP.inc_), size_(BTP.size_)
  {
    Py_INCREF(base_arr_);
  }
  
  ~BallTree_Point(){
    if(base_arr_==NULL){
      if(size_>0)
	delete[] data_ptr_;
    }else{
      Py_DECREF(base_arr_);
    }
  }

  BallTree_Point& operator=(const BallTree_Point& BTP){
    //store in temporary variable in case we're
    // assigning self - data cannot accidentally be deleted
    PyObject* tmp = base_arr_;
    base_arr_ = BTP.base_arr_;
    if(base_arr_ != NULL)
      Py_INCREF(base_arr_);
    if(tmp != NULL)
      Py_DECREF(tmp);

    size_ = BTP.size_;
    inc_ = BTP.inc_;
    data_ptr_ = BTP.data_ptr_;

    return *this;
  }

  const int& size() const{return size_;}
  const int& inc() const{return inc_;}
  const double* arr() const{return data_ptr_;}
  double* arr(){return data_ptr_;}

  double& operator[](int i){return data_ptr_[inc_*i];}
  const double& operator[](int i) const{return data_ptr_[inc_*i];}
  double& at(int i){return data_ptr_[inc_*i];}
  const double& at(int i) const{return data_ptr_[inc_*i];}
 private:
  PyObject* base_arr_;
  double* data_ptr_;
  int inc_;
  int size_;
};

//Manhattan Distance
double P1_Dist(const BallTree_Point& p1,
	       const BallTree_Point& p2){
  int D = p1.size();
  if(p2.size() != D){
    std::cerr << "P1_Dist : point sizes must match\n";
    std::exit(-1);
  }
  double dist = 0;
  for(int i=0;i<D;i++){
    dist += fabs(p1[i] - p2[i]);
  }
  return dist;
}

//Euclidean Distance
double P2_Dist(const BallTree_Point& p1,
	       const BallTree_Point& p2){
  int D = p1.size();
  if(p2.size() != D){
    std::cerr << "P2_Dist : point sizes must match\n";
    std::exit(-1);
  }
  double dist = 0;
  double diff;
  for(int i=0;i<D;i++){
    diff = p1[i] - p2[i];
    dist += diff*diff;
  }
  return sqrt(dist);
}

//Infinity Distance
double Pinf_Dist(const BallTree_Point& p1,
	       const BallTree_Point& p2){
  int D = p1.size();
  if(p2.size() != D){
    std::cerr << "PN_Dist : point sizes must match\n";
    std::exit(-1);
  }
  double dist = 0;
  for(int i=0; i<D; i++){
    dist = std::max( dist, fabs(p1[i]-p2[i]) );
  }
  return dist;
}

#endif //BALLTREE_POINT_H

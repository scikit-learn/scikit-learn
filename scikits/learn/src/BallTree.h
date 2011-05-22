#ifndef BALL_TREE_H
#define BALL_TREE_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

/************************************************************
 * templated Ball Tree class
 *  written by Jake VanderPlas, October 2009
 *   vanderplas@astro.washington.edu
 *
 *  this will work with any user-defined Point type.
 *  Point type must define the following methods:
 *
 *  class Point{
 *     typedef double value_type;          //type of the data
 *                                         // (can be any numerical type)
 *
 *     Point::Point(size_t size);          //constructor of new point
 *
 *     Point::Point(const Point&);         //copy constructor
 *     Point& operator=(const Point&);     //assignment operator
 *                                         //  may either copy or view
 *                                         //   the other data
 *
 *     value_type& operator[](size_t i);   //element access
 *     value_type& at(size_t i);
 *     const value_type& operator[](size_t i) const;
 *     const value_type& at(size_t i) const;
 *
 *     size_t size() const;                //size (dimension) of the point
 *  }
 *
 *  Note that all these requirements are satisfied by the
 *   std::vector container
 ************************************************************/


/* Custom exception to allow Python to catch C++ exceptions */
class BallTreeException : public std::runtime_error {
  public:
    BallTreeException(std::string const &msg) : std::runtime_error(msg) {}
    ~BallTreeException() throw() {};
};


/************************************************************
 * Euclidean distance between two points
 *  this is used by default, but a function pointer to a
 *   user-defined distance function can be passed to BallTree
 *   below.
 ************************************************************/
template<class P1_Type,class P2_Type>
 typename P1_Type::value_type Euclidean_Dist(const P1_Type& p1,
                                             const P2_Type& p2){
    size_t D = p1.size();
    if(p2.size() != D){
        std::stringstream oss;
        oss << "Euclidean_Dist : point sizes must match (" << D << " != " << p2.size() << ").\n";
        throw BallTreeException(oss.str());
    }
    typename P1_Type::value_type dist = 0;
    typename P1_Type::value_type diff;
    for(size_t i=0;i<D;i++){
        diff = p1[i] - p2[i];
        dist += diff*diff;
    }
    return std::sqrt(dist);
}

template<class P1_Type,class P2_Type>
 typename P1_Type::value_type Euclidean_Dist(const P1_Type* p1,
                                             const P2_Type* p2){
  return Euclidean_Dist(*p1, *p2);
}


/************************************************************
 * VectorView object
 *  allows views of a portion of a std::vector object.
 *  without copying the data.
 *  this is used for keeping track of indices while
 *  constructing the Ball Tree
 ************************************************************/
template<typename T, class Alloc=std::allocator<T> >
class VectorView{
public:
    typedef typename std::vector<T,Alloc>::iterator iterator;

    VectorView(std::vector<T,Alloc>& data)
        : data_(data), istart_(0), size_(data.size()){}
    VectorView(std::vector<T,Alloc>& data, size_t istart, size_t size)
        : data_(data), istart_(istart), size_(size){}

    VectorView(VectorView<T,Alloc>& V, size_t istart, size_t size)
        : data_(V.data_), istart_(V.istart_+istart), size_(size) {}

    size_t size() const{ return size_; }
    const T& operator[](size_t i) const { return data_[istart_+i]; }
    T& operator[](size_t i) { return data_[istart_+i]; }

    iterator begin(){return data_.begin() + istart_;}
    iterator end(){return data_.begin() + istart_ + size_;}

private:
    std::vector<T,Alloc>& data_;
    size_t istart_;
    size_t size_;
};

/************************************************************
 * argsort:
 *  wrapper of std::sort which is used to quickly sort the
 *   points across a single dimension
 ************************************************************/
template <class Point>
class LT_Indices{
public:
    LT_Indices(std::vector<Point*> const &x, size_t elem=0) : points(x), e(elem) {}
    bool operator()(size_t i, size_t j) const
    {
        return points.at(i)->at(e) < points.at(j)->at(e);
    }
private:
    std::vector<Point*> const &points;
    size_t e;
};

template<class Point,class IndexVec>
void argsort(const std::vector<Point*>& x ,
	     IndexVec& indices,
	     size_t elem = 0)
{
    std::sort(indices.begin(), indices.end(), LT_Indices<Point>(x,elem));
}

/**********************************************************************
 *  Node
 *   This is an object representing a node of the Ball Tree, with
 *   the following data members:
 *
 *    Points   : a reference to the full data set
 *    sub{1,2} : subnodes to this node
 *    indices  : pointer to the BallTree index array.  Contains
 *               the list of indices in this node
 *    centroid : the centroid of all points within this node
 *    radius   : the maximum distance from the centroid to any
 *               point contained within this node
 *    dist_LB  : temporary variable which holds the lower-bound
 *               distance from a test point to any point in the node
 **********************************************************************/
template<class Point>
class Node{
public:
  //----------------------------------------------------------------------
  // Node : typedefs
  typedef typename Point::value_type value_type;
  typedef typename std::vector<value_type>::iterator value_iterator;
  typedef value_type (*DistFunc)(const Point&, const Point&);

private:
  //----------------------------------------------------------------------
  // Node : Data Members
  const std::vector<Point*>& Points;
  VectorView<size_t> indices;
  std::auto_ptr<Node<Point> > sub1, sub2;
  bool is_leaf;
  value_type radius;
  Point centroid;
  DistFunc Dist;

public:
  //----------------------------------------------------------------------
  // Node::Node
  //  this is called recursively to construct the ball tree
    Node(const std::vector<Point*>& Points,
        VectorView<size_t> indices,
        size_t leaf_size,
        DistFunc Dist_Func = &(Euclidean_Dist<Point,Point>),
        size_t depth=0)
        : Points(Points), indices(indices),
        is_leaf(false), radius(-1),
        centroid(Points[0]->size()), Dist(Dist_Func)
    {
        size_t D = Points[0]->size();
        size_t N = indices.size();
        if(N==0){
          throw BallTreeException("Node : zero-sized node\n   Abort\n");
        }else if (N==1){
            radius = 0;
            is_leaf = true;
            centroid = *(Points[indices[0]]);
        }else{
      //determine centroid
            for(size_t j=0;j<D;++j){
                centroid[j] = 0;
                for(size_t i=0;i<N;++i){
                    centroid[j] += Points[indices[i]]->at(j);
                }
                centroid[j] /= N;
            }

      //determine radius of node
            radius = 0;
            for(size_t i=0;i<N;i++){
                value_type this_R = Dist( *(Points[indices[i]]),centroid);
                if(this_R > radius)
                    radius = this_R;
            }

            if( N <= leaf_size ){
                is_leaf = true;
            }else{
    //sort the indices
                value_type range_max = 0;
                size_t imax = 0;

    //determine which dimension of M has the largest spread
                for(size_t i=0; i<D; ++i){
                    value_type rmin = Points[indices[0]]->at(i);
                    value_type rmax = Points[indices[0]]->at(i);
                    for(size_t j=1; j<N; ++j){
                        value_type Mij =  Points[indices[j]]->at(i);
                        if(Mij < rmin){
                            rmin = Mij;
                        }else if(Mij > rmax){
                            rmax = Mij;
                        }
                    }
                    if(rmax-rmin >= range_max){
                        range_max = rmax-rmin;
                        imax = i;
                    }
                }
    //sort along this dimension
                argsort(Points,indices,imax);

    //divide the points into two groups along this dimension
                sub1.reset(new Node<Point>(Points,
                                           VectorView<size_t>(indices,0,N/2),
                                           leaf_size, Dist, depth+1));
                sub2.reset(new Node<Point>(Points,
                                           VectorView<size_t>(indices,
                                                              N/2,N-(N/2)),
                                           leaf_size, Dist, depth+1));
            }
        }
    }

  //----------------------------------------------------------------------
  // Node::calc_dist_LB
  //   calculate the lower bound of distances from a test point to the
  //   points in this node
    value_type calc_dist_LB(const Point& P){
        return std::max( 0.0, Dist(P,centroid)-radius );
    }

  //----------------------------------------------------------------------
  // Node::calc_dist_UB
  //   calculate the upper bound of distances from a test point to the
  //   points in this node
    value_type calc_dist_UB(const Point& P){
        return Dist(P,centroid) + radius;
    }

  //----------------------------------------------------------------------
  // Node::query
  //  recursively called to query the tree
  //   for the k nearest neighbors of a point
  //
  //  Parameters:
  //    pt       : the point being queried
  //    k        : the number of nearest neighbors desired
  //    Pts_dist : distances to the points found so far (length k)
  //    Pts_ind  : indices of the points found so far   (length k)
  //    dist_LB  : the lower-bound distance.  This is passed as a parameter
  //                to keep from calculating it multiple times.
  void query(const Point& pt,
	     size_t k,
	     std::vector<value_type>& Pts_dist,
	     std::vector<size_t>& Pts_ind,
	     value_type dist_LB = -1){
    
    if(dist_LB < 0)
      dist_LB = calc_dist_LB(pt);
    
    
    //--------------------------------------------------
    //Case 1: all points in the node are further than
    //        the furthest point
    // We don't need to search this node or its subnodes
    if(dist_LB >= Pts_dist[k-1]){return;}
    
    //--------------------------------------------------
    //Case 2: this is a leaf node
    // we've found a leaf node containing points that are
    //  closer than those we've found so far.
    // Update the Pts_dist, Pts_ind and Dsofar
    else if(is_leaf)
      {
	value_iterator iter1;
	size_t position;
	value_type dist;
	//Do a brute-force search through the points.
	for(size_t i_pt=0; i_pt<indices.size(); i_pt++){
	  
	  //determine distance to this point
	  if( indices.size()==1 )
	    dist = dist_LB;
	  else
	    dist = Dist(pt, *(Points[ indices[i_pt] ]));
	  if(dist > Pts_dist[k-1])
	    continue;

	  //Determine where this point fits in the list
	  iter1 = std::lower_bound(Pts_dist.begin(),
				   Pts_dist.end(),
				   dist);
	  position = iter1 - Pts_dist.begin();

	  //Insert point into the list
	  for(size_t i=k-1; i>position; i--){
	    Pts_dist[i] = Pts_dist[i-1];
	    Pts_ind[i]  = Pts_ind[i-1];
	  }
	  Pts_dist[position] = dist;
	  Pts_ind[position]  = indices[i_pt];  
	}
      }
    
    
    //--------------------------------------------------
    //Case 3:
    //   Node is not a leaf, so we recursively query each
    //    of its subnodes
    else
      {
	value_type dist_LB_0 = sub1->calc_dist_LB(pt);
	value_type dist_LB_1 = sub2->calc_dist_LB(pt);
	if(dist_LB_0 <= dist_LB_1){
	  sub1->query(pt,k,Pts_dist,Pts_ind,dist_LB_0);
	  sub2->query(pt,k,Pts_dist,Pts_ind,dist_LB_1);
	}else{
	  sub2->query(pt,k,Pts_dist,Pts_ind,dist_LB_1);
	  sub1->query(pt,k,Pts_dist,Pts_ind,dist_LB_0);
	}
      }
  }

  //query all points within a radius r of pt
  size_t query_radius(const Point& pt,
                      value_type r,
                      std::vector<size_t>& nbrs){
    value_type dist_LB = calc_dist_LB(pt);

    //--------------------------------------------------
    //  Case 1: balls do not intersect.  return
    if( dist_LB > r ){}

    //--------------------------------------------------
    //  Case 2: this node is inside ball.  Add all pts to nbrs
    else if(dist_LB + 2*radius <= r){
      for(size_t i=0; i<indices.size(); i++)
	nbrs.push_back(indices[i]);
    }

    //--------------------------------------------------
    //  Case 3: node is a leaf.  Go through points and
    //          determine if they're in the ball
    else if(is_leaf){
      for(size_t i=0; i<indices.size(); i++){
	value_type D = Dist( pt,*(Points[indices[i]]) );
	if( D<=r )
	  nbrs.push_back(indices[i]);
      }
      return nbrs.size();
    }
    
    //--------------------------------------------------
    //  Case 4: node is not a leaf.  Recursively search
    //          subnodes
    else{
      sub1->query_radius(pt,r,nbrs);
      sub2->query_radius(pt,r,nbrs);
    }
    return nbrs.size();
  }
  
  //count all points within a radius r of p
  size_t query_radius(const Point& pt, value_type r)
  {
    value_type dist_LB = calc_dist_LB(pt);
    
    //--------------------------------------------------
    //  Case 1: balls do not intersect. return 0
    if( dist_LB > r )
      return 0;
    
    //--------------------------------------------------
    //  Case 2: this node is inside ball.  Add all pts to nbrs
    else if(dist_LB + 2*radius <= r){
      return indices.size();
    }
    
    //--------------------------------------------------
    //  Case 3: node is a leaf.  Go through points and
    //          determine if they're in the ball
    else if(is_leaf){
      size_t count = 0;
      for(size_t i=0; i<indices.size(); i++){
	value_type D = Dist( pt,*(Points[indices[i]]) );
	if( D<=r )
	  ++count;
      }
      return count;
    }
    
    //--------------------------------------------------
    //  Case 4: node is not a leaf.  Recursively search
    //          subnodes
    else{
      size_t count = 0;
      count += sub1->query_radius(pt,r);
      count += sub2->query_radius(pt,r);
      return count;
    }
  }
};



/**********************************************************************
 *  BallTree
 *   a class structure interface to a Ball Tree
 *   contains two data members:
 *    Points_    : a reference to the Matrix of input data
 *    head_node_ : a pointer to the root node of the Ball Tree
 *    indices_   : a vector which is referenced by all nodes in the tree
 **********************************************************************/
template<class Point>
class BallTree{
public:
    typedef typename Point::value_type value_type;
    typedef typename Node<Point>::DistFunc DistFunc;

  //----------------------------------------------------------------------
  // BallTree::BallTree
    BallTree(std::vector<Point*>* Points,
             size_t leaf_size = 1,
             DistFunc Dist_Func = &(Euclidean_Dist<Point,Point>))
        : Points_(*Points), indices_(Points_.size()),
    Dist(Dist_Func), leaf_size_(leaf_size){
    //initialize indices
        for(size_t i=0;i<indices_.size();i++)
            indices_[i] = i;

    //construct head node: this will recursively
    // construct the entire tree
        head_node_ = new Node<Point>(Points_,
            VectorView<size_t>(indices_),
            leaf_size,
            Dist);
    }

  //----------------------------------------------------------------------
  // BallTree::~BallTree
    ~BallTree(){
        delete head_node_; // this will recursively delete all subnodes
    }

    size_t PointSize() const{
        return Points_[0]->size();
    }

  //----------------------------------------------------------------------
  // BallTree::query
  //   on return, nbrs contains the unsorted indices of the k
  //     nearest neighbors.  The distance to the furthest neighbor
  //     is also returned.
    void query(const Point& pt,
	       std::vector<size_t>& nbrs) const{
      query(pt, nbrs.size(), &nbrs[0]);
    }
    
    void query(const Point* pt,
	       std::vector<size_t>& nbrs) const{
      query(*pt, nbrs.size(), &nbrs[0]);
    }
    void query(const Point& pt,
	       std::vector<size_t>& nbrs,
	       std::vector<value_type>& dist) const{
      query(pt, nbrs.size(), &nbrs[0], &dist[0]);
    }
    
    void query(const Point* pt,
	       std::vector<size_t>& nbrs,
	       std::vector<value_type>& dist) const{
      query(*pt, nbrs.size(), &nbrs[0], &dist[0]);
    }
    
    void query(const Point& pt,
	       size_t num_nbrs,
	       size_t* nbrs,
	       double* dist = 0) const{
      if(num_nbrs > Points_.size() ){
	std::stringstream oss;
	oss << "query: k must be less than or equal to N Points (" 
	    << num_nbrs << " > " << Points_.size() << ")\n";
	throw BallTreeException(oss.str());
      }

      std::vector<value_type> 
	Pts_dist(num_nbrs, std::numeric_limits<value_type>::max());
      
      std::vector<size_t> Pts_ind(num_nbrs, 0);
      
      head_node_->query(pt,num_nbrs,Pts_dist,Pts_ind);
      
      
      for(size_t i=0;i<num_nbrs;i++){
	nbrs[i] = Pts_ind[i];
      }
      
      if(dist != 0)
	for(size_t i=0;i<num_nbrs;i++)
	  dist[i] = Pts_dist[i];
      
    }
    
    //count number of points within a distance r of the point
    // return number of points.
    // on return, nbrs is an array of indices of nearest points
    // if nbrs is not supplied, just count points within radius
    size_t query_radius(const Point& pt, value_type r, std::vector<size_t>& nbrs)
    {
        return head_node_->query_radius(pt,r,nbrs);
    }

    size_t query_radius(const Point& pt, value_type r)
    {
        return head_node_->query_radius(pt,r);
    }
    
    size_t query_radius(const Point* pt, value_type r, std::vector<size_t>& nbrs)
    {
        return head_node_->query_radius(*pt,r,nbrs);
    }

    size_t query_radius(const Point* pt, value_type r)
    {
        return head_node_->query_radius(*pt,r);
    }
    
private:
  //----------------------------------------------------------------------
  // BallTree : Data Members
    Node<Point>* head_node_;
    const std::vector<Point*>& Points_;
    std::vector<size_t> indices_;
    DistFunc Dist;
    size_t leaf_size_;
};

#endif //BALL_TREE_H

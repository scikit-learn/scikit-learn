#ifndef BALL_TREE_H
#define BALL_TREE_H

#include <math.h>
#include <vector>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <exception>

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
 *     Point::Point(int size);             //constructor of new point
 *
 *     Point::Point(const Point&);         //copy constructor
 *     Point& operator=(const Point&);     //assignment operator
 *                                         //  may either copy or view
 *                                         //   the other data
 *
 *     value_type& operator[](int i);      //element access
 *     value_type& at(int i);
 *     const value_type& operator[](int i) const;
 *     const value_type& at(int i) const;
 *
 *     int size() const;                   //size (dimension) of the point
 *  }
 *
 *  Note that all these requirements are satisfied by the
 *   std::vector container
 ************************************************************/


/* Custom exception to allow Python to catch C++ exceptions */
class BallTreeException : public std::exception {
public:
    BallTreeException(std::string msg) { message = msg; }
    ~BallTreeException() throw() {};
    virtual const char* what() const throw() { return message.c_str(); }
private:
    std::string message;
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
    int D = p1.size();
    if(p2.size() != D){
        std::stringstream oss;
        oss << "Euclidean_Dist : point sizes must match (" << D << " != " << p2.size() << ").\n";
        throw BallTreeException(oss.str());
    }
    typename P1_Type::value_type dist = 0;
    typename P1_Type::value_type diff;
    for(int i=0;i<D;i++){
        diff = p1[i] - p2[i];
        dist += diff*diff;
    }
    return sqrt(dist);
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
    VectorView(std::vector<T,Alloc>& data,
        int istart,
        int size)
        : data_(data), istart_(istart), size_(size){}

    VectorView(VectorView<T,Alloc>& V,
        int istart,
        int size) : data_(V.data_), istart_(V.istart_+istart),
        size_(size){}

    int size() const{return size_;}
    const T& operator[](int i) const{return data_[istart_+i];}
    T& operator[](int i){return data_[istart_+i];}

    iterator begin(){return data_.begin() + istart_;}
    iterator end(){return data_.begin() + istart_ + size_;}

private:
    std::vector<T,Alloc>& data_;
    int istart_;
    int size_;
};

/************************************************************
 * argsort:
 *  wrapper of std::sort which is used to quickly sort the
 *   points across a single dimension
 ************************************************************/
template <class Point>
class LT_Indices{
public:
    LT_Indices(const std::vector<Point*>&x,
        int elem=0) : x_ptr(&x), e(elem){}
    bool operator()(int i1, int i2) const{
        return ( x_ptr->at(i1)->at(e) < x_ptr->at(i2)->at(e) );
    }
private:
    const std::vector<Point*>* x_ptr;
    int e;
};

template<class Point,class IndexVec>
void argsort(const std::vector<Point*>& x ,
    IndexVec& indices,
    int elem = 0)
{
    std::sort(indices.begin(), indices.end(), LT_Indices<Point>(x,elem));
}

/**********************************************************************
 * pd_tuple
 *  struct used to hold point indices and distances
 *   while querying the ball tree
 **********************************************************************/
template<typename value_type>
struct pd_tuple{
    pd_tuple(int i, value_type d) : index(i), dist(d){}
    int index;
    value_type dist;

    bool operator< (const pd_tuple<value_type>& RHS) const{
        return ( this->dist < RHS.dist );
    }
};

/**********************************************************************
 *  Node
 *   This is an object representing a node of the Ball Tree, with
 *   the following data members:
 *
 *    Points   : a reference to the full data set
 *    SubNodes : an array of subnodes to this node
 *    indices  : pointer to the BallTree index array.  Contains
 *               the list of indices in this node
 *    centroid : the centroid of all points within this node
 *    radius   : the maximum distance from the centroid to any
 *               point contained within this node
 *    dist_LB  : temporary variable which holds the lower-bound
 *               distance from a test point to any point in the node
 **********************************************************************/
template<class Point>
struct Node{
  //----------------------------------------------------------------------
  // Node : typedefs
    typedef typename Point::value_type value_type;
    typedef value_type (*DistFunc)(const Point&, const Point&);

  //----------------------------------------------------------------------
  // Node : Data Members
    const std::vector<Point*>& Points;
    VectorView<int> indices;
    std::vector<Node<Point>*> SubNodes;
    bool is_leaf;
    value_type radius;
    Point centroid;
    DistFunc Dist;

  //----------------------------------------------------------------------
  // Node::Node
  //  this is called recursively to construct the ball tree
    Node(const std::vector<Point*>& Points,
        VectorView<int> indices,
        int leaf_size,
        DistFunc Dist_Func = &(Euclidean_Dist<Point,Point>),
        int depth=0)
        : Points(Points), indices(indices),
        SubNodes(0), is_leaf(false), radius(-1),
        centroid(Points[0]->size()), Dist(Dist_Func)
    {
        int D = Points[0]->size();
        int N = indices.size();
        if(N==0){
          throw BallTreeException("Node : zero-sized node\n   Abort\n");
        }else if (N==1){
            radius = 0;
            is_leaf = true;
            centroid = *(Points[indices[0]]);
        }else{
      //determine centroid
            for(int j=0;j<D;++j){
                centroid[j] = 0;
                for(int i=0;i<N;++i){
                    centroid[j] += Points[indices[i]]->at(j);
                }
                centroid[j] /= N;
            }

      //determine radius of node
            radius = 0;
            for(int i=0;i<N;i++){
                value_type this_R = Dist( *(Points[indices[i]]),centroid);
                if(this_R > radius)
                    radius = this_R;
            }

            if( N <= leaf_size ){
                is_leaf = true;
            }else{
    //sort the indices
                value_type range_max = 0;
                int imax = -1;

    //determine which dimension of M has the largest spread
                for(int i=0; i<D; ++i){
                    value_type rmin = Points[indices[0]]->at(i);
                    value_type rmax = Points[indices[0]]->at(i);
                    for(int j=1; j<N; ++j){
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
                SubNodes.resize(2);
                SubNodes[0] = new Node<Point> (Points,
                                        VectorView<int>(indices,0,N/2),
                                        leaf_size, Dist, depth+1);
                SubNodes[1] = new Node<Point> (Points,
                                        VectorView<int>(indices,N/2,N-(N/2)),
                                        leaf_size, Dist, depth+1);
            }
        }
    }

  //----------------------------------------------------------------------
  // Node::~Node
    ~Node(){
    //recursively delete all lower nodes
        for(size_t i=0; i<SubNodes.size(); ++i)
            delete SubNodes[i];
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
  //    PointSet : the points found so far
  //    Dsofar   : the distance within which new points are of interest
  //                that is, if k points have already been found, then
  //                Dsofar is the distance to the furthest point.
  //    dist_LB  : the lower-bound distance.  This is passed as a parameter
  //                to keep from calculating it multiple times.
    value_type query(const Point& pt,
                int k,
                std::vector< pd_tuple<value_type> >& PointSet,
                value_type Dsofar = std::numeric_limits<value_type>::max(),
                value_type dist_LB = -1){
        if(dist_LB < 0){
            dist_LB = calc_dist_LB(pt);
        }

    //--------------------------------------------------
    //Case 1: all points in the node are further than
    //        the furthest point in PointSet
    // We don't need to search this node or its subnodes
        if(dist_LB >= Dsofar)
        {
            return Dsofar;
        }


    //--------------------------------------------------
    //Case 2: this is a leaf node
    // we've found a leaf node containing points that are
    //  closer than those we've found so far.
    // Update the PointSet and Dsofar
        else if(is_leaf)
        {
    //Do a brute-force search through the points.
            for(int i_pt=0; i_pt<indices.size(); i_pt++){

      //determine distance to this point
                value_type dist;
                if( indices.size()==1 )
                    dist = dist_LB;
                else
                    dist = Dist(pt, *(Points[ indices[i_pt] ]));
                if(dist > Dsofar)
                    continue;

      //If there are less than k points in PointSet,
      //  add this point to PointSet
                if(PointSet.size() < (size_t)k){
                    PointSet.push_back( pd_tuple<value_type>(indices[i_pt],
                                                             dist) );
        //determine new Dsofar, if necessary
                    if(PointSet.size() == (size_t)k){
                        std::sort( PointSet.begin(),PointSet.end() );
                        Dsofar = PointSet[k-1].dist;
                    }
                }
      //If there are already k points in PointSet, we
      //  must replace the most distant with this
                else{
                    PointSet[k-1] = pd_tuple<value_type>(indices[i_pt],dist);
                    std::sort( PointSet.begin(),PointSet.end() );
                    Dsofar = PointSet[k-1].dist;
                }
            }
            return Dsofar;
        }


    //--------------------------------------------------
    //Case 3:
    //   Node is not a leaf, so we recursively query each
    //    of its subnodes
        else
        {
            value_type dist_LB_0 = SubNodes[0]->calc_dist_LB(pt);
            value_type dist_LB_1 = SubNodes[1]->calc_dist_LB(pt);
            if(dist_LB_0 <= dist_LB_1){
                Dsofar = SubNodes[0]->query(pt,k,PointSet,Dsofar,dist_LB_0);
                Dsofar = SubNodes[1]->query(pt,k,PointSet,Dsofar,dist_LB_1);
            }else{
                Dsofar = SubNodes[1]->query(pt,k,PointSet,Dsofar,dist_LB_1);
                Dsofar = SubNodes[0]->query(pt,k,PointSet,Dsofar,dist_LB_0);
            }
            return Dsofar;
        }
    }

  //query all points within a radius r of pt
    int query_ball(const Point& pt,
        value_type r,
        std::vector<long int>& nbrs){
        value_type dist_LB = calc_dist_LB(pt);

    //--------------------------------------------------
    //  Case 1: balls do not intersect.  return
        if( dist_LB > r ){}

    //--------------------------------------------------
    //  Case 2: this node is inside ball.  Add all pts to nbrs
        else if(dist_LB + 2*radius <= r){
            for(int i=0; i<indices.size(); i++)
                nbrs.push_back(indices[i]);
        }

    //--------------------------------------------------
    //  Case 3: node is a leaf.  Go through points and
    //          determine if they're in the ball
        else if(is_leaf){
            for(int i=0; i<indices.size(); i++){
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
            SubNodes[0]->query_ball(pt,r,nbrs);
            SubNodes[1]->query_ball(pt,r,nbrs);
        }
        return nbrs.size();
    }

  //count all points within a radius r of p
    int query_ball(const Point& pt,
    value_type r){
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
            int count = 0;
            for(int i=0; i<indices.size(); i++){
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
            int count = 0;
            count += SubNodes[0]->query_ball(pt,r);
            count += SubNodes[1]->query_ball(pt,r);
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
        int leaf_size = 1,
        DistFunc Dist_Func = &(Euclidean_Dist<Point,Point>))
        : Points_(*Points), indices_(Points_.size()),
    Dist(Dist_Func), leaf_size_(leaf_size){
    //initialize indices
        for(size_t i=0;i<indices_.size();i++)
            indices_[i] = i;

    //construct head node: this will recursively
    // construct the entire tree
        head_node_ = new Node<Point>(Points_,
            VectorView<int>(indices_),
            leaf_size,
            Dist);
    }

  //----------------------------------------------------------------------
  // BallTree::~BallTree
    ~BallTree(){
        delete head_node_; // this will recursively delete all subnodes
    }

    int PointSize() const{
        return Points_[0]->size();
    }

  //----------------------------------------------------------------------
  // BallTree::query
  //   on return, nbrs contains the unsorted indices of the k
  //     nearest neighbors.  The distance to the furthest neighbor
  //     is also returned.
    value_type query(const Point& pt,
                     std::vector<long int>& nbrs) const{
        query(pt,nbrs.size(),&(nbrs[0]));
    }

    value_type query(const Point* pt,
                     std::vector<long int>& nbrs) const{
        query(*pt,nbrs.size(),&(nbrs[0]));
    }

    value_type query(const Point& pt,
                     int num_nbrs,
                     long int* nbrs,
                     double* dist = NULL) const{
        if(num_nbrs > (int)(Points_.size()) ){
            std::stringstream oss;
            oss << "query: k must be less than or equal to N Points (" << num_nbrs << " > " << (int)(Points_.size()) << ")\n";
            throw BallTreeException(oss.str());
        }

        std::vector<pd_tuple<value_type> > PointSet;
        value_type Dmax = head_node_->query(pt,num_nbrs,PointSet);


        for(int i=0;i<num_nbrs;i++){
            nbrs[i] = PointSet[i].index;
        }

        if(dist != NULL)
            for(int i=0;i<num_nbrs;i++)
            dist[i] = PointSet[i].dist;
        return Dmax;
    }

  //count number of points within a distance r of the point
  // return number of points.
  // on return, nbrs is an array of indices of nearest points
  // if nbrs is not supplied, just count points within radius
    int query_ball(const Point& pt,
                   value_type r,
                   std::vector<long int>& nbrs){
        return head_node_->query_ball(pt,r,nbrs);
    }

    int query_ball(const Point& pt,
    value_type r){
        return head_node_->query_ball(pt,r);
    }
private:
  //----------------------------------------------------------------------
  // BallTree : Data Members
    Node<Point>* head_node_;
    const std::vector<Point*>& Points_;
    std::vector<int> indices_;
    DistFunc Dist;
    int leaf_size_;
};

/************************************************************
 * BruteForceNeighbors
 *  determine neighbors using a brute-force algorithm.  This
 *  is for comparison with results from BallTree
 ************************************************************/

template<class Point>
void BruteForceNeighbors(std::vector<Point*> Points,
             const Point& pt,
             long int k,
             long int* neighbors){
  int N = Points.size();
  typedef pd_tuple<typename Point::value_type> PD;
  std::vector<PD> distances;

  for(int i=0;i<N;i++)
    distances.push_back( PD(i,Euclidean_Dist(pt,*(Points[i])) ) );

  std::partial_sort( distances.begin(),
             distances.begin()+k,
             distances.end()      );

  for(int i=0; i<k; i++)
    neighbors[i] = distances[i].index;
}

template<class Point>
void BruteForceNeighbors(std::vector<Point*> Points,
             const Point& pt,
             std::vector<long int>& neighbors){
  long int k = neighbors.size();
  BruteForceNeighbors(Points,pt,k,&(neighbors[0]));
}


#endif //BALL_TREE_H

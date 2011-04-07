#ifndef NODE_H
#define NODE_H

#include <vector>
#include <utility>
#include "Histogram.h"

using namespace std;

enum label_t
{
    SIGNAL,
    BACKGROUND
};

struct Object
{
    double* attrs;
    double weight;
    unsigned int dim;
    label_t label;
};

class Node
{
    public:
        
        Node():
            left_child(0),
            right_child(0),
            attribute(-1),
            cut(0),
            response(0)
        {}

        ~Node()
        {
            delete this->left_child;
            delete this->right_child;
        }

        void calc_response();

        //void update_classification();

        void add_signal(Object* s) { this->signal.push_back(s); }

        void add_background(Object* s) { this->background.push_back(s); }

        void recursive_split(unsigned int minleafsize, unsigned int nbins, unsigned int maxdepth = 0, unsigned int depth = 0);
        
        bool split(unsigned int minleafsize, unsigned int nbins);

        const Node* get_left_child() const
        {
            return left_child;
        }

        const Node* get_right_child() const
        {
            return right_child;
        }

        void set_left_child(Node* child)
        {
            if (this->left_child) delete this->left_child;
            this->left_child = child;
        }

        void set_right_child(Node* child)
        {
            if (this->right_child) delete this->right_child;
            this->right_child = child;
        }

        bool leaf() const
        {
            return !this->left_child && !this->right_child;
        }

        bool complete() const
        {
            return this->left_child && this->right_child;
        }

        double predict(const double* attrs) const;

        pair<double,double> minmax(unsigned int attribute) const;
    
    private:

        vector<Object*> signal;
        vector<Object*> background;
        Node* left_child;
        Node* right_child;
        int attribute;
        double cut;
        double response;
};

#endif

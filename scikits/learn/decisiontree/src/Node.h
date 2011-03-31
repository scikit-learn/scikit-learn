#ifndef NODE_H
#define NODE_H

#include <vector>
#include <utility>
#include "Object.h"
#include "Histogram.h"

using namespace std;

class Node
{
    public:
        
        Node():
            left_child(0),
            right_child(0),
            attribute(-1),
            cut(0),
            purity(0)
        {}

        ~Node(){ delete left_child; delete right_child; }

        void calc_purity();

        void update_classification();

        void add_signal(Object* s) { signal.push_back(s); }

        void add_background(Object* s) { background.push_back(s); }

        void set_signal(vector<Object*> s) { signal = s; }

        void set_background(vector<Object*> s) { background = s; }

        void recursive_split(unsigned int minLeafSize, unsigned int numCuts);
        
        bool split(unsigned int minSize, unsigned int resolution);

        const Node* get_left_child() const
        {
            return left_child;
        }

        const Node* get_right_child() const
        {
            return right_child;
        }

        void set_left_child(const Node* child)
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

        float predict(const double* attrs) const;

        pair<float,float> get_extrema(unsigned int attribute);

        void drop_objects()
        {
            this->signal.clear();
            this->background.clear();
        }

    private:

        vector<Object*> signal;
        vector<Object*> background;
        Node* left_child;
        Node* right_child;
        int attribute;
        double cut;
        double purity;
};

#endif

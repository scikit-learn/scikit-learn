#ifndef NODE_H
#define NODE_H

#include <vector>
#include <utility>
#include "Object.h"
#include "Histogram.h"

using namespace std;

enum Classes
{
    SIGNAL,
    BACKGROUND
};

class Node
{
    public:
        
        Node():
            leftChild(0),
            rightChild(0),
            splitAttribute(-1),
            splitAttributeValue(0),
            purity(0)
        {}

        ~Node(){ delete leftChild; delete rightChild; }

        float calc_purity();

        void update_classification();

        void add_signal(Object* s) { signal.push_back(s); }

        void add_background(Object* s) { background.push_back(s); }

        void set_signal(vector<Object*> s) { signal = s; }

        void set_background(vector<Object*> s) { background = s; }

        void recursive_split(unsigned int minLeafSize, unsigned int numCuts);
        
        bool split(unsigned int minSize, unsigned int resolution);

        Node* get_left_child() { return leftChild; }

        Node* get_right_child() { return rightChild; }

        void set_left_child(Node* child)
        {
            if (leftChild) delete leftChild;
            leftChild = child;
        }

        void set_right_child(Node* child)
        {
            if (rightChild) delete rightChild;
            rightChild = child;
        }

        bool leaf()
        {
            return !leftChild && !rightChild;
        }

        bool complete()
        {
            return leftChild&&rightChild;
        }

        float response(const Object* object);

        pair<float,float> get_extrema(unsigned int attribute);

    private:

        vector<Object*> signal;
        vector<Object*> background;
        Node* leftChild;
        Node* rightChild;
        int splitAttribute;
        float splitAttributeValue;
        float purity;
};

#endif

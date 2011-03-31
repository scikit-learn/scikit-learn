#ifndef BDT_H
#define BDT_H

#include "Node.h"
#include "Object.h"

using namespace std;

class BDT
{
    public:
        
        BDT(
        unsigned int numTrees = 10,
        float adaBoostBeta = .5,
        unsigned int maxDepth = 10,
        unsigned int minLeafSize = 100,
        unsigned int numCuts = 10):
            numTrees(numTrees),
            beta(adaBoostBeta),
            maxDepth(maxDepth),
            minLeafSize(minLeafSize),
            numCuts(numCuts)
        {}

        ~BDT()
        {
            vector<pair<float,Node*> >::iterator it(trees.begin());
            for (; it != trees.end(); ++it)
            {
                delete it->second;
            }
        }

        void fit(vector<Object*>& signal, vector<Object*>& background);

        float predict_single(const Object* object);

        void predict(const vector<Object*>&, vector<float>& scores);

    private:

        void adaboost(vector<Object*>& signal, vector<Object*>& background,float&,float&);

        unsigned int numTrees;
        float beta;
        unsigned int maxDepth;
        unsigned int minLeafSize;
        unsigned int numCuts;
        vector<pair<float,Node*> > trees;
};

#endif

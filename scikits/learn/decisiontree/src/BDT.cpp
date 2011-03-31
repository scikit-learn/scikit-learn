#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <math.h>

#include "BDT.h"

using namespace std;

void BDT::adaboost(vector<Object*>& signal, vector<Object*>& background, float& boostWeight, float& errorFraction) {

    float totalWeighted = 0;
    float errWeighted = 0;
    vector<Object*>::const_iterator it(signal.begin());
    for (; it != signal.end(); ++it)
    {
        totalWeighted += (*it)->get_weight();
        if ((*it)->get_class() == BACKGROUND) errWeighted += (*it)->get_weight();
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        totalWeighted += (*it)->get_weight();
        if ((*it)->get_class() == SIGNAL) errWeighted += (*it)->get_weight();
    }
    errorFraction = totalWeighted > 0 ? errWeighted / totalWeighted : 0.;
    float alpha;
    if (errorFraction > 0 && errorFraction < .5) alpha = beta*log((1.-errorFraction)/errorFraction);
    else alpha = 0.;
    float boostU = exp(alpha);
    float boostD = exp(-alpha);
    boostWeight = alpha;
    
    it = signal.begin();
    for (; it != signal.end(); ++it)
    {
        if ((*it)->get_class() == BACKGROUND) (*it)->boost_weight(boostU);
        else (*it)->boost_weight(boostD);
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        if ((*it)->get_class() == SIGNAL) (*it)->boost_weight(boostU);
        else (*it)->boost_weight(boostD);
    }
}

void BDT::fit(vector<Object*>& signal, vector<Object*>& background)
{
    float boostWeight;
    float errorFraction;
    Node* node;
    for (int i(0); i < numTrees; ++i)
    {
        node = new Node();
        node->set_signal(signal);
        node->set_background(background);
        node->recursive_split(minLeafSize,numCuts);
        adaboost(signal,background,boostWeight,errorFraction);
        cout << "AdaBoost alpha: " << boostWeight << endl;
        cout << "Error fraction: " << errorFraction << endl;
        cout << "Boost: " << i << endl;
        trees.push_back(pair<float,Node*>(boostWeight,node));
        if (boostWeight == 0 || errorFraction == 0) break;
    }
}

float BDT::predict_single(const Object* object)
{
    if (!object) return -1.;
    float norm(0);
    float score(0);
    vector<pair<float,Node*> >::iterator it(trees.begin());
    for (; it != trees.end(); ++it)
    {
           norm += it->first;
           score += it->first * it->second->response(object);
    }
    return norm > 0 ? score / norm : 0;
}

void BDT::predict(const vector<Object*>& objects, vector<float>& scores)
{
    scores.assign(objects.size(),0);
    vector<Object*>::const_iterator objects_it(objects.begin());
    vector<float>::iterator scores_it(scores.begin());
    for (; objects_it != objects.end(); ++objects_it)
    {
        *scores_it = predict_single(*objects_it);
    }
}

int main(void)
{
	srand ( time(NULL) );

    vector<Object*> signal;
    vector<Object*> background;
    Object* obj;
    int numFeatures = 5;
    vector<float> features(numFeatures,0);
    for (int i(0); i< 100; ++i)
    {
        obj = new Object();
        for (int j(0); j < numFeatures; ++j) features[j] = float(rand()%50)+10;
        obj->set_attributes(features);
        signal.push_back(obj);
        obj = new Object();
        for (int j(0); j < numFeatures; ++j) features[j] = float(rand()%50);
        obj->set_attributes(features);
        background.push_back(obj);
    }

    BDT* bdt = new BDT(50,.5,10,10,100);

    bdt->fit(signal,background);

    vector<Object*>::const_iterator it(signal.begin());
    /*
    for (; it != signal.end(); ++it)
    {
        cout << bdt->predict_single(*it) << endl;
    }*/

    it = signal.begin();
    for (; it != signal.end(); ++it) delete *it;
    it = background.begin();
    for (; it != background.end(); ++it) delete *it;
	delete bdt;
    return EXIT_SUCCESS;
}

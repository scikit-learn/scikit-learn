#include "Node.h"
#include <iostream>

using namespace std;

double Node::predict(const double* attrs) const
{
    if (this->leaf())
    {
        return this->response;
    }
    else if (attrs[this->attribute] > this->cut)
    {
        return this->right_child->predict(attrs);
    }
    else
    {
        return this->left_child->predict(attrs);
    }
}

void Node::recursive_split(unsigned int minleafsize, unsigned int nbins, unsigned int maxdepth, unsigned int depth)
{
    if ( (depth < maxdepth || maxdepth == 0) && this->split(minleafsize, nbins))
    {
        this->left_child->recursive_split(minleafsize, nbins, maxdepth, depth+1);
        this->right_child->recursive_split(minleafsize, nbins, maxdepth, depth+1);
    }
    else
    {
        this->calc_response();
    }
    this->signal.clear();
    this->background.clear();
}

void Node::calc_response()
{
    vector<Object*>::const_iterator it(this->signal.begin());
    double weightedSignal(0.);
    for (; it != this->signal.end(); ++it)
    {
        weightedSignal += (*it)->weight;
    }
    it = this->background.begin();
    double weightedBackground(0.);
    for (; it != this->background.end(); ++it)
    {
        weightedBackground += (*it)->weight;
    }
    this->response = (weightedBackground + weightedSignal) > 0 ? \
         -1. + 2. * weightedSignal / (weightedBackground + weightedSignal) : -1.;
}

/*
void Node::update_classification()
{
    vector<Object*>::const_iterator it(this->signal.begin());
    if (this->response > (this->sig_score - this->bkg_score)/2.)
    {
        for (; it != this->signal.end(); ++it)
            (*it)->label = SIGNAL;
        it = this->background.begin();
        for (; it != this->background.end(); ++it)
            (*it)->label = SIGNAL;
    }
    else
    {
        for (; it != this->signal.end(); ++it)
            (*it)->label = BACKGROUND;
        it = this->background.begin();
        for (; it != this->background.end(); ++it)
            (*it)->label = BACKGROUND;
    }
}
*/

bool Node::split(unsigned int minleafsize, unsigned int nbins)
{
    if (signal.size()==0 || background.size()==0 || nbins < 2 || signal.size() + background.size() <= minleafsize)
    {
        return false;
    }
    Object* sample = signal.at(0);
    unsigned int numAttributes = sample->dim;
    int bestAttribute(-1);
    double bestSplit(0.);
    pair<double,double> extrema;
    Histogram<double,double>* sigHist;
    Histogram<double,double>* bkgHist;
    double bestGini(-1.);
    double step;
    double min;
    double max;
    double Gini_left;
    double Gini_right;
    double sig_left, sig_right;
    double bkg_left, bkg_right;
    double purity_left, purity_right;

    for (unsigned int i(0); i < numAttributes; ++i)
    {
        // get max and min value of this attribute over signal and background combined
        extrema = this->minmax(i);
        
        // create histograms for signal and background
        step = (extrema.second - extrema.first)/nbins;
        if (step==0) continue;
        min = extrema.first-step/2;
        max = extrema.second+step/2;
        sigHist = new Histogram<double,double>(nbins, min, max);
        bkgHist = new Histogram<double,double>(nbins, min, max);
        step = (max - min)/nbins;
        
        // fill histograms
        std::vector<Object*>::const_iterator it(signal.begin());
        for (; it != signal.end(); ++it)
        {
            sigHist->fill((*it)->attrs[i],(*it)->weight);
        }
        
        it = background.begin();
        for (; it != background.end(); ++it)
        {
            bkgHist->fill((*it)->attrs[i],(*it)->weight);
        }

        // calculate Gini_left + Gini_right for a split at each internal bin boundary and find minimum where minleafsize is respected
        for (unsigned int binCut(1); binCut < nbins; ++binCut)
        {
            // This can be optimized
            if (sigHist->integral(0,binCut,false) + bkgHist->integral(0,binCut,false) >= minleafsize && \
                sigHist->integral(binCut,nbins,false) + bkgHist->integral(binCut,nbins,false) >= minleafsize)
            {
                sig_left = sigHist->integral(0,binCut);
                bkg_left = bkgHist->integral(0,binCut);
                sig_right = sigHist->integral(binCut,nbins);
                bkg_right = bkgHist->integral(binCut,nbins);
                purity_left = sig_left / (sig_left + bkg_left);
                purity_right = sig_right / (sig_right + bkg_right);
                Gini_left = purity_left * (1 - purity_left) * (sig_left + bkg_left);
                Gini_right = purity_right * (1 - purity_right) * (sig_right + bkg_right);
                
                // if a possible split is found and if this split is the best so far update bestAttribute and bestSplit
                if (Gini_left + Gini_right < bestGini || bestGini == -1)
                {
                    bestGini = Gini_left + Gini_right;
                    bestAttribute = i;
                    bestSplit = min + binCut * step;
                }
            }
        }
        delete sigHist;
        delete bkgHist;
    }

    cout << bestGini << endl;
    
    if (bestAttribute == -1) return false;
    
    // create left and right Nodes and add them as children
    Node* left = new Node();
    Node* right = new Node();
    
    std::vector<Object*>::const_iterator it(signal.begin());
    for (; it != signal.end(); ++it)
    {
        if ((*it)->attrs[bestAttribute] > bestSplit) right->add_signal(*it);
        else left->add_signal(*it);
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        if ((*it)->attrs[bestAttribute] > bestSplit) right->add_background(*it);
        else left->add_background(*it);
    }
    
    this->attribute = bestAttribute;
    this->cut = bestSplit;
    cout << "best attr: " << bestAttribute << endl;
    cout << "best cut:  " << bestSplit << endl;
    this->set_left_child(left);
    this->set_right_child(right);
    return true;
}

pair<double,double> Node::minmax(unsigned int attribute) const
{
    vector<Object*>::const_iterator it(signal.begin());
    double value = (*it++)->attrs[attribute];
    pair<double,double> extrema(value,value); // min,max
    for (; it != signal.end(); ++it)
    {
        value = (*it)->attrs[attribute];
        if (value < extrema.first) extrema.first = value;
        if (value > extrema.second) extrema.second = value;
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        value = (*it)->attrs[attribute];
        if (value < extrema.first) extrema.first = value;
        if (value > extrema.second) extrema.second = value;
    }
    return extrema;
}

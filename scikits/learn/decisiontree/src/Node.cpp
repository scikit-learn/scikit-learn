#include "Node.h"

#include <iostream>

float Node::response(const Object* object)
{
    if (leaf())
    {
        return purity;
    }
    else if (!complete())
    {
        return -100.;
    }
    else if (object->attribute(splitAttribute) > splitAttributeValue)
    {
        return get_right_child()->response(object);
    }
    else
    {
        return get_left_child()->response(object);
    }
}


void Node::recursive_split(unsigned int minLeafSize, unsigned int numCuts)
{
    cout << calc_purity() << endl;
    if (split(minLeafSize,numCuts))
    {
        get_left_child()->recursive_split(minLeafSize,numCuts);
        get_right_child()->recursive_split(minLeafSize,numCuts);
    }
    else
    {
        purity = calc_purity();
        update_classification();
    }
}

float Node::calc_purity()
{
    vector<Object*>::const_iterator it(signal.begin());
    float weightedSignal = 0.;
    for (; it != signal.end(); ++it)
        weightedSignal += (*it)->get_weight();
    it = background.begin();
    float weightedBackground = 0.;
    for (; it != background.end(); ++it)
        weightedBackground += (*it)->get_weight();
    return (weightedSignal + weightedBackground) > 0 ? weightedSignal / (weightedSignal + weightedBackground) : -1.;
}

void Node::update_classification()
{
    vector<Object*>::const_iterator it(signal.begin());
    if (purity > .5)
    {
        for (; it != signal.end(); ++it)
            (*it)->set_class(SIGNAL);
        it = background.begin();
        for (; it != background.end(); ++it)
            (*it)->set_class(SIGNAL);
    }
    else
    {
        for (; it != signal.end(); ++it)
            (*it)->set_class(BACKGROUND);
        it = background.begin();
        for (; it != background.end(); ++it)
            (*it)->set_class(BACKGROUND);
    }
}

bool Node::split(unsigned int minSize, unsigned int resolution)
{
    if (signal.size()==0 || background.size()==0 || resolution < 2 || signal.size() + background.size() <= minSize)
    {
        return false;
    }
    Object* sample = signal.at(0);
    unsigned int numAttributes = sample->num_attributes();
    int bestAttribute = -1;
    float bestSplit = 0.;
    pair<float,float> extrema;
    Histogram<double,float>* sigHist;
    Histogram<double,float>* bkgHist;
    float bestGini = -1;
    for (unsigned int i(0); i < numAttributes; ++i)
    {
        // get max and min value of this attribute over signal and background combined
        extrema = get_extrema(i);
        // create histograms for signal and background
        float step = (extrema.second - extrema.first)/resolution;
        if (step==0) continue;
        sigHist = new Histogram<double,float>(resolution,extrema.first,extrema.second+step);
        bkgHist = new Histogram<double,float>(resolution,extrema.first,extrema.second+step);
        // fill histograms
        std::vector<Object*>::const_iterator it(signal.begin());
        for (; it != signal.end(); ++it)
            sigHist->fill((*it)->attribute(i),(*it)->get_weight());
        it = background.begin();
        for (; it != background.end(); ++it)
            bkgHist->fill((*it)->attribute(i),(*it)->get_weight());
        // calculate Gini_left + Gini_right for a split at each internal bin boundary and find minimum where minSize is respected
        float Gini_left;
        float Gini_right;
        float sig_left, sig_right;
        float bkg_left, bkg_right;
        float purity_left, purity_right;
        for (unsigned int binCut(1); binCut < resolution; ++binCut)
        {
            if (sigHist->integral(0,binCut,false) + bkgHist->integral(0,binCut,false) >= minSize && \
                    sigHist->integral(binCut,resolution,false) + bkgHist->integral(binCut,resolution,false) >= minSize)
            {
                sig_left = sigHist->integral(0,binCut);
                bkg_left = bkgHist->integral(0,binCut);
                sig_right = sigHist->integral(binCut,resolution);
                bkg_right = bkgHist->integral(binCut,resolution);
                purity_left = sig_left / (sig_left + bkg_left);
                purity_right = sig_right / (sig_right + bkg_right);
                Gini_left = purity_left * (1 - purity_left) * (sig_left + bkg_left);
                Gini_right = purity_right * (1 - purity_right) * (sig_right + bkg_right);
                // if a possible split is found and if this split is the best so far update bestAttribute and bestSplit
                if (Gini_left + Gini_right < bestGini || bestGini == -1)
                {
                    bestGini = Gini_left + Gini_right;
                    bestAttribute = i;
                    bestSplit = extrema.first + binCut * (extrema.second + step - extrema.first) / resolution;
                }
            }
        }
        delete sigHist;
        delete bkgHist;
    }
    if (bestAttribute == -1) return false;
    // create left and right partitions and add them as children
    Node* left = new Node();
    Node* right = new Node();
    std::vector<Object*>::const_iterator it(signal.begin());
    for (; it != signal.end(); ++it)
    {
        if ((*it)->attribute(bestAttribute) > bestSplit) right->add_signal(*it);
        else left->add_signal(*it);
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        if ((*it)->attribute(bestAttribute) > bestSplit) right->add_background(*it);
        else left->add_background(*it);
    }
    splitAttribute = bestAttribute;
    set_left_child(left);
    set_right_child(right);
    return true;
}

pair<float,float> Node::minmax(unsigned int attribute)
{
    vector<Object*>::const_iterator it(signal.begin());
    float value = (*it++)->attribute(attribute);
    pair<float,float> extrema(value,value); // min,max
    for (; it != signal.end(); ++it)
    {
        value = (*it)->attribute(attribute);
        if (value < extrema.first) extrema.first = value;
        if (value > extrema.second) extrema.second = value;
    }
    it = background.begin();
    for (; it != background.end(); ++it)
    {
        value = (*it)->attribute(attribute);
        if (value < extrema.first) extrema.first = value;
        if (value > extrema.second) extrema.second = value;
    }
    return extrema;
}

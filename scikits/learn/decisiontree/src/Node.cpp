#include "Node.h"

double Node::predict(const double* attrs)
{
    if (this->leaf())
    {
        return this->purity;
    }
    else if (*attrs[this->attribute] > this->cut)
    {
        return this->get_right_child()->predict(attrs);
    }
    else
    {
        return this->get_left_child()->predict(attrs);
    }
}

void Node::recursive_split(unsigned int minLeafSize, unsigned int numCuts)
{
    if (this->split(minLeafSize,numCuts))
    {
        this->get_left_child()->recursive_split(minLeafSize,numCuts);
        this->get_right_child()->recursive_split(minLeafSize,numCuts);
    }
    else
    {
        this->calc_purity();
        this->update_classification();
    }
}

void Node::calc_purity()
{
    vector<Object*>::const_iterator it(this->signal.begin());
    double weightedSignal(0.);
    for (; it != this->signal.end(); ++it)
        weightedSignal += (*it)->weight;
    it = this->background.begin();
    double weightedBackground(0.);
    for (; it != this->background.end(); ++it)
        weightedBackground += (*it)->weight;
    this->purity = (weightedSignal + weightedBackground) > 0 ? weightedSignal / (weightedSignal + weightedBackground) : -1.;
}

void Node::update_classification()
{
    vector<Object*>::const_iterator it(this->signal.begin());
    if (this->purity > .5)
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

bool Node::split(unsigned int minSize, unsigned int resolution)
{
    if (signal.size()==0 || background.size()==0 || resolution < 2 || signal.size() + background.size() <= minSize)
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
    float bestGini = -1;
    for (unsigned int i(0); i < numAttributes; ++i)
    {
        // get max and min value of this attribute over signal and background combined
        extrema = this->minmax(i);
        // create histograms for signal and background
        double step = (extrema.second - extrema.first)/resolution;
        if (step==0) continue;
        sigHist = new Histogram<double,double>(resolution, extrema.first, extrema.second+step);
        bkgHist = new Histogram<double,double>(resolution, extrema.first, extrema.second+step);
        // fill histograms
        std::vector<Object*>::const_iterator it(signal.begin());
        for (; it != signal.end(); ++it)
            sigHist->fill((*it)->attribute(i),(*it)->get_weight());
        it = background.begin();
        for (; it != background.end(); ++it)
            bkgHist->fill((*it)->attribute(i),(*it)->get_weight());
        // calculate Gini_left + Gini_right for a split at each internal bin boundary and find minimum where minSize is respected
        double Gini_left;
        double Gini_right;
        double sig_left, sig_right;
        double bkg_left, bkg_right;
        double purity_left, purity_right;
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
    this->set_left_child(left);
    this->set_right_child(right);
    return true;
}

pair<double,double> Node::minmax(unsigned int attribute)
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

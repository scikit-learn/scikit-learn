#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <utility>

using namespace std;

template <typename T, typename U>
class Histogram
{
    public:

        Histogram(unsigned int nbins, U min, U max):
            nbins(nbins),
            min(min),
            max(max)
        {
            if (min >= max || nbins == 0) throw 0;
            bins = new vector<pair<T,T> >(nbins,pair<T,T>(0,0));
        }

        ~Histogram() { delete bins; }

        void fill(U value, T weight = 1)
        {
            if (value < this->min || value >= this->max) return;
            int bin = (this->nbins-1) * float(value - this->min) / float(this->max - this->min);
            this->bins->at(bin).first += weight;
            this->bins->at(bin).second += 1;
        }

        T integral(unsigned int lowBin = 0, unsigned int highBin = 0, bool weighted = true)
        {
            if (lowBin < 0) lowBin = 0;
            if (lowBin >= this->nbins) return 0;
            if ((highBin <= lowBin) || (highBin > this->nbins)) highBin = this->nbins;
            T sum = 0;
            typename vector<pair<T,T> >::const_iterator it(this->bins->begin()+lowBin);
            if (weighted)
                for (; it != this->bins->begin()+highBin; ++it) sum += it->first;
            else
                for (; it != this->bins->begin()+highBin; ++it) sum += it->second;
            return sum;
        }

        void reset()
        {
            this->bins->assign(this->nbins,pair<T,T>(0,0));
        }

    private:

        vector<pair<T,T> >* bins;
        unsigned int nbins;
        U min;
        U max;
};

#endif

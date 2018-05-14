/*
Passing variables / arrays between cython and cpp
Example from 
http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html

Adapted to include passing of multidimensional arrays

*/

#include "aode_helper.h"
#include <vector>

AODEHelper::AODEHelper(std::vector< std::vector< std::vector<int> > > ft, 
                       std::vector< std::vector< std::vector< std::vector< std::vector<int> > > > > jft, std::vector<int> uniq,
                       int na, int nc, int ns, int m, int m_lmt)
{
    freq_table = ft;
    joint_freq_table = jft;
    num_uniq_xvals = uniq;
    num_attr = na;
    num_classes = nc;
    num_samples = ns;
    m_val = m;
    m_limit = m_lmt;
}

AODEHelper::~AODEHelper()
{
}


std::vector< std::vector<double> > AODEHelper::calculate_naive_bayes(const std::vector< std::vector<int> > & preds, int & num_preds)
{
        std::vector< std::vector<double> > ret(num_preds, std::vector<double>(num_classes, 0));

        for (int pr = 0; pr < num_preds; pr++)
        {
                for (int c = 0; c < num_classes; c++)
                {
                        int count_y = 0;
                        for (int i = 0; i < num_uniq_xvals[0]; i++)
                        {
                                count_y += freq_table[c][0][i];
                        }

                        ret[pr][c] = (double)(count_y + ((double)m_val/(double)num_classes)) / (double)(num_samples+m_val);

                        for (int xx = 0; xx < num_attr; xx++)
                        {
                                int count_x_and_y = 0;
                                if (preds[pr][xx] != -1)
                                {
                                        count_x_and_y += freq_table[c][xx][preds[pr][xx]];
                                }

                                ret[pr][c] *= (double)(count_x_and_y + ((double)m_val/(double)num_uniq_xvals[xx])) / (double)(count_y+m_val);
                        }
                }

                // Normalize probabilities
                double sum_prob = 0;
                for (int c = 0; c < num_classes; c++)
                {
                        sum_prob += ret[pr][c];
                }
                for (int c = 0; c < num_classes; c++)
                {
                        ret[pr][c] = ret[pr][c] / sum_prob;
                }
        }

        return ret;
}

std::vector< std::vector<double> > AODEHelper::calculate_probabilities(const std::vector< std::vector<int> > & preds, int & num_preds)
{
        std::vector< std::vector<double> > ret(num_preds, std::vector<double>(num_classes, 0));

        for (int pr = 0; pr < num_preds; pr++)
        {
                std::vector< std::vector<double> > spode_probs(num_classes, std::vector<double>(num_attr, 0));
                int parent_count = 0;

                for (int xx = 0; xx < num_attr; xx++)
                {
                        int count_x1_and_y = 0;
                        for (int c = 0; c < num_classes; c++)
                        {
                                if (preds[pr][xx] != -1)
                                {
                                        count_x1_and_y += freq_table[c][xx][preds[pr][xx]];
                                }
                        }

                        // Check that attribute value has a frequency of m_limit or greater
                        if (count_x1_and_y < m_limit)
                        {
                                continue;
                        }
                        parent_count++;

                        for (int c = 0; c < num_classes; c++)
                        {
                                int xy_count = 0;
                                if (preds[pr][xx] != -1)
                                {
                                        xy_count = freq_table[c][xx][preds[pr][xx]];
                                }
                                spode_probs[c][xx] = (double)(xy_count + ((double)m_val/(double)(num_uniq_xvals[xx]*num_classes))) / (double)(num_samples+m_val);
                        }
                }

                // check parent_count, default to NB if 0
                if (parent_count <= 0)
                {
                        return calculate_naive_bayes(preds, num_preds);
                }

                for (int up = 0; up < num_attr; up++)
                {
                        for (int uc = 0; uc < up; uc++)
                        {
                                for (int c = 0; c < num_classes; c++)
                                {
                                        int xp_count = 0;
                                        if (preds[pr][up] != -1)
                                        {
                                                xp_count = freq_table[c][up][preds[pr][up]];
                                        }

                                        int xc_count = 0;
                                        if (preds[pr][uc] != -1)
                                        {
                                                xc_count = freq_table[c][uc][preds[pr][uc]];
                                        }

                                        int joint_count = 0;
                                        if (preds[pr][up] != -1 && preds[pr][uc] != -1)
                                        {
                                                joint_count = joint_freq_table[c][up][uc][preds[pr][up]][preds[pr][uc]];
                                        }

                                        spode_probs[c][uc] *= ((double)(joint_count + (m_val/num_uniq_xvals[uc]))/(double)(xc_count+m_val));
                                        spode_probs[c][up] *= ((double)(joint_count + (m_val/num_uniq_xvals[up]))/(double)(xp_count+m_val));
                                }
                        }
                }

                double total_probability = 0;
                for (int c = 0; c < num_classes; c++)
                {
                        for (int u = 0; u < num_attr; u++)
                        {
                                ret[pr][c] += spode_probs[c][u];
                        }

                        total_probability += ret[pr][c];
                }

                // Normalize probabilities
                for (int c = 0; c < num_classes; c++)
                {
                        ret[pr][c] = ret[pr][c] / total_probability;
                }

        }

        return ret;
}

std::vector<int> AODEHelper::calculate_predictions(const std::vector< std::vector<int> > & preds, int & num_preds)
{
        std::vector< std::vector<double> > probs = calculate_probabilities(preds, num_preds);

        std::vector<int> ret(num_preds, 0);

        for (int pr = 0; pr < num_preds; pr++)
        {
                int max_idx = 0;
                for (int c = 0; c < num_classes; c++)
                {
                        if (probs[pr][c] > probs[pr][max_idx])
                        {
                                max_idx = c;
                        }
                }

                ret[pr] = max_idx;
        }

        return ret;     
}


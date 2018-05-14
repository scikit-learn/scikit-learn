#include <vector>
#include <iostream>


class AODEHelper {
public:
    std::vector< std::vector< std::vector<int> > > freq_table;
    std::vector< std::vector< std::vector< std::vector< std::vector<int> > > > > joint_freq_table;
    std::vector<int> num_uniq_xvals;
    int num_attr, num_classes, num_samples, m_val, m_limit;
    AODEHelper(std::vector< std::vector< std::vector<int> > > ft, 
    std::vector< std::vector< std::vector< std::vector< std::vector<int> > > > > jft, std::vector<int> uniq,
    int na, int nc, int ns, int m, int m_lmt);
    ~AODEHelper();
    
    std::vector<int> calculate_predictions(const std::vector< std::vector<int> > & preds, int & num_preds);
    std::vector< std::vector<double> > calculate_probabilities(const std::vector< std::vector<int> > & preds, int & num_preds);
    std::vector< std::vector<double> > calculate_naive_bayes(const std::vector< std::vector<int> > & preds, int & num_preds);
};
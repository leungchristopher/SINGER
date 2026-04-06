//
//  Coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 4/17/23.
//  Modified for ChrisSINGER: flat arrays, optional lazy mode.
//

#ifndef Coalescent_calculator_hpp
#define Coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include <vector>
#include <algorithm>
#include "Branch.hpp"

class Coalescent_calculator {

public:

    double cut_time;
    double min_time, max_time;
    bool lazy_mode = false;

    // Flat sorted vectors replace std::map/std::set
    vector<pair<double, int>> rate_changes = {};
    vector<pair<double, int>> rates = {};
    vector<pair<double, double>> probs = {};      // sorted by .first (time)
    vector<pair<double, double>> quantiles = {};   // sorted by .second (cum_prob)

    Coalescent_calculator(double t);

    ~Coalescent_calculator();

    void compute(set<Branch> &branches);

    double weight(double lb, double ub);

    double time(double lb, double ub);

    // Lazy incremental update: modify rate_changes in-place
    void incremental_update(Branch removed, Branch inserted);

// private:

    void compute_rate_changes(set<Branch> &branches);

    void compute_rates();

    void compute_probs_quantiles();

    double prob(double x);

    double quantile(double p);

};

#endif /* Coalescent_calculator_hpp */

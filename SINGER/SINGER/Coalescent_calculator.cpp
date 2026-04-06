//
//  Coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 4/17/23.
//  Modified for ChrisSINGER: flat sorted vectors, optional lazy mode.
//

#include "Coalescent_calculator.hpp"

Coalescent_calculator::Coalescent_calculator(double t) {
    cut_time = t;
}

Coalescent_calculator::~Coalescent_calculator() {}

void Coalescent_calculator::compute(set<Branch> &branches) {
    compute_rate_changes(branches);
    compute_rates();
    assert(rates.size() == rate_changes.size());
    compute_probs_quantiles();
    assert(rates.size() == rate_changes.size());
}

double Coalescent_calculator::weight(double lb, double ub) {
    double p = prob(ub) - prob(lb);
    assert(!isnan(p));
    return p;
}

double Coalescent_calculator::time(double lb, double ub) {
    if (isinf(ub)) {
        return lb + log(2);
    }
    if (ub - lb < 1e-3) {
        return 0.5 * (lb + ub);
    }
    double lq = prob(lb);
    double uq = prob(ub);
    double mid = 0;
    double t;
    if (uq - lq < 1e-3) {
        t = 0.5 * (lb + ub);
    } else {
        mid = 0.5 * (lq + uq);
        t = quantile(mid);
    }
    assert(t >= lb and t <= ub);
    return t;
}

void Coalescent_calculator::compute_rate_changes(set<Branch> &branches) {
    // Build into a temporary map, then flatten to sorted vector
    map<double, int> rc_map;
    double lb = 0;
    double ub = 0;
    for (Branch b : branches) {
        lb = max(cut_time, b.lower_node->time);
        ub = b.upper_node->time;
        rc_map[lb] += 1;
        rc_map[ub] -= 1;
    }
    rate_changes.assign(rc_map.begin(), rc_map.end());
    min_time = rate_changes.front().first;
    max_time = rate_changes.back().first;
}

void Coalescent_calculator::compute_rates() {
    rates.clear();
    rates.reserve(rate_changes.size());
    int curr_rate = 0;
    for (auto &x : rate_changes) {
        curr_rate += x.second;
        rates.emplace_back(x.first, curr_rate);
    }
    assert(rates.size() == rate_changes.size());
}

void Coalescent_calculator::compute_probs_quantiles() {
    probs.clear();
    quantiles.clear();
    probs.reserve(rates.size());
    quantiles.reserve(rates.size());

    int curr_rate = 0;
    double prev_time = 0.0;
    double next_time = 0.0;
    double prev_prob = 1.0;
    double next_prob = 1.0;
    double cum_prob = 0.0;

    for (size_t i = 0; i + 1 < rates.size(); i++) {
        curr_rate = rates[i].second;
        prev_time = rates[i].first;
        next_time = rates[i + 1].first;
        if (curr_rate > 0) {
            next_prob = prev_prob * exp(-curr_rate * (next_time - prev_time));
            cum_prob += (prev_prob - next_prob) / curr_rate;
        } else {
            next_prob = prev_prob;
        }
        assert(cum_prob <= 1.0001);
        probs.push_back({next_time, cum_prob});
        quantiles.push_back({next_time, cum_prob});
        prev_prob = next_prob;
    }
    probs.push_back({min_time, 0});
    quantiles.push_back({min_time, 0});

    // Sort probs by time (first), quantiles by cum_prob (second) then time
    sort(probs.begin(), probs.end(), [](const pair<double,double> &a, const pair<double,double> &b) {
        return a.first < b.first;
    });
    sort(quantiles.begin(), quantiles.end(), [](const pair<double,double> &a, const pair<double,double> &b) {
        if (a.second != b.second) return a.second < b.second;
        return a.first < b.first;
    });

    assert(probs.size() == rates.size());
    assert(quantiles.size() == rates.size());
}

double Coalescent_calculator::prob(double x) {
    if (x > max_time) {
        x = max_time;
    } else if (x < min_time) {
        x = min_time;
    }
    // Binary search in probs (sorted by .first = time)
    auto it = lower_bound(probs.begin(), probs.end(), make_pair(x, -1e300),
                          [](const pair<double,double> &a, const pair<double,double> &b) {
                              return a.first < b.first;
                          });
    // Exact match
    if (it != probs.end() && it->first == x) {
        return it->second;
    }
    // Interpolate between surrounding entries
    auto u_it = it; // first entry with time > x
    auto l_it = prev(it);
    double base_prob = l_it->second;
    // Find rate at l_it->first using binary search in rates
    auto rate_it = lower_bound(rates.begin(), rates.end(), make_pair(l_it->first, 0),
                               [](const pair<double,int> &a, const pair<double,int> &b) {
                                   return a.first < b.first;
                               });
    int rate = (rate_it != rates.end() && rate_it->first == l_it->first) ? rate_it->second : 0;
    if (rate == 0) {
        return base_prob;
    }
    double delta_t = u_it->first - l_it->first;
    double delta_p = u_it->second - l_it->second;
    double new_delta_t = x - l_it->first;
    double new_delta_p = delta_p * expm1(-rate * new_delta_t) / expm1(-rate * delta_t);
    assert(!isnan(new_delta_p));
    return base_prob + new_delta_p;
}

double Coalescent_calculator::quantile(double p) {
    // Binary search in quantiles (sorted by .second = cum_prob)
    auto u_it = upper_bound(quantiles.begin(), quantiles.end(), make_pair(-1.0, p),
                            [](const pair<double,double> &a, const pair<double,double> &b) {
                                if (a.second != b.second) return a.second < b.second;
                                return a.first < b.first;
                            });
    auto l_it = prev(u_it);
    double base_time = l_it->first;
    // Find rate at l_it->first
    auto rate_it = lower_bound(rates.begin(), rates.end(), make_pair(l_it->first, 0),
                               [](const pair<double,int> &a, const pair<double,int> &b) {
                                   return a.first < b.first;
                               });
    int rate = (rate_it != rates.end() && rate_it->first == l_it->first) ? rate_it->second : 0;
    double delta_t = u_it->first - l_it->first;
    double delta_p = u_it->second - l_it->second;
    double new_delta_p = p - l_it->second;
    double new_delta_t = 1 - new_delta_p / delta_p * (1 - exp(-rate * delta_t));
    new_delta_t = -log(new_delta_t) / rate;
    assert(!isnan(new_delta_t));
    return base_time + new_delta_t;
}

void Coalescent_calculator::incremental_update(Branch removed, Branch inserted) {
    if (!lazy_mode) {
        return;
    }
    // Remove the old branch's contribution from rate_changes
    double old_lb = max(cut_time, removed.lower_node->time);
    double old_ub = removed.upper_node->time;
    // Add the new branch's contribution
    double new_lb = max(cut_time, inserted.lower_node->time);
    double new_ub = inserted.upper_node->time;

    // Use a temporary map for modifications, then re-flatten
    map<double, int> rc_map(rate_changes.begin(), rate_changes.end());
    rc_map[old_lb] -= 1;
    rc_map[old_ub] += 1;
    rc_map[new_lb] += 1;
    rc_map[new_ub] -= 1;

    // Clean up zero entries
    for (auto it = rc_map.begin(); it != rc_map.end(); ) {
        if (it->second == 0) {
            it = rc_map.erase(it);
        } else {
            ++it;
        }
    }

    rate_changes.assign(rc_map.begin(), rc_map.end());
    if (!rate_changes.empty()) {
        min_time = rate_changes.front().first;
        max_time = rate_changes.back().first;
    }

    compute_rates();
    compute_probs_quantiles();
}

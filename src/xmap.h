#ifndef XMAP_H
#define XMAP_H

#include <Rcpp.h>
#include <iostream>
#include <random>
#include "forecast_machine.h"

using namespace Rcpp;

class Xmap: public ForecastMachine
{
public:
    // *** constructors *** //
    Xmap();
    
    // *** methods *** //
    void set_time(const NumericVector time);
    void set_block(const NumericMatrix new_block);
    void set_norm_type(const int norm_type);
    void set_lib(const NumericMatrix lib);
    void set_pred(const NumericMatrix pred);
    void set_lib_sizes(const NumericVector new_lib_sizes);
    void set_exclusion_radius(const double new_exclusion_radius);
    void set_epsilon(const double new_epsilon);
    void set_p(const double new_p);
    void set_lib_column(const size_t new_lib_col);
    void set_target_column(const size_t new_target);
    void set_params(const size_t new_E, const size_t new_tau, const int new_tp, 
                    const size_t new_nn, const bool new_random_libs, 
                    const size_t new_num_samples, const bool new_replace);
    void set_seed(const size_t new_seed);
    void suppress_warnings();
    void gpr();
    void run();
    DataFrame get_output();
    
private:
    void prepare_forecast();
    void make_vectors();
    void make_targets();
    
    // *** local parameters *** //
    std::vector<vec> block;
    std::vector<size_t> lib_sizes;
    int tp;
    size_t E, tau;
    size_t lib_col, target;
    bool random_libs;
    size_t num_samples;
    size_t seed;
    bool replace;
    bool remake_vectors;
    bool remake_targets;
    bool remake_ranges;
    
    // *** output data structures *** //
    std::vector<PredStats> predicted_stats;
    std::vector<size_t> predicted_lib_sizes;
};

#endif

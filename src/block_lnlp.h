#ifndef BLOCK_LNLP_H
#define BLOCK_LNLP_H

#include <Rcpp.h>
#include <iostream>
#include "forecast_machine.h"

using namespace Rcpp;

class BlockLNLP: public ForecastMachine
{
public:
    // *** constructors *** //
    BlockLNLP();
    
    // *** methods *** //
    void set_time(const NumericVector time);
    void set_block(const NumericMatrix new_block);
    void set_norm_type(const int norm_type);
    void set_pred_type(const int pred_type);
    void set_lib(const NumericMatrix lib);
    void set_pred(const NumericMatrix pred);
    void set_exclusion_radius(const double new_exclusion_radius);
    void set_epsilon(const double new_epsilon);
    void set_embedding(const NumericVector new_embedding);
    void set_target_columns(const NumericVector new_targets);
    void set_params(const int new_tp, const size_t new_nn);
    void set_theta(const double new_theta);
    void suppress_warnings();
    void save_smap_coefficients();
    void run();
    DataFrame get_output();
    List get_smap_coefficients();
    DataFrame get_short_output(const size_t target_idx);
    DataFrame get_stats(const size_t target_idx);
    
private:
    void prepare_forecast();
    void make_vectors();
    void make_targets();
    std::vector<vec> block;
    
    // *** local parameters *** //
    int tp;
    size_t E;
    std::vector<size_t> embedding;
    std::vector<size_t> target_cols;
    bool remake_vectors;
    bool remake_targets;
    bool remake_ranges;
};

#endif

#ifndef FORECAST_MACHINE_H
#define FORECAST_MACHINE_H

#include <iostream>
#include <vector>
#include <numeric>
#include <thread>
#include <stdexcept>
#include <math.h>
#include <Rcpp.h>
#include "data_types.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace Rcpp;

class ForecastMachine
{
protected:
    // *** constructors *** //
    ForecastMachine();
    
    // *** computational methods *** //
    void init_distances();
    void compute_distances();
    //void sort_neighbors();
    std::vector<size_t> find_nearest_neighbors(const size_t curr_pred, const std::vector<bool>& valid_lib_indices);

    void forecast();
    void set_indices_from_range(std::vector<bool>& indices, const std::vector<time_range>& range, 
                                  int start_shift, int end_shift, bool check_target);
    void check_cross_validation();
    bool is_vec_valid(const size_t vec_index);
    bool is_target_valid(const size_t vec_index);
    PredStats compute_stats(const vec& obs, const vec& pred);
    PredStats make_stats(const size_t target_idx);
    PredStats make_const_stats(const size_t target_idx);
    void LOG_WARNING(const char* warning_text);
    
    // *** variables *** //
    std::vector<bool> lib_indices;
    std::vector<bool> pred_indices;
    std::vector<size_t> which_lib;
    std::vector<size_t> which_pred;
    
    vec time;
    vec target_time;
    std::vector<vec> data_vectors;
    std::vector<vec> smap_coefficients;
    size_t num_targets;
    std::vector<vec> targets;
    std::vector<vec> predicted;
    std::vector<vec> predicted_var;
    std::vector<vec> const_targets;
    std::vector<vec> const_predicted;
    size_t num_vectors;
    double (*dist_func)(const vec&, const vec&);
    std::vector<vec> distances;
    
    // *** parameters *** //
    bool CROSS_VALIDATION;
    bool SUPPRESS_WARNINGS;
    bool SAVE_SMAP_COEFFICIENTS;
    PredEnum pred_mode;
    NormEnum norm_mode;
    size_t nn;
    double theta;
    double exclusion_radius;
    double epsilon;
    std::vector<time_range> lib_ranges;
    std::vector<time_range> pred_ranges;
    static const double qnan;
    
private:
    // *** methods *** //
    void simplex_forecast();
    void smap_forecast();
    void simplex_prediction(const size_t start, const size_t end);
    void smap_prediction(const size_t start, const size_t end);
    void const_prediction(const size_t start, const size_t end);
    std::vector<bool> adjust_lib(const size_t curr_pred);
    
    //int num_threads;
};

std::vector<size_t> which_indices_true(const std::vector<bool>& indices);
double l1_distance_func(const vec& A, const vec& B);
double l2_distance_func(const vec& A, const vec& B);
std::vector<size_t> sort_indices(const std::vector<double>& v, const std::vector<size_t> idx);

#endif

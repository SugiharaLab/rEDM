#ifndef FORECAST_MACHINE_H
#define FORECAST_MACHINE_H

#include <iostream>
#include <vector>
#include <numeric>
#include <thread>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>
#include "data_types.h"
//#include <Eigen/Dense>
#include <RcppEigen.h>

//using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
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
    std::vector<size_t> find_nearest_neighbors(const vec& dist);

    void forecast();
    void set_indices_from_range(std::vector<bool>& indices, const std::vector<time_range>& range, 
                                  int start_shift, int end_shift, bool check_target);
    void check_cross_validation();
    bool is_vec_valid(const size_t vec_index);
    bool is_target_valid(const size_t vec_index);
    PredStats make_stats();
    PredStats make_const_stats();
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
    vec targets;
    vec predicted;
    vec predicted_var;
    vec const_targets;
    vec const_predicted;
    size_t num_vectors;
    std::function<double (const vec&, const vec&)> dist_func;
    std::vector<vec > distances;
    
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
    double p;
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
    void adjust_lib(const size_t curr_pred);
    
    //int num_threads;
};

std::vector<size_t> which_indices_true(const std::vector<bool>& indices);
std::vector<size_t> sort_indices(const std::vector<double>& v, const std::vector<size_t> idx);
PredStats compute_stats_internal(const vec& obs, const vec& pred);
DataFrame get_stats(const vec& obs, const vec& pred);

#endif

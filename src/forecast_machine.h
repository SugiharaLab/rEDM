#ifndef FORECAST_MACHINE_H
#define FORECAST_MACHINE_H

#include <iostream>
#include <vector>
#include <numeric>
#include <thread>
#include <stdexcept>
#include <math.h>
#include "data_types.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class ForecastMachine
{
protected:
    // *** constructors *** //
    ForecastMachine();    
    
    // *** computational methods *** //
    void init_distances();
    void compute_distances();
    void sort_neighbors();
    vector<size_t> find_nearest_neighbors(const size_t curr_pred, const vector<bool>& valid_lib_indices);

    void forecast();
    void set_indices_from_range(vector<bool>& indices, const vector<time_range>& range, 
                                  int start_shift, int end_shift, bool check_target);
    void check_cross_validation();
    bool is_vec_valid(const size_t vec_index);
    bool is_target_valid(const size_t vec_index);
    PredStats compute_stats();
    void LOG_WARNING(const char* warning_text);
    
    // *** variables *** //
    vector<bool> lib_indices;
    vector<bool> pred_indices;
    vector<size_t> which_lib;
    vector<size_t> which_pred;
    
    vector<double> time;
    vector<vec> data_vectors;
    vector<double> observed;
    vector<double> predicted;
    size_t num_vectors;
    double (*dist_func)(const vec&, const vec&);
    vector<vector<double> > distances;
    vector<vector<size_t> > neighbors;
    
    // *** parameters *** //
    bool CROSS_VALIDATION;
    bool SUPPRESS_WARNINGS;
    PredEnum pred_mode;
    NormEnum norm_mode;
    int nn;
    double theta;
    double exclusion_radius;
    vector<time_range> lib_ranges;
    vector<time_range> pred_ranges;
    static const double qnan;
    
private:
    // *** methods *** //
    void simplex_forecast();
    void smap_forecast();
    void simplex_prediction(const size_t start, const size_t end);
    void smap_prediction(const size_t start, const size_t end);
    vector<bool> adjust_lib(const size_t curr_pred);
    
    int num_threads;
};

vector<size_t> which_indices_true(const vector<bool>& indices);
double l1_distance_func(const vec& A, const vec& B);
double l2_distance_func(const vec& A, const vec& B);
vector<size_t> sort_indices(const vector<double>& v, const vector<size_t> idx);

#endif
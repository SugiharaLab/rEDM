#ifndef FORECAST_MACHINE_H
#define FORECAST_MACHINE_H

#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>
#include "data_types.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class ForecastMachine
{
protected:
    // *** constructors *** //
    ForecastMachine();
    
    // *** methods *** //
    void init_distances();
    void compute_distances();
    void sort_neighbors();
    vector<size_t> find_nearest_neighbors(const size_t curr_pred);
    void forecast();
    bool is_vec_valid(const size_t vec_index);
    void LOG_WARNING(const char* warning_text);
    
    // *** variables *** //
    vector<bool> valid_lib_indices;
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
    void simplex_prediction(const size_t curr_pred);
    void smap_prediction(const size_t curr_pred);
    
    void adjust_lib();
    void adjust_lib(const size_t curr_pred);
};

vector<size_t> which_indices_true(const vector<bool>& indices);
double l1_distance_func(const vec& A, const vec& B);
double l2_distance_func(const vec& A, const vec& B);
vector<size_t> sort_indices(const vector<double>& v, const vector<size_t> idx);

#endif
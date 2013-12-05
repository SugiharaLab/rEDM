#ifndef FORECAST_MACHINE_H
#define FORECAST_MACHINE_H

#include <iostream>
#include <vector>
#include <math.h>
#include "data_types.h"

using namespace std;

class ForecastMachine
{
protected:
    // *** constructors *** //
    ForecastMachine();
    
    // *** methods *** //
    void setup_lib_and_pred();
    void forecast();
    
    
    // *** variables *** //
    vector<time_range> lib;
    vector<time_range> pred;
    vector<bool> valid_lib_indices;
    vector<bool> lib_indices;
    vector<bool> pred_indices;
    vector<vec> data_vectors;
    vector<double> target_vals;
    vector<double> predicted;
    int num_vectors;
    
    // *** parameters *** //
    bool CROSS_VALIDATION;
    PredEnum pred_mode;
    NormEnum norm_mode;
    
private:
    // *** methods *** //
    void set_indices_from_range(vector<bool>& indices, const vector<time_range>& range);
    void compute_distances();
    
    // forecasting
    void simplex_forecast();
    void smap_forecast();
    
    void make_lib();
    void make_lib(const int curr_pred);
    
    void simplex_prediction(const int curr_pred);
    void smap_prediction(const int curr_pred);
    
    bool is_vec_valid(const int vec_index);
    void LOG_WARNING(const char* warning_text);
    
};

vector<int> which_indices_true(const vector<bool>& indices);
double l1_distance_func(const vec& A, const vec& B);
double l2_distance_func(const vec& A, const vec& B);

#endif
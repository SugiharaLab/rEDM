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
    void forecast();
    void LOG_WARNING(const char* warning_text);
    
    // *** variables *** //
    vector<bool> valid_lib_indices;
    vector<bool> lib_indices;
    vector<bool> pred_indices;
    vector<int> which_lib;
    vector<int> which_pred;
    
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
    void compute_distances();
    
    // forecasting
    void simplex_forecast();
    void smap_forecast();
    
    void make_lib();
    void make_lib(const int curr_pred);
    
    void simplex_prediction(const int curr_pred);
    void smap_prediction(const int curr_pred);
    
    bool is_vec_valid(const int vec_index);
    
    // *** data *** //
    vector<vector<double> > distances;
    
};

vector<int> which_indices_true(const vector<bool>& indices);
double l1_distance_func(const vec& A, const vec& B);
double l2_distance_func(const vec& A, const vec& B);

#endif
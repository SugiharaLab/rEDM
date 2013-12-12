#ifndef LNLP_H
#define LNLP_H

#include <Rcpp.h>
#include "forecast_machine.h"

using namespace std;
using namespace Rcpp;

class LNLP: public ForecastMachine
{
public:
    // *** constructors *** //
    LNLP();
    LNLP(const NumericVector data, const NumericMatrix lib, const NumericMatrix pred,
         const int norm_type, const int pred_type, const double exclusion_radius,
         const int E, const int tau, const int tp, const int nn);
    
    // *** methods *** //
    void set_time(const NumericVector time);
    void set_time_series(const NumericVector data);
    void set_norm_type(const int norm_type);
    void set_pred_type(const int pred_type);
    void set_lib(const NumericMatrix lib);
    void set_pred(const NumericMatrix pred);
    void set_exclusion_radius(const double new_exclusion_radius);
    void set_params(const int new_E, const int new_tau, const int new_tp, const int new_nn);
    void run();
    DataFrame get_output();
    DataFrame get_stats();
    
private:
    void prepare_forecast();
    void make_vectors();
    void make_targets();
    void set_indices_from_range(vector<bool>& indices, const vector<time_range>& range);
    void check_cross_validation();
    vector<double> time_series;
    
    // *** local parameters *** //
    int E, tau, max_lag, tp;
    bool remake_vectors;
    bool remake_targets;
    bool remake_ranges;
};

#endif
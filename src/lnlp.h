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
    
    // *** methods *** //
    void set_time(const NumericVector new_time);
    void set_time_series(const NumericVector data);
    void set_norm_type(const int norm_type);
    void set_pred_type(const int pred_type);
    void set_lib(const NumericMatrix lib);
    void set_pred(const NumericMatrix pred);
    void set_exclusion_radius(const double new_exclusion_radius);
    void set_params(const int new_E, const int new_tau, const int new_tp, const int new_nn);
    void set_theta(const double new_theta);
    void run();
    DataFrame get_output();
    DataFrame get_stats();
    
private:
    void prepare_forecast();
    void make_vectors();
    void make_targets();
    void check_cross_validation();
    vector<double> time_series;
    
    // *** local parameters *** //
    int E, tau, tp;
    bool remake_vectors;
    bool remake_targets;
    bool remake_ranges;
};

#endif
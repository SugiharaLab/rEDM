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
           const NumericVector E, const NumericVector tau, const NumericVector tp, 
           const NumericVector num_neighbors);
    
    // *** methods *** //
    void make_vectors();
    void run();
    void set_time(const NumericVector time);
    void set_time_series(const NumericVector data);
    void set_norm_type(const int norm_type);
    void set_lib(const NumericMatrix lib);
    void set_pred(const NumericMatrix pred);
    void set_params(const int E, const int tau, const int tp, const int nn);
    
  private:
    vector<double> time_series;
    
    
};

#endif
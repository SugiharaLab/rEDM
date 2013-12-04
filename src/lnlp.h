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
    LNLP(const NumericVector ts, const NumericMatrix lib, const NumericMatrix pred, 
           const int norm_type, const int pred_type, const double exclusion_radius, 
           const NumericVector E, const NumericVector tau, const NumericVector tp, 
           const NumericVector num_neighbors);
    
    // *** methods *** //
    void make_vectors();
    void run();
    void set_time_series(NumericVector ts);
    
  private:
    vector<double> time_series;
    
    
};

#endif
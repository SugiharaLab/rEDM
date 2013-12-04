#include "lnlp.h"

/*** Constructors ***/
LNLP::LNLP()
{
}

LNLP::LNLP(const NumericVector ts, const NumericMatrix lib, const NumericMatrix pred, 
           const int norm_type, const int pred_type, const double exclusion_radius, 
           const NumericVector E, const NumericVector tau, const NumericVector tp, 
           const NumericVector num_neighbors)
{
}

void LNLP::make_vectors()
{
  return;
}

void LNLP::run()
{
  
  return;
}

void LNLP::set_time_series(const NumericVector ts)
{
  return;
}



RCPP_MODULE(lnlp_module)
{
  class_<LNLP>("LNLP")
  
  .constructor()
  ;
}
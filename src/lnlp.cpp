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

void LNLP::set_time(const NumericVector time)
{
    return;
}

void LNLP::set_time_series(const NumericVector data)
{
  return;
}

void LNLP::set_norm_type(const int norm_type)
{
    return;
}

void LNLP::set_lib(const NumericMatrix lib)
{
    return;
}

void LNLP::set_pred(const NumericMatrix pred)
{
    return;
}

void LNLP::set_params(const int E, const int tau, const int tp, const int nn)
{
    return;
}


RCPP_MODULE(lnlp_module)
{
  class_<LNLP>("LNLP")
  
  .constructor()
  
  .method( "run", &LNLP::run)
  .method( "set_time", &LNLP::set_time)
  .method( "set_time_series", &LNLP::set_time_series)
  .method( "set_norm_type", &LNLP::set_norm_type)
  .method( "set_lib", &LNLP::set_lib)
  .method( "set_pred", &LNLP::set_pred)
  .method( "set_params", &LNLP::set_params)
;
}
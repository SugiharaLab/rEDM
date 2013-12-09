#include "lnlp.h"

/*** Constructors ***/
LNLP::LNLP()
{
}

LNLP::LNLP(const NumericVector ts, const NumericMatrix lib, const NumericMatrix pred, 
           const int norm_type, const int pred_type, const double exclusion_radius, 
           const int E, const int tau, const int tp, const int num_neighbors)
{
}

void LNLP::make_vectors()
{
  return;
}

void LNLP::run()
{
    // check parameters
    Rcpp::stop("check parameters are valid.");  
    // forecast
    forecast();
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
    set_indices_from_range(lib_indices, lib);
    which_lib = which_indices_true(lib_indices);
    return;
}

void LNLP::set_pred(const NumericMatrix pred)
{
    set_indices_from_range(pred_indices, pred);
    which_pred = which_indices_true(pred_indices);
    return;
}

void LNLP::set_params(const int E, const int tau, const int tp, const int nn)
{
    return;
}

void LNLP::set_indices_from_range(vector<bool>& indices, const NumericMatrix range)
{
    int start_of_range, end_of_range;
    indices.assign(num_vectors, false); // initialize indices
    for(int i = 0; i < range.nrow(); ++i)
    {
        start_of_range = range(i, 1) - 1 - max_lag;
        end_of_range = range(i, 2) - 1 - max(0, tp);
        if(start_of_range < 0) // check start of range
        {
            LOG_WARNING("beginning of time_range was less than 1; corrected.");
            start_of_range = 0;
        }
        if(end_of_range > (num_vectors-1)) // check end of range
        {
            LOG_WARNING("end of time_range was greater than the number of vectors; corrected");
            end_of_range = num_vectors-1;
        }
        
        for(int j = start_of_range; j <= end_of_range; j++)
            indices[j] = true;
	}
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
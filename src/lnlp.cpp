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
    prepare_forecast();
    
    // forecast
    forecast(); // forecast code is in forecast_machine
    return;
}

void LNLP::prepare_forecast()
{
    num_vectors = time_series.size();

    compute_distances();

    set_indices_from_range(lib_indices, lib_ranges);
    which_lib = which_indices_true(lib_indices);

    set_indices_from_range(pred_indices, pred_ranges);
    which_pred = which_indices_true(pred_indices);


    return;
}


void LNLP::set_time(const NumericVector time)
{
    this->time = as<std::vector<double> >(time);
    return;
}

void LNLP::set_time_series(const NumericVector data)
{
    time_series = as<std::vector<double> >(data);
    return;
}

void LNLP::set_norm_type(const int norm_type)
{
    switch(norm_type)
    {
        case 1:
            norm_mode = L1_NORM;
            break;
        case 2:
            norm_mode = L2_NORM;
            break;
        default:
            Rcpp::stop("unknown norm type selected");
    }
    return;
}

void LNLP::set_pred_type(const int pred_type)
{
    switch(pred_type)
    {
        case 1:
            pred_mode = SMAP;
            break;
        case 2:
            pred_mode = SIMPLEX;
            break;
        case 3:
            pred_mode = FAST_LINEAR;
            break;
        default:
            Rcpp::stop("unknown pred type selected");
    }
    return;
}

void LNLP::set_lib(const NumericMatrix lib)
{
    lib_ranges.resize(lib.nrow());
    for(int i = 0; i < lib.nrow(); ++i)
    {
        lib_ranges[i].first = lib(i,1) - 1; // convert 1-index to 0-index
        lib_ranges[i].second = lib(i,2) - 1;
    }
    return;
}

void LNLP::set_pred(const NumericMatrix pred)
{
    pred_ranges.resize(pred.nrow());
    for(int i = 0; i < pred.nrow(); ++i)
    {
        pred_ranges[i].first = pred(i,1) - 1; // convert 1-index to 0-index
        pred_ranges[i].second = pred(i,2) - 1;
    }
    return;
}

void LNLP::set_params(const int E, const int tau, const int tp, const int nn)
{
    return;
}

void LNLP::set_indices_from_range(vector<bool>& indices, const vector<time_range>& range)
{
    int start_of_range, end_of_range;
    indices.assign(num_vectors, false); // initialize indices
    for(auto& range_iter: range)
    {
        start_of_range = range_iter.first - max_lag;
        end_of_range = range_iter.second - max(0, tp);
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
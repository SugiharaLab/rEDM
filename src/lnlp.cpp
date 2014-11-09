#include "lnlp.h"

/*** Constructors ***/
LNLP::LNLP(): remake_vectors(true), remake_targets(true), remake_ranges(true)
{
}

void LNLP::set_time(const NumericVector new_time)
{
    time = as<std::vector<double> >(new_time);
    return;
}

void LNLP::set_time_series(const NumericVector data)
{
    time_series = as<std::vector<double> >(data);
    num_vectors = time_series.size();
    init_distances();
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
            throw(std::domain_error("unknown norm type selected"));
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
            throw(std::domain_error("unknown pred type selected"));
    }
    return;
}

void LNLP::set_lib(const NumericMatrix lib)
{
    size_t num_rows = size_t(lib.nrow());
    lib_ranges.resize(num_rows);
    for(size_t i = 0; i < num_rows; ++i)
    {
        lib_ranges[i].first = lib(i,0) - 1; // convert 1-index to 0-index
        lib_ranges[i].second = lib(i,1) - 1;
    }
    remake_ranges = true;
    return;
}

void LNLP::set_pred(const NumericMatrix pred)
{
    size_t num_rows = size_t(pred.nrow());
    pred_ranges.resize(num_rows);
    for(size_t i = 0; i < num_rows; ++i)
    {
        pred_ranges[i].first = pred(i,0) - 1; // convert 1-index to 0-index
        pred_ranges[i].second = pred(i,1) - 1;
    }
    remake_ranges = true;
    return;
}

void LNLP::set_exclusion_radius(const double new_exclusion_radius)
{
    exclusion_radius = new_exclusion_radius;
    if(exclusion_radius >= 0)
        CROSS_VALIDATION = true;
    return;
}

void LNLP::set_params(const int new_E, const int new_tau, const int new_tp, const int new_nn)
{
    if(E != new_E || tau != new_tau)
        remake_vectors = true;
    if(tp != new_tp)
        remake_targets = true;
    if(remake_vectors || remake_targets)
        remake_ranges = true;
    
    E = new_E;
    tau = new_tau;
    tp = new_tp;
    nn = new_nn;    
    return;
}

void LNLP::set_theta(const double new_theta)
{
    theta = new_theta;
    return;
}

void LNLP::suppress_warnings()
{
    SUPPRESS_WARNINGS = true;
    return;
}


void LNLP::run()
{
    prepare_forecast(); // check parameters
    forecast(); // forecast code is in forecast_machine
    return;
}

DataFrame LNLP::get_output()
{
    return DataFrame::create( Named("time") = time, 
                              Named("obs") = observed, 
                              Named("pred") = predicted, 
                              Named("pred_var") = predicted_var);
}

DataFrame LNLP::get_short_output()
{
    vec short_time(which_pred.size(), qnan);
    vec short_obs(which_pred.size(), qnan);
    vec short_pred(which_pred.size(), qnan);
    
    for(size_t i = 0; i < which_pred.size(); ++i)
    {
        short_time[i] = time[which_pred[i]];
        short_obs[i] = observed[which_pred[i]];
        short_pred[i] = predicted[which_pred[i]];
    }
    
    return DataFrame::create( Named("time") = short_time, 
                              Named("obs") = short_obs, 
                              Named("pred") = short_pred);
}

DataFrame LNLP::get_stats()
{
    PredStats output = compute_stats();
    return DataFrame::create( Named("num_pred") = output.num_pred, 
                              Named("rho") = output.rho, 
                              Named("mae") = output.mae, 
                              Named("rmse") = output.rmse );
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void LNLP::prepare_forecast()
{
    if(remake_vectors)
    {
        make_vectors();
        init_distances();
    }
    
    if(remake_targets)
        make_targets();
    
    if(remake_ranges)
    {
        set_indices_from_range(lib_indices, lib_ranges, (E-1)*tau, -std::max(0, tp), true);
        set_indices_from_range(pred_indices, pred_ranges, (E-1)*tau, -std::max(0, tp), false);

        check_cross_validation();

        which_lib = which_indices_true(lib_indices);
        which_pred = which_indices_true(pred_indices);
        
        remake_ranges = false;
    }
    
    compute_distances();
    return;
}

void LNLP::make_vectors()
{
    data_vectors.assign(num_vectors, vec(E, qnan));

    // beginning of lagged vectors cannot lag before start of time series
    for(size_t i = 0; i < (unsigned int)((E-1)*tau); ++i)
        for(size_t j = 0; j < (unsigned int)(E); ++j)
            if(i >= j*tau)
                data_vectors[i][j] = time_series[i - j * tau];
    
    // remaining lagged vectors
    for(size_t i = (unsigned int)((E-1)*tau); i < num_vectors; ++i)
        for(size_t j = 0; j < (unsigned int)(E); ++j)
            data_vectors[i][j] = time_series[i - j * tau];
            
    remake_vectors = false;
    return;
}

void LNLP::make_targets()
{
    observed.assign(time_series.begin()+tp, time_series.end());
    observed.insert(observed.end(), tp, qnan);
    
    remake_targets = false;
    return;
}

RCPP_MODULE(lnlp_module)
{
    class_<LNLP>("LNLP")
    
    .constructor()
    
    .method("set_time", &LNLP::set_time)
    .method("set_time_series", &LNLP::set_time_series)
    .method("set_norm_type", &LNLP::set_norm_type)
    .method("set_pred_type", &LNLP::set_pred_type)
    .method("set_lib", &LNLP::set_lib)
    .method("set_pred", &LNLP::set_pred)
    .method("set_exclusion_radius", &LNLP::set_exclusion_radius)
    .method("set_params", &LNLP::set_params)
    .method("set_theta", &LNLP::set_theta)
    .method("suppress_warnings", &LNLP::suppress_warnings)
    .method("run", &LNLP::run)
    .method("get_output", &LNLP::get_output)
    .method("get_short_output", &LNLP::get_short_output)
    .method("get_stats", &LNLP::get_stats)
    ;
}
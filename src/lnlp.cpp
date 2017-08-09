#include "lnlp.h"

/*** Constructors ***/
LNLP::LNLP(): 
    time_series(vec()), tp(1), E(1), tau(1), 
    remake_vectors(true), remake_targets(true), remake_ranges(true)
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
        case 3:
            norm_mode = P_NORM;
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

void LNLP::set_epsilon(const double new_epsilon)
{
    epsilon = new_epsilon;
    return;
}

void LNLP::set_params(const size_t new_E, const size_t new_tau, const int new_tp, const size_t new_nn)
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

void LNLP::set_p(const double new_p)
{
    p = new_p;
    return;
}

void LNLP::suppress_warnings()
{
    SUPPRESS_WARNINGS = true;
    return;
}

void LNLP::save_smap_coefficients()
{
    SAVE_SMAP_COEFFICIENTS = true;
    return;
}

void LNLP::glm()
{
    GLM = true;
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
    return DataFrame::create( Named("time") = target_time, 
                              Named("obs") = targets, 
                              Named("pred") = predicted, 
                              Named("pred_var") = predicted_var);
}

List LNLP::get_smap_coefficients()
{     
    return(wrap(smap_coefficients));
}

DataFrame LNLP::get_short_output()
{
    vec short_time(which_pred.size(), qnan);
    vec short_obs(which_pred.size(), qnan);
    vec short_pred(which_pred.size(), qnan);
    
    for(size_t i = 0; i < which_pred.size(); ++i)
    {
        short_time[i] = target_time[which_pred[i]];
        short_obs[i] = targets[which_pred[i]];
        short_pred[i] = predicted[which_pred[i]];
    }
    
    return DataFrame::create( Named("time") = short_time, 
                              Named("obs") = short_obs, 
                              Named("pred") = short_pred);
}

DataFrame LNLP::get_stats()
{
    PredStats output = make_stats();
    PredStats const_output = make_const_stats();
    return DataFrame::create( Named("num_pred") = output.num_pred, 
                              Named("rho") = output.rho, 
                              Named("mae") = output.mae, 
                              Named("rmse") = output.rmse,
                              Named("perc") = output.perc, 
                              Named("p_val") = output.p_val, 
                              Named("const_pred_num_pred") = const_output.num_pred, 
                              Named("const_pred_rho") = const_output.rho, 
                              Named("const_pred_mae") = const_output.mae, 
                              Named("const_pred_rmse") = const_output.rmse, 
                              Named("const_pred_perc") = const_output.perc, 
                              Named("const_p_val") = const_output.p_val);
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
    if(tp >= 0)
    {
        targets.assign(time_series.begin()+tp, time_series.end());
        targets.insert(targets.end(), tp, qnan);
        target_time.assign(time.begin()+tp, time.end());
        target_time.insert(target_time.end(), tp, qnan);
    }
    else
    {
        targets.assign(time_series.begin(), time_series.end()+tp);
        targets.insert(targets.begin(), -tp, qnan);
        target_time.assign(time.begin(), time.end()+tp);
        target_time.insert(target_time.begin(), -tp, qnan);
    }
    const_targets = time_series;
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
    .method("set_epsilon", &LNLP::set_epsilon)
    .method("set_params", &LNLP::set_params)
    .method("set_theta", &LNLP::set_theta)
    .method("set_p", &LNLP::set_p)
    .method("suppress_warnings", &LNLP::suppress_warnings)
    .method("save_smap_coefficients", &LNLP::save_smap_coefficients)
    .method("glm", &LNLP::glm)
    .method("run", &LNLP::run)
    .method("get_output", &LNLP::get_output)
    .method("get_smap_coefficients", &LNLP::get_smap_coefficients)
    .method("get_short_output", &LNLP::get_short_output)
    .method("get_stats", &LNLP::get_stats)
    ;
}

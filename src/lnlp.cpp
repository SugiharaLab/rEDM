#include "lnlp.h"

/*** Constructors ***/
LNLP::LNLP(): remake_vectors(true), remake_targets(true), remake_ranges(true)
{
}

void LNLP::set_time(const NumericVector time)
{
    this->time = as<std::vector<double> >(time);
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
    for(size_t i = 0; i < lib.nrow(); ++i)
    {
        lib_ranges[i].first = lib(i,0) - 1; // convert 1-index to 0-index
        lib_ranges[i].second = lib(i,1) - 1;
    }
    return;
}

void LNLP::set_pred(const NumericMatrix pred)
{
    pred_ranges.resize(pred.nrow());
    for(size_t i = 0; i < pred.nrow(); ++i)
    {
        pred_ranges[i].first = pred(i,0) - 1; // convert 1-index to 0-index
        pred_ranges[i].second = pred(i,1) - 1;
    }
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

void LNLP::run()
{
    // check parameters
    prepare_forecast();
    
    // forecast
    forecast(); // forecast code is in forecast_machine
    return;
}

DataFrame LNLP::get_output()
{
    return DataFrame::create( Named("obs") = observed, 
                              Named("pred") = predicted);
}

DataFrame LNLP::get_stats()
{
    size_t num_pred = 0;
    double sum_errors = 0;
    double sum_squared_errors = 0;
    double sum_obs = 0;
    double sum_pred = 0;
    double sum_squared_obs = 0;
    double sum_squared_pred = 0;
    double sum_prod = 0;
    
    for(size_t k = 0; k < num_vectors; ++k)
    {
        if(!isnan(observed[k]) && !isnan(predicted[k]))
        {
            ++ num_pred;
            sum_errors += fabs(observed[k] - predicted[k]);
            sum_squared_errors += (observed[k] - predicted[k]) * (observed[k] - predicted[k]);
            sum_obs += observed[k];
            sum_pred += predicted[k];
            sum_squared_obs += observed[k] * observed[k];
            sum_squared_pred += predicted[k] * predicted[k];
            sum_prod += observed[k] * predicted[k];
        }
    }
    double rho = (sum_prod * num_pred - sum_obs * sum_pred) / 
            sqrt((sum_squared_obs * num_pred - sum_obs * sum_obs) * 
                 (sum_squared_pred * num_pred - sum_pred * sum_pred));

    return DataFrame::create( Named("num_pred") = num_pred, 
                              Named("rho") = rho, 
                              Named("mae") = sum_errors / num_pred, 
                              Named("rmse") = sqrt(sum_squared_errors / num_pred));
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
        set_indices_from_range(lib_indices, lib_ranges);
        set_indices_from_range(pred_indices, pred_ranges);

        check_cross_validation();

        which_lib = which_indices_true(lib_indices);
        which_pred = which_indices_true(pred_indices);
        
        remake_ranges = false;
    }
    
    compute_distances();
    sort_neighbors();
    return;
}

void LNLP::make_vectors()
{
    data_vectors.assign(num_vectors, vector<double>(E, qnan));

    // beginning of lagged vectors cannot lag before start of time series
    for(int i = 0; i < (E-1)*tau; ++i)
        for(int j = 0; j < E; ++j)
            if((i - j*tau) >= 0)
                data_vectors[i][j] = time_series[i - j * tau];
    
    // remaining lagged vectors
    for(size_t i = (E-1)*tau; i < num_vectors; ++i)
        for(int j = 0; j < E; ++j)
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

void LNLP::set_indices_from_range(vector<bool>& indices, const vector<time_range>& range)
{
    size_t start_of_range, end_of_range;
    indices.assign(num_vectors, false); // initialize indices
    for(auto& range_iter: range)
    {
        start_of_range = range_iter.first + (E-1)*tau;
        end_of_range = range_iter.second - max(0, tp);
        if(end_of_range >= num_vectors) // check end of range
        {
            cerr << "end_of_range = " << end_of_range << ", but num_vectors = " << num_vectors << "\n";
            LOG_WARNING("end of time_range was greater than the number of vectors; corrected");
            end_of_range = num_vectors-1;
        }
        
        for(size_t j = start_of_range; j <= end_of_range; ++j)
            indices[j] = true;
	}
	return;
}

void LNLP::check_cross_validation()
{
    CROSS_VALIDATION = true;
    if (exclusion_radius < 0) // if exclusion_radius is set, always do cross_validation
    {
        for (size_t i = 0; i < num_vectors; ++i) // see if all lib indices == pred_indices
        {
            if(lib_indices[i] != pred_indices[i])
            {
                CROSS_VALIDATION = false;
                break;
            }
        }
    }
    
    if(!CROSS_VALIDATION) // some difference -> resolve any equal cases
    {
        for(size_t i = 0; i < num_vectors; ++i)
        {
            if(lib_indices[i] && pred_indices[i])
            {
                lib_indices[i] = true;
                pred_indices[i] = false;
            }
        }
    }
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
    .method("run", &LNLP::run)
    .method("get_output", &LNLP::get_output)
    .method("get_stats", &LNLP::get_stats)
    ;
}
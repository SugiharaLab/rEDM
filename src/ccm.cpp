#include "ccm.h"

/*** Constructors ***/
CCM::CCM(): remake_vectors(true), remake_targets(true), remake_ranges(true)
{
    pred_mode = SIMPLEX;
}

void CCM::set_time(const NumericVector new_time)
{
    time = as<std::vector<double> >(new_time);
    return;
}

void CCM::set_block(const NumericMatrix new_block)
{
    size_t num_cols = size_t(new_block.ncol());
    block.resize(num_cols);
    num_vectors = size_t(new_block.nrow());
    for(size_t i = 0; i < num_cols; ++i)
    {
        block[i].resize(num_vectors);
        for(size_t j = 0; j < num_vectors; ++j)
            block[i][j] = new_block(j,i);
    }
    init_distances();
    return;
}

void CCM::set_norm_type(const int norm_type)
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

void CCM::set_lib(const NumericMatrix lib)
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

void CCM::set_pred(const NumericMatrix pred)
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

void CCM::set_lib_sizes(const NumericVector new_lib_sizes)
{
    lib_sizes = as<std::vector<size_t> >(new_lib_sizes);
    return;
}

void CCM::set_exclusion_radius(const double new_exclusion_radius)
{
    exclusion_radius = new_exclusion_radius;
    if(exclusion_radius >= 0)
        CROSS_VALIDATION = true;
    return;
}

void CCM::set_lib_column(const size_t new_lib_col)
{
    lib_col = new_lib_col;
    remake_vectors = true;
    return;
}

void CCM::set_target_column(const size_t new_target)
{
    target = new_target;
    remake_targets = true;
    return;
}

void CCM::set_params(const size_t new_E, const size_t new_tau, const int new_tp, 
                    const size_t new_nn, const bool new_random_libs, 
                    const size_t new_num_samples, const bool new_replace)
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
    random_libs = new_random_libs;
    num_samples = new_num_samples;
    replace = new_replace;
    
    return;
}

void CCM::suppress_warnings()
{
    SUPPRESS_WARNINGS = true;
    return;
}

void CCM::run()
{
    prepare_forecast(); // check parameters
    
    // setup data structures and compute maximum lib size
    stats.clear();
    predicted_lib_sizes.clear();
    size_t max_lib_size = full_lib.size();
    std::mt19937 rng(42); // init mersenne twister with 42 as seed
    std::uniform_int_distribution<uint32_t> lib_sampler(0, max_lib_size); 
    vector<int> idx;
 
    for(auto lib_size: lib_sizes)
    {
        if(lib_size == max_lib_size && (!random_libs || !replace))
        // no possible lib variation if using all vectors and
        // [no random libs OR (random_libs and sampling without replacement)]
        {
            which_lib = full_lib; // use all lib vectors
            forecast();
            stats.push_back(compute_stats());
            predicted_lib_sizes.push_back(lib_size);
        }
        else if(random_libs)
        {
            for(size_t k = 0; k < num_samples; ++k)
            {
                if(replace)
                {
                    which_lib.resize(lib_size);
                    for(auto& lib: which_lib)
                    {
                        lib = full_lib[lib_sampler(rng)];
                    }
                }
                else
                {
                    A:LSKJD:ALKSJ
                }
                forecast();
                stats.push_back(compute_stats());
                predicted_lib_sizes.push_back(lib_size);
            }
        }
        else
        // no random libs and using contiguous segments
        {
            for(size_t k = 0; k < max_lib_size; ++k)
            {
            // setup libs as contiguous segments ***
                forecast();
                stats.push_back(compute_stats());
                predicted_lib_sizes.push_back(lib_size);
            }
        }
    }
    return;
}

DataFrame CCM::get_output()
{
        PredStats output = compute_stats();
    return DataFrame::create( Named("num_pred") = output.num_pred, 
                              Named("rho") = output.rho, 
                              Named("mae") = output.mae, 
                              Named("rmse") = output.rmse );


    return DataFrame::create( Named("time") = time, 
                              Named("obs") = observed, 
                              Named("pred") = predicted, 
                              Named("pred_var") = predicted_var);
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void CCM::prepare_forecast()
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

        full_lib = which_indices_true(lib_indices);
        which_pred = which_indices_true(pred_indices);
        
        remake_ranges = false;
    }
    
    compute_distances();
    sort_neighbors();
    
    return;
}

void CCM::make_vectors()
{
    auto time_series = block[lib_col];
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

void CCM::make_targets()
{
    if((target < 1) || (target-1 >= block.size()))
    {
        throw std::domain_error("invalid target column");
    }
    
    observed.clear();
    if(tp >= 0)
    {
        observed.assign(block[target-1].begin()+tp, block[target-1].end());
        observed.insert(observed.end(), tp, qnan);
    }
    else
    {
        observed.assign(block[target-1].begin(), block[target-1].end()+tp);
        observed.insert(observed.begin(), -tp, qnan);
    }
    remake_targets = false;
    return;
}

RCPP_MODULE(ccm_module)
{
    class_<CCM>("CCM")
    
    .constructor()
    
    .method("set_time", &CCM::set_time)
    .method("set_block", &CCM::set_block)
    .method("set_norm_type", &CCM::set_norm_type)
    .method("set_lib", &CCM::set_lib)
    .method("set_pred", &CCM::set_pred)
    .method("set_lib_sizes", &CCM::set_lib_sizes)
    .method("set_exclusion_radius", &CCM::set_exclusion_radius)
    .method("set_lib_column", &CCM::set_lib_column)
    .method("set_target_column", &CCM::set_target_column)
    .method("set_params", &CCM::set_params)
    .method("suppress_warnings", &CCM::suppress_warnings)
    .method("run", &CCM::run)
    .method("get_output", &CCM::get_output)
    ;
}
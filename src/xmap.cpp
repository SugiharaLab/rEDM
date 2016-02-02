#include "xmap.h"

/*** Constructors ***/
Xmap::Xmap(): remake_vectors(true), remake_targets(true), remake_ranges(true)
{
    pred_mode = SIMPLEX;
}

void Xmap::set_time(const NumericVector new_time)
{
    time = as<std::vector<double> >(new_time);
    return;
}

void Xmap::set_block(const NumericMatrix new_block)
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

void Xmap::set_norm_type(const int norm_type)
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

void Xmap::set_lib(const NumericMatrix lib)
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

void Xmap::set_pred(const NumericMatrix pred)
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

void Xmap::set_lib_sizes(const NumericVector new_lib_sizes)
{
    lib_sizes = as<std::vector<size_t> >(new_lib_sizes);
    return;
}

void Xmap::set_exclusion_radius(const double new_exclusion_radius)
{
    exclusion_radius = new_exclusion_radius;
    if(exclusion_radius >= 0)
        CROSS_VALIDATION = true;
    return;
}

void Xmap::set_epsilon(const double new_epsilon)
{
    epsilon = new_epsilon;
    return;
}

void Xmap::set_lib_column(const size_t new_lib_col)
{
    lib_col = new_lib_col;
    remake_vectors = true;
    return;
}

void Xmap::set_target_columns(const NumericVector new_targets)
{
    target_cols = as<std::vector<size_t> >(new_targets);
    remake_targets = true;
    return;
}

void Xmap::set_params(const size_t new_E, const size_t new_tau, const int new_tp, 
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

void Xmap::suppress_warnings()
{
    SUPPRESS_WARNINGS = true;
    return;
}

void Xmap::run()
{
    prepare_forecast(); // check parameters
    // setup data structures and compute maximum lib size
    predicted_stats.clear();
    predicted_lib_sizes.clear();
    predicted_targets.clear();
    std::vector<size_t> full_lib = which_lib;
    size_t max_lib_size = full_lib.size();
    std::mt19937 rng(42); // init mersenne twister with 42 as seed
    // need to update with true random seed
    std::uniform_int_distribution<uint32_t> lib_sampler(0, max_lib_size-1);
    std::uniform_real_distribution<double> unif_01(0, 1);
    std::vector<int> idx;
 
    size_t m;
    size_t t;

    for(auto lib_size: lib_sizes)
    {
        if(lib_size >= max_lib_size && (!random_libs || !replace))
        // no possible lib variation if using all vectors and
        // [no random libs OR (random_libs and sampling without replacement)]
        {
            which_lib = full_lib; // use all lib vectors
            forecast();
            for(size_t target_idx = 0; target_idx < num_targets; ++target_idx)
            {
                predicted_stats.push_back(make_stats(target_idx));
                predicted_lib_sizes.push_back(lib_size);
                predicted_targets.push_back(target_cols[target_idx]);
            }
            break;
        }
        else if(random_libs)
        {
            which_lib.resize(lib_size, 0);
            for(size_t k = 0; k < num_samples; ++k)
            {
                if(replace)
                {
                    for(auto& lib: which_lib)
                    {
                        lib = full_lib[lib_sampler(rng)];
                    }
                }
                else
                {
                    // sample without replacement (algorithm from Knuth)
                    m = 0;
                    t = 0;
                    while(m < lib_size)
                    {
                        if(double(max_lib_size - t) * unif_01(rng) >= double(lib_size - m))
                        {
                            ++t;
                        }
                        else
                        {
                            which_lib[m] = full_lib[t];
                            ++t; ++m;
                        }
                    }
                }
                forecast();
                for(size_t target_idx = 0; target_idx < num_targets; ++target_idx)
                {
                    predicted_stats.push_back(make_stats(target_idx));
                    predicted_lib_sizes.push_back(lib_size);
                    predicted_targets.push_back(target_cols[target_idx]);
                }
            }
        }
        else
        // no random libs and using contiguous segments
        {
            for(size_t k = 0; k < max_lib_size; ++k)
            {
                if((k + lib_size) > max_lib_size) // need to loop around
                {
                    which_lib.assign(full_lib.begin()+k, full_lib.end()); // k to end
                    which_lib.insert(which_lib.begin(), 
                                     full_lib.begin(), 
                                     full_lib.begin() + lib_size - (max_lib_size-k));
                }
                else
                {
                    which_lib.assign(full_lib.begin()+k, full_lib.begin()+k+lib_size);
                }
                forecast();
                for(size_t target_idx = 0; target_idx < num_targets; ++target_idx)
                {
                    predicted_stats.push_back(make_stats(target_idx));
                    predicted_lib_sizes.push_back(lib_size);
                    predicted_targets.push_back(target_cols[target_idx]);
                }
            }
        }
    }
    which_lib = full_lib;
    return;
}

DataFrame Xmap::get_output()
{
    std::vector<size_t> num_pred;
    std::vector<double> rho;
    std::vector<double> mae;
    std::vector<double> rmse;

    for(auto& stats: predicted_stats)
    {
        num_pred.push_back(stats.num_pred);
        rho.push_back(stats.rho);
        mae.push_back(stats.mae);
        rmse.push_back(stats.rmse);
    }

    return DataFrame::create( Named("target_column") = predicted_targets, 
                              Named("lib_size") = predicted_lib_sizes, 
                              Named("num_pred") = num_pred, 
                              Named("rho") = rho, 
                              Named("mae") = mae, 
                              Named("rmse") = rmse);
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void Xmap::prepare_forecast()
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
    //sort_neighbors();
    
    return;
}

void Xmap::make_vectors()
{
    if((lib_col < 1) || (lib_col-1 >= block.size()))
    {
        throw std::domain_error("invalid target column");
    }
    
    auto time_series = block[lib_col-1];
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

void Xmap::make_targets()
{
    targets.clear();
    const_targets.clear();
    vec curr_target;
    
    if(tp >= 0)
    {
        target_time.assign(time.begin()+tp, time.end());
        target_time.insert(target_time.end(), tp, qnan);
    }
    else
    {
        target_time.assign(time.begin(), time.end()+tp);
        target_time.insert(target_time.begin(), -tp, qnan);
    }
    
    for(auto& curr_col: target_cols)
    {
        if((curr_col < 1) || (curr_col-1 >= block.size()))
        {
            throw std::domain_error("invalid target column");
        }
        curr_target.clear();
        if(tp >= 0)
        {
            curr_target.assign(block[curr_col-1].begin()+tp, block[curr_col-1].end());
            curr_target.insert(curr_target.end(), tp, qnan);
        }
        else
        {
            curr_target.assign(block[curr_col-1].begin(), block[curr_col-1].end()+tp);
            curr_target.insert(curr_target.begin(), -tp, qnan);
        }
        targets.push_back(curr_target);
        const_targets.push_back(block[curr_col-1]);
    }
    remake_targets = false;
    return;
}

RCPP_MODULE(xmap_module)
{
    class_<Xmap>("Xmap")
    
    .constructor()
    
    .method("set_time", &Xmap::set_time)
    .method("set_block", &Xmap::set_block)
    .method("set_norm_type", &Xmap::set_norm_type)
    .method("set_lib", &Xmap::set_lib)
    .method("set_pred", &Xmap::set_pred)
    .method("set_lib_sizes", &Xmap::set_lib_sizes)
    .method("set_exclusion_radius", &Xmap::set_exclusion_radius)
    .method("set_epsilon", &Xmap::set_epsilon)
    .method("set_lib_column", &Xmap::set_lib_column)
    .method("set_target_columns", &Xmap::set_target_columns)
    .method("set_params", &Xmap::set_params)
    .method("suppress_warnings", &Xmap::suppress_warnings)
    .method("run", &Xmap::run)
    .method("get_output", &Xmap::get_output)
    ;
}
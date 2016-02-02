#include "block_lnlp.h"

/*** Constructors ***/
BlockLNLP::BlockLNLP(): remake_vectors(true), remake_targets(true), remake_ranges(true)
{
}

void BlockLNLP::set_time(const NumericVector new_time)
{
    time = as<std::vector<double> >(new_time);
    return;
}

void BlockLNLP::set_block(const NumericMatrix new_block)
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

void BlockLNLP::set_norm_type(const int norm_type)
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

void BlockLNLP::set_pred_type(const int pred_type)
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

void BlockLNLP::set_lib(const NumericMatrix lib)
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

void BlockLNLP::set_pred(const NumericMatrix pred)
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

void BlockLNLP::set_exclusion_radius(const double new_exclusion_radius)
{
    exclusion_radius = new_exclusion_radius;
    if(exclusion_radius >= 0)
        CROSS_VALIDATION = true;
    return;
}

void BlockLNLP::set_epsilon(const double new_epsilon)
{
    epsilon = new_epsilon;
    return;
}

void BlockLNLP::set_embedding(const NumericVector new_embedding)
{
    embedding = as<std::vector<size_t> >(new_embedding);
    E = embedding.size();
    remake_vectors = true;
    return;
}

void BlockLNLP::set_target_columns(const NumericVector new_targets)
{
    target_cols = as<std::vector<size_t> >(new_targets);
    remake_targets = true;
    return;
}

void BlockLNLP::set_params(const int new_tp, const size_t new_nn)
{
    if(tp != new_tp)
        remake_targets = true;
    if(remake_vectors || remake_targets)
        remake_ranges = true;
    
    tp = new_tp;
    nn = new_nn;    
    return;
}

void BlockLNLP::set_theta(const double new_theta)
{
    theta = new_theta;
    return;
}

void BlockLNLP::suppress_warnings()
{
    SUPPRESS_WARNINGS = true;
    return;
}

void BlockLNLP::save_smap_coefficients()
{
    SAVE_SMAP_COEFFICIENTS = true;
    return;
}

void BlockLNLP::run()
{
    prepare_forecast(); // check parameters
    forecast(); // forecast code is in forecast_machine
    return;
}

DataFrame BlockLNLP::get_output()
{
    return DataFrame::create( Named("time") = target_time, 
                              Named("obs") = targets, 
                              Named("pred") = predicted, 
                              Named("pred_var") = predicted_var);
}

List BlockLNLP::get_smap_coefficients()
{     
    return(wrap(smap_coefficients));
}

DataFrame BlockLNLP::get_short_output(size_t target_idx)
{
    vec short_time(which_pred.size(), qnan);
    vec short_obs(which_pred.size(), qnan);
    vec short_pred(which_pred.size(), qnan);
    
    for(size_t i = 0; i < which_pred.size(); ++i)
    {
        short_time[i] = target_time[which_pred[i]];
        short_obs[i] = targets[target_idx][which_pred[i]];
        short_pred[i] = predicted[target_idx][which_pred[i]];
    }
    
    return DataFrame::create( Named("time") = short_time, 
                              Named("obs") = short_obs, 
                              Named("pred") = short_pred);
}

DataFrame BlockLNLP::get_stats()
{
    std::vector<size_t> target(num_targets, 0);
    std::vector<size_t> num_pred(num_targets, 0);
    std::vector<double> rho(num_targets, 0);
    std::vector<double> mae(num_targets, 0);
    std::vector<double> rmse(num_targets, 0);
    std::vector<double> perc(num_targets, 0);
    std::vector<size_t> const_num_pred(num_targets, 0);
    std::vector<double> const_rho(num_targets, 0);
    std::vector<double> const_mae(num_targets, 0);
    std::vector<double> const_rmse(num_targets, 0);
    std::vector<double> const_perc(num_targets, 0);
    
    // assemble each column
    for(size_t target_idx = 0; target_idx < num_targets; ++target_idx)
    {
        target[target_idx] = target_cols[target_idx];

        PredStats output = make_stats(target_idx);
        num_pred[target_idx] = output.num_pred;
        rho[target_idx] = output.rho;
        mae[target_idx] = output.mae;
        rmse[target_idx] = output.rmse;
        perc[target_idx] = output.perc;
        
        PredStats const_output = make_const_stats(target_idx);
        const_num_pred[target_idx] = const_output.num_pred;
        const_rho[target_idx] = const_output.rho;
        const_mae[target_idx] = const_output.mae;
        const_rmse[target_idx] = const_output.rmse;
        const_perc[target_idx] = const_output.perc;
    }
    
    return DataFrame::create( Named("target") = target, 
                              Named("num_pred") = num_pred, 
                              Named("rho") = rho, 
                              Named("mae") = mae, 
                              Named("rmse") = rmse,
                              Named("perc") = perc, 
                              Named("const_pred_num_pred") = const_num_pred, 
                              Named("const_pred_rho") = const_rho, 
                              Named("const_pred_mae") = const_mae, 
                              Named("const_pred_rmse") = const_rmse, 
                              Named("const_pred_perc") = const_perc);
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void BlockLNLP::prepare_forecast()
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
        set_indices_from_range(lib_indices, lib_ranges, 0, -std::max(0, tp), true);
        set_indices_from_range(pred_indices, pred_ranges, 0, -std::max(0, tp), false);

        check_cross_validation();

        which_lib = which_indices_true(lib_indices);
        which_pred = which_indices_true(pred_indices);
        
        remake_ranges = false;
    }
    
    compute_distances();
    //sort_neighbors();
    
    return;
}

void BlockLNLP::make_vectors()
{
    data_vectors.assign(num_vectors, vec(E, qnan));
    for(size_t i = 0; i < num_vectors; ++i)
    {
        for(size_t j = 0; j < E; ++j)
        {
            data_vectors[i][j] = block[embedding[j]-1][i];
        }
    }

    remake_vectors = false;
    return;
}

void BlockLNLP::make_targets()
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

RCPP_MODULE(block_lnlp_module)
{
    class_<BlockLNLP>("BlockLNLP")
    
    .constructor()
    
    .method("set_time", &BlockLNLP::set_time)
    .method("set_block", &BlockLNLP::set_block)
    .method("set_norm_type", &BlockLNLP::set_norm_type)
    .method("set_pred_type", &BlockLNLP::set_pred_type)
    .method("set_lib", &BlockLNLP::set_lib)
    .method("set_pred", &BlockLNLP::set_pred)
    .method("set_exclusion_radius", &BlockLNLP::set_exclusion_radius)
    .method("set_epsilon", &BlockLNLP::set_epsilon)
    .method("set_embedding", &BlockLNLP::set_embedding)
    .method("set_target_columns", &BlockLNLP::set_target_columns)
    .method("set_params", &BlockLNLP::set_params)
    .method("set_theta", &BlockLNLP::set_theta)
    .method("suppress_warnings", &BlockLNLP::suppress_warnings)
    .method("save_smap_coefficients", &BlockLNLP::save_smap_coefficients)
    .method("run", &BlockLNLP::run)
    .method("get_output", &BlockLNLP::get_output)
    .method("get_smap_coefficients", &BlockLNLP::get_smap_coefficients)
    .method("get_short_output", &BlockLNLP::get_short_output)
    .method("get_stats", &BlockLNLP::get_stats)
    ;
}

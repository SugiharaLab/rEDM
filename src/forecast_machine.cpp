#include "forecast_machine.h"

static const double min_weight = 0.000001;
static const double qnan = std::numeric_limits<double>::quiet_NaN();

ForecastMachine::ForecastMachine()
{
}

void ForecastMachine::forecast()
{    
    switch(pred_mode)
    {
        case SIMPLEX:
            simplex_forecast();
            break;
        case SMAP:
            smap_forecast();
            break;
        default:
            throw std::domain_error("Unknown pred type");
    }
    return;
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void ForecastMachine::compute_distances()
{
    double (*dist_func)(const vec&, const vec&);
    
    // select distance function
    switch(norm_mode)
    {
        case L1_NORM:
            dist_func = &l1_distance_func;
            break;
        case L2_NORM:
            dist_func = &l2_distance_func;
            break;
        default:
            throw std::domain_error("Unknown norm type");
    }
    
    // compute distances
    distances.resize(num_vectors, vector<double>(num_vectors, qnan)); // initialize distance matrix
    
    for(auto& curr_pred: which_pred)
    {
        for(auto& curr_lib: which_lib)
        {
            if(isnan(distances[curr_pred][curr_lib]))
                distances[curr_pred][curr_lib] = dist_func(data_vectors[curr_pred], 
                                                           data_vectors[curr_lib]);
        }
    }
    return;
}

void ForecastMachine::sort_neighbors()
{
    neighbors.resize(num_vectors);
    for(auto& curr_pred: which_pred)
        neighbors[curr_pred] = sort_indices(distances[curr_pred], which_lib);
    return;
}

void ForecastMachine::simplex_forecast()
{
    if(CROSS_VALIDATION)
    {
        for(auto& curr_pred: which_pred)
        {
            make_lib(curr_pred); // adjust library for cross-validation
            simplex_prediction(curr_pred);
        }
    }
    else
    {
        make_lib();
        for(auto& curr_pred: which_pred)
			simplex_prediction(curr_pred);
    }
    return;
}

void ForecastMachine::smap_forecast()
{
    for(auto& curr_pred: which_pred)
    {
        make_lib(curr_pred); // adjust library for cross-validation
        smap_prediction(curr_pred);
    }
    return;
}

void ForecastMachine::make_lib()
{
    
    return;
}

void ForecastMachine::make_lib(const size_t curr_pred)
{
    // go through lib
    // remove lib vectors that are within exclusion
    
    return;
}

void ForecastMachine::simplex_prediction(const size_t curr_pred)
{
    vector<size_t> nearest_neighbors; // cleared by default
        
    // find nearest neighbors
    size_t j = 0;
    double tie_distance;
    while(int(nearest_neighbors.size()) < nn)
    {
        if(valid_lib_indices[neighbors[curr_pred][j]]) // is this vector in library?
        {
            nearest_neighbors.push_back(j);
            ++j;
        }
    }
            
    // check for ties
    double prev_distance = distances[curr_pred][j-1];
    bool done_checking_ties = false;
    int effective_lib_size;
    while(!done_checking_ties && nearest_neighbors.size() < effective_lib_size)
    {
        if(lib_indices[neighbors[curr_pred][j]]) // is this vector in library?
        {
            if(distances[curr_pred][j] == prev_distance) // distance is the same
            {
                tie_distance = prev_distance;
                nearest_neighbors.push_back(j);
            }
            else
                done_checking_ties = true;
        }
        ++j;
    }
    
    // compute weights
    double min_distance = distances[curr_pred][nearest_neighbors[0]];
    vector<double> weights;
    if(min_distance == 0)
    {
        weights.assign(nearest_neighbors.size(), min_weight);
        for(size_t k = 0; k < nearest_neighbors.size(); k++)
        {
            if(distances[curr_pred][nearest_neighbors[k]] == 0)
                weights[k] = 1;
            else
                break;
        }
    }
    else
    {
        weights.assign(nearest_neighbors.size(), 0);
        for(size_t k = 0; k < nearest_neighbors.size(); k++)
        {
            weights[k] = exp(-distances[curr_pred][nearest_neighbors[k]] / min_distance);
            if(weights[k] < min_weight)
                weights[k] = min_weight;
        }
    }
        
    // identify ties and adjust weights
    int num_ties = 0;
    int neighbor_slots_used_for_ties;
    double tie_adj_factor;
    if(nearest_neighbors.size() > nn) // ties exist
    {
        // count ties
        for(int k = 0; k < nearest_neighbors.size(); k++)
        {
            if(distances[curr_pred][nearest_neighbors[k]] == tie_distance)
                num_ties++;
        }
        neighbor_slots_used_for_ties = nn - (nearest_neighbors.size() - num_ties);
        tie_adj_factor = double(neighbor_slots_used_for_ties) / double(num_ties);
        
        // adjust weights
        for(int k = 0; k < nearest_neighbors.size(); k++)
        {
            if(distances[curr_pred][nearest_neighbors[k]] == tie_distance)
                weights[k] *= tie_adj_factor;
        }
    }
        
    // make prediction
    double total_weight = 0;
    predicted[curr_pred] = 0;
    for(int k = 0; k < nearest_neighbors.size(); k++)
    {
        predicted[curr_pred] += weights[k] * target_vals[neighbors[curr_pred][nearest_neighbors[k]]];
        total_weight += weights[k];
    }
    predicted[curr_pred] = predicted[curr_pred] / total_weight;
    return;
}

void ForecastMachine::smap_prediction(const size_t curr_pred)
{
    double pred = 0;
    // compute smap
    // compute predicted value
    
    // save prediction
    predicted[curr_pred] = pred;
    return;
}

bool ForecastMachine::is_vec_valid(const size_t vec_index)
{
    // check dims and target vals
    return true;
}

void ForecastMachine::LOG_WARNING(const char* warning_text)
{
    cerr << "WARNING: " << warning_text << "\n";
}

vector<size_t> which_indices_true(const vector<bool>& indices)
{
    vector<size_t> which;
    int index = 0;
    for(vector<bool>::const_iterator iter = indices.begin();
        iter != indices.end(); ++iter, ++ index)
    {
        if(*iter)
            which.push_back(index);
    }
    return which;
}

double l1_distance_func(const vec& A, const vec& B)
{
    double dist = 0;
	for (vector<double>::const_iterator a_iter = A.begin(), b_iter = B.begin();
         a_iter != A.end(); ++a_iter, ++b_iter)
	{
		dist += fabs(*a_iter - *b_iter);
	}
	return dist;
}

double l2_distance_func(const vec& A, const vec& B)
{
    double dist = 0;
	for (vector<double>::const_iterator a_iter = A.begin(), b_iter = B.begin();
         a_iter != A.end(); ++a_iter, ++b_iter)
	{
		dist += (*a_iter - *b_iter) * (*a_iter - *b_iter);
	}
	return sqrt(dist);
}

vector<size_t> sort_indices(const vector<double>& v, vector<size_t> idx)
{
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}
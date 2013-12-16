#include "forecast_machine.h"

static const double min_weight = 0.000001;
const double ForecastMachine::qnan = std::numeric_limits<double>::quiet_NaN();

ForecastMachine::ForecastMachine(): valid_lib_indices(vector<bool>()), 
lib_indices(vector<bool>()), pred_indices(vector<bool>()), 
which_lib(vector<size_t>()), which_pred(vector<size_t>()), 
time(vector<double>()), data_vectors(vector<vec>()), 
observed(vector<double>()), predicted(vector<double>()), 
num_vectors(0), 
distances(vector<vector<double> >()), neighbors(vector<vector<size_t> >()), 
CROSS_VALIDATION(false), pred_mode(SIMPLEX), norm_mode(L2_NORM), 
nn(0), exclusion_radius(-1), 
lib_ranges(vector<time_range>()), pred_ranges(vector<time_range>())
{
}

void ForecastMachine::init_distances()
{
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
    
    // initialize distance matrix
    distances.assign(num_vectors, vector<double>(num_vectors, qnan));
    return;
}

void ForecastMachine::compute_distances()
{
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

void ForecastMachine::forecast()
{
    predicted.assign(num_vectors, qnan); // initialize predictions
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

bool ForecastMachine::is_vec_valid(const size_t vec_index)
{
    // check data vector
    for(auto& val: data_vectors[vec_index])
        if(isnan(val)) return false;
    
    // check target value
    if(isnan(observed[vec_index])) return false;

    // if all is good, then:
    return true;
}

void ForecastMachine::LOG_WARNING(const char* warning_text)
{
    cerr << "WARNING: " << warning_text << "\n";
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void ForecastMachine::simplex_forecast()
{
    if(CROSS_VALIDATION)
    {
        for(auto& curr_pred: which_pred)
        {
            adjust_lib(curr_pred); // adjust library for cross-validation
            simplex_prediction(curr_pred);
        }
    }
    else
    {
        adjust_lib();
        for(auto& curr_pred: which_pred)
			simplex_prediction(curr_pred);
    }
    return;
}

void ForecastMachine::smap_forecast()
{
    if(CROSS_VALIDATION)
    {
        for(auto& curr_pred: which_pred)
        {
            adjust_lib(curr_pred); // adjust library for cross-validation
            smap_prediction(curr_pred);
        }
    }
    else
    {
        adjust_lib();
        for(auto& curr_pred: which_pred)
    		smap_prediction(curr_pred);
    }
    return;
}

void ForecastMachine::simplex_prediction(const size_t curr_pred)
{
    vector<size_t> nearest_neighbors(nn, 0); // cleared by default
    
    // find nearest neighbors
    size_t i;
    size_t j = 0;
    for(i = 0; i < neighbors[curr_pred].size(); ++i)
    {
        if(valid_lib_indices[neighbors[curr_pred][i]])
        {
            nearest_neighbors[j] = neighbors[curr_pred][i];
            ++j;
            if(j >= nn)
                break;
        }
    }
    double tie_distance = distances[curr_pred][nearest_neighbors.back()];
    
    // check for ties
    for(i; i < neighbors[curr_pred].size(); ++i)
    {
        if(distances[curr_pred][i] > tie_distance) // distance is bigger
            break;
        if(valid_lib_indices[neighbors[curr_pred][i]]) // valid lib
            nearest_neighbors.push_back(j); // add to nearest neighbors
    }
    
    // compute weights
    double min_distance = distances[curr_pred][nearest_neighbors[0]];
    size_t effective_nn = nearest_neighbors.size();
    vector<double> weights(effective_nn, min_weight);
    if(min_distance == 0)
    {
        for(size_t k = 0; k < effective_nn; ++k)
        {
            if(distances[curr_pred][nearest_neighbors[k]] == min_distance)
                weights[k] = 1;
            else
                break;
        }
    }
    else
    {
        for(size_t k = 0; k < effective_nn; ++k)
        {
            weights[k] = max(exp(-distances[curr_pred][nearest_neighbors[k]] / min_distance), 
                             min_weight);
        }
    }
        
    // identify ties and adjust weights
    if(effective_nn > nn) // ties exist
    {
        // count ties
        int num_ties = 0;
        for(auto& neighbor_index: nearest_neighbors)
            if(distances[curr_pred][neighbor_index] == tie_distance)
                num_ties++;

        int neighbor_slots_used_for_ties = nn - (effective_nn - num_ties);
        double tie_adj_factor = double(neighbor_slots_used_for_ties) / double(num_ties);
        
        // adjust weights
        for(size_t k = 0; k < nearest_neighbors.size(); ++k)
            if(distances[curr_pred][nearest_neighbors[k]] == tie_distance)
                weights[k] *= tie_adj_factor;
    }
    
    // make prediction
    double total_weight = accumulate(weights.begin(), weights.end(), 0.0);
    predicted[curr_pred] = 0;    
    for(size_t k = 0; k < effective_nn; ++k)
        predicted[curr_pred] += weights[k] * observed[nearest_neighbors[k]];
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

void ForecastMachine::adjust_lib()
{
    valid_lib_indices = lib_indices;
    return;
}

void ForecastMachine::adjust_lib(const size_t curr_pred)
{
    valid_lib_indices = lib_indices;
    
    // go through lib and remove lib vectors that are within exclusion radius
    if(exclusion_radius >= 0)
    {
        double start_time = time[curr_pred] - exclusion_radius;
        double end_time = time[curr_pred] + exclusion_radius;
        for(size_t i = 0; i < num_vectors; ++i)
        {
            if(lib_indices[i] && time[i] >= start_time && time[i] <= end_time)
                valid_lib_indices[i] = false;
        }
    }
    
    valid_lib_indices[curr_pred] = false; // always remove vector that we are predicting
    return;
}

vector<size_t> which_indices_true(const vector<bool>& indices)
{
    vector<size_t> which;
    int index = 0;
    for(vector<bool>::const_iterator iter = indices.begin();
        iter != indices.end(); ++iter, ++index)
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
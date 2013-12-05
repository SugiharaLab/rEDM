#include "forecast_machine.h"

ForecastMachine::ForecastMachine()
{
}

void ForecastMachine::setup_lib_and_pred()
{
    return;
}

void ForecastMachine::forecast()
{
    setup_lib_and_pred();
    compute_distances();
    
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

void ForecastMachine::set_indices_from_range(vector<bool>& indices,
                                             const vector<time_range>& range)
{
    int start_of_range, end_of_range;
    indices.assign(num_vectors, false); // initialize indices
    for(vector<time_range>::const_iterator curr_range = range.begin();
        curr_range != range.end(); ++curr_range)
    {
        start_of_range = curr_range->first-1;
        end_of_range = curr_range->last-1;
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
        
		for(int vec_index = start_of_range; vec_index <= end_of_range; ++vec_index) // loop through time_range
		{
            indices[vec_index] = is_vec_valid(vec_index);
		}
	}
	return;
}

void ForecastMachine::compute_distances()
{
    double (*dist)(const vec&, const vec&);
    
    // select distance function
    switch(norm_mode)
    {
        case L1_NORM:
            dist = &l1_distance_func;
            break;
        case L2_NORM:
            dist = &l2_distance_func;
            break;
        default:
            throw std::domain_error("Unknown norm type");
    }
    
    // compute distances
    
    return;
}

void ForecastMachine::simplex_forecast()
{
    vector<int> which_preds = which_indices_true(pred_indices);
    
    if(CROSS_VALIDATION)
    {
        for(vector<int>::iterator curr_pred = which_preds.begin();
            curr_pred != which_preds.end(); ++ curr_pred)
        {
            make_lib(*curr_pred); // adjust library for cross-validation
            simplex_prediction(*curr_pred);
        }
    }
    else
    {
        make_lib();
        for(vector<int>::iterator curr_pred = which_preds.begin();
            curr_pred != which_preds.end(); ++ curr_pred)
        {
			simplex_prediction(*curr_pred);
		}
    }
    return;
}

void ForecastMachine::smap_forecast()
{
    vector<int> which_preds = which_indices_true(pred_indices);
    
    for(vector<int>::iterator curr_pred = which_preds.begin();
        curr_pred != which_preds.end(); ++ curr_pred)
    {
        make_lib(*curr_pred); // adjust library for cross-validation
        smap_prediction(*curr_pred);
    }
    return;
}

void ForecastMachine::make_lib()
{
    return;
}

void ForecastMachine::make_lib(const int curr_pred)
{
    return;
}

void ForecastMachine::simplex_prediction(const int curr_pred)
{
    double pred = 0;
    // find nearest neighbors
    // compute weights
    // compute prediction
    
    // save prediction
    predicted[curr_pred] = pred;
    return;
}

void ForecastMachine::smap_prediction(const int curr_pred)
{
    double pred = 0;
    // compute smap
    // compute predicted value
    
    // save prediction
    predicted[curr_pred] = pred;
    return;
}

bool ForecastMachine::is_vec_valid(const int vec_index)
{
    // check dims and target vals
    return true;
}

void ForecastMachine::LOG_WARNING(const char* warning_text)
{
    cerr << "WARNING (ForecastMachine): " << warning_text << "\n";
}

vector<int> which_indices_true(const vector<bool>& indices)
{
    vector<int> which;
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



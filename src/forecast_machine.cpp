#include "forecast_machine.h"

ForecastMachine::ForecastMachine()
{
}

void ForecastMachine::forecast()
{
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
    distances.resize(which_pred.size());
    for(int i = 0; i < int(which_pred.size()); ++i)
    {
        distances[i].resize(which_lib.size());
        for(int j = 0; j < int(which_lib.size()); ++j)
        {
            distances[i][j] = dist_func(data_vectors[which_pred[i]], 
                                        data_vectors[which_lib[j]]);
        }
    }
    
    return;
}

void ForecastMachine::simplex_forecast()
{
    if(CROSS_VALIDATION)
    {
        for(vector<int>::iterator curr_pred = which_pred.begin();
            curr_pred != which_pred.end(); ++ curr_pred)
        {
            make_lib(*curr_pred); // adjust library for cross-validation
            simplex_prediction(*curr_pred);
        }
    }
    else
    {
        make_lib();
        for(vector<int>::iterator curr_pred = which_pred.begin();
            curr_pred != which_pred.end(); ++ curr_pred)
        {
			simplex_prediction(*curr_pred);
		}
    }
    return;
}

void ForecastMachine::smap_forecast()
{
    for(vector<int>::iterator curr_pred = which_pred.begin();
        curr_pred != which_pred.end(); ++ curr_pred)
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
    cerr << "WARNING: " << warning_text << "\n";
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



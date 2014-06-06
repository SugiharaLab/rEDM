#include "forecast_machine.h"

static const double min_weight = 0.000001;
const double ForecastMachine::qnan = std::numeric_limits<double>::quiet_NaN();

ForecastMachine::ForecastMachine():
lib_indices(std::vector<bool>()), pred_indices(std::vector<bool>()),
which_lib(std::vector<size_t>()), which_pred(std::vector<size_t>()),
time(vec()), data_vectors(std::vector<vec>()),
observed(vec()), predicted(vec()),
num_vectors(0),
distances(std::vector<vec>()), neighbors(std::vector<std::vector<size_t> >()),
CROSS_VALIDATION(false), SUPPRESS_WARNINGS(false), pred_mode(SIMPLEX), norm_mode(L2_NORM),
nn(0), exclusion_radius(-1),
lib_ranges(std::vector<time_range>()), pred_ranges(std::vector<time_range>()),
num_threads(1)
{
    num_threads = std::thread::hardware_concurrency();
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
    distances.assign(num_vectors, vec(num_vectors, qnan));
    return;
}

void ForecastMachine::compute_distances()
{
    size_t rows = which_pred.size() / num_threads;
    size_t extra = which_pred.size() % num_threads;
    size_t start = 0;
    size_t end = rows;
    
    std::vector<std::thread> workers;
    for(int t = 1; t <= num_threads; ++t)
    {
        if(t == num_threads)
            end += extra;
        
        // set up calculations to be done
        workers.push_back(std::thread([start, end, this]()
                                 {
                                     size_t curr_pred;
                                     for(size_t i = start; i < end; ++i)
                                     {
                                         curr_pred = which_pred[i];
                                         for(auto& curr_lib: which_lib)
                                         {
                                             if(std::isnan(distances[curr_pred][curr_lib]))
                                                 distances[curr_pred][curr_lib] = dist_func(data_vectors[curr_pred],
                                                                                            data_vectors[curr_lib]);
                                         }
                                     }
                                 }));
        
        // set up rows for next calc
        start = end;
        end = start + rows;
    }
    
    for(auto& tt: workers)
        tt.join();
    return;
}

void ForecastMachine::sort_neighbors()
{
    neighbors.resize(num_vectors);
    for(auto& curr_pred: which_pred)
    {
        neighbors[curr_pred] = sort_indices(distances[curr_pred], which_lib);
    }
    return;
}

std::vector<size_t> ForecastMachine::find_nearest_neighbors(const size_t curr_pred, const std::vector<bool>& valid_lib_indices)
{
    std::vector<size_t> nearest_neighbors;
    
    if(nn < 1)
    {
        for(auto& curr_lib: neighbors[curr_pred])
            if(valid_lib_indices[curr_lib])
                nearest_neighbors.push_back(curr_lib);
        return nearest_neighbors;
    }
    // else
    nearest_neighbors.assign(nn, 0);
    std::vector<size_t>::iterator curr_lib;
    
    // find nearest neighbors
    size_t j = 0;
    for(curr_lib = neighbors[curr_pred].begin(); curr_lib != neighbors[curr_pred].end(); ++curr_lib)
    {
        if(valid_lib_indices[*curr_lib])
        {
            nearest_neighbors[j] = *curr_lib;
            ++j;
            if(j >= nn)
                break;
        }
    }
    double tie_distance = distances[curr_pred][nearest_neighbors.back()];
    
    // check for ties
    for(++curr_lib; curr_lib != neighbors[curr_pred].end(); ++curr_lib)
    {
        if(distances[curr_pred][*curr_lib] > tie_distance) // distance is bigger
            break;
        if(valid_lib_indices[*curr_lib]) // valid lib
            nearest_neighbors.push_back(*curr_lib); // add to nearest neighbors
    }
    
    return nearest_neighbors;
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

void ForecastMachine::set_indices_from_range(std::vector<bool>& indices, const std::vector<time_range>& range,
                                             int start_shift, int end_shift, bool check_target)
{
    size_t start_of_range, end_of_range;
    indices.assign(num_vectors, false); // initialize indices
    for(auto& range_iter: range)
    {
        if(start_shift < 0) // check beginning of range
        {
            LOG_WARNING("adjustment to beginning of ts was < 0; corrected");
            start_shift = 0;
        }
        start_of_range = range_iter.first + start_shift;
        end_of_range = range_iter.second + end_shift;
        if(end_of_range >= num_vectors) // check end of range
        {
            std::string temp = "end_of_range = ";
            temp += std::to_string(end_of_range);
            temp += ", but num_vectors = ";
            temp += std::to_string(num_vectors);
            LOG_WARNING(temp.c_str());
            LOG_WARNING("end of time_range was greater than the number of vectors; corrected");
            end_of_range = num_vectors-1;
        }
        
        for(size_t j = start_of_range; j <= end_of_range; ++j)
        {
            if(is_vec_valid(j) && (!check_target || is_target_valid(j)))
                indices[j] = true;
        }
    }
	return;
}

void ForecastMachine::check_cross_validation()
{
    if (exclusion_radius >= 0) // if exclusion_radius is set, always do cross_validation
    {
        CROSS_VALIDATION = true;
        return;
    }
    // else
    for (size_t i = 0; i < num_vectors; ++i) // see if lib indices overlap with pred_indices
    {
        if(lib_indices[i] && pred_indices[i])
        {
            CROSS_VALIDATION = true; // don't change flags, just enable cross-validation to deal with overlaps
            exclusion_radius = 0; // exclusion_radius < 0 due to above check
            LOG_WARNING("Found overlap between lib and pred. Enabling cross-validation with exclusion radius = 0.");
            return;
        }
    }
    
    /*
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
        // lib == pred   
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
    */
    return;
}

bool ForecastMachine::is_vec_valid(const size_t vec_index)
{
    // check data vector
    for(auto& val: data_vectors[vec_index])
        if(std::isnan(val)) return false;
    
    // if all is good, then:
    return true;
}

bool ForecastMachine::is_target_valid(const size_t vec_index)
{
    // check target value
    if(std::isnan(observed[vec_index])) return false;
    
    // if all is good, then:
    return true;
}

PredStats ForecastMachine::compute_stats()
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
        if(!std::isnan(observed[k]) && !std::isnan(predicted[k]))
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
    
    PredStats output;
    output.num_pred = num_pred;
    output.rho = (sum_prod * num_pred - sum_obs * sum_pred) /
    sqrt((sum_squared_obs * num_pred - sum_obs * sum_obs) *
         (sum_squared_pred * num_pred - sum_pred * sum_pred));
    output.mae = sum_errors / double(num_pred);
    output.rmse = sqrt(sum_squared_errors / double(num_pred));
    
    return output;
}

void ForecastMachine::LOG_WARNING(const char* warning_text)
{
    if(!SUPPRESS_WARNINGS)
        std::cerr << "WARNING: " << warning_text << "\n";
    return;
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void ForecastMachine::simplex_forecast()
{
    size_t rows = which_pred.size() / num_threads;
    size_t extra = which_pred.size() % num_threads;
    size_t start = 0;
    size_t end = rows;
    std::vector<std::thread> workers;
    
    for(int t = 1; t <= num_threads; ++t)
    {
        if(t == num_threads)
            end += extra;
        
        // set up calculations to be done
        workers.push_back(std::thread(&ForecastMachine::simplex_prediction, this, start, end));
        
        // set up rows for next calc
        start = end;
        end = start + rows;
    }
    
    // wait for threads to finish
    for(auto& tt: workers)
        tt.join();
    
    return;
}

void ForecastMachine::smap_forecast()
{
    size_t rows = which_pred.size() / num_threads;
    size_t extra = which_pred.size() % num_threads;
    size_t start = 0;
    size_t end = rows;
    std::vector<std::thread> workers;
    
    for(int t = 1; t <= num_threads; ++t)
    {
        if(t == num_threads)
            end += extra;
        
        // set up calculations to be done
        workers.push_back(std::thread(&ForecastMachine::smap_prediction, this, start, end));
        
        // set up rows for next calc
        start = end;
        end = start + rows;
    }
    
    // wait for threads to finish
    for(auto& tt: workers)
        tt.join();
        
    return;
}

void ForecastMachine::simplex_prediction(const size_t start, const size_t end)
{
    size_t curr_pred, effective_nn, num_ties;
    double min_distance, tie_distance;
    vec weights;
    std::vector<size_t> nearest_neighbors;
    double tie_adj_factor;
    double total_weight;
    
    for(size_t i = start; i < end; ++i)
    {
        curr_pred = which_pred[i];
        
        // find nearest neighbors
        if(CROSS_VALIDATION)
        {
            nearest_neighbors = find_nearest_neighbors(curr_pred, adjust_lib(curr_pred));
        }
        else
        {
            nearest_neighbors = find_nearest_neighbors(curr_pred, lib_indices);
        }
        effective_nn = nearest_neighbors.size();
        
        // compute weights
        min_distance = distances[curr_pred][nearest_neighbors[0]];
        weights.assign(effective_nn, min_weight);
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
                weights[k] = fmax(exp(-distances[curr_pred][nearest_neighbors[k]] / min_distance),
                                 min_weight);
            }
        }
        
        // identify ties and adjust weights
        if(effective_nn > nn) // ties exist
        {
            tie_distance = distances[curr_pred][nearest_neighbors.back()];
            
            // count ties
            num_ties = 0;
            for(auto& neighbor_index: nearest_neighbors)
                if(distances[curr_pred][neighbor_index] == tie_distance)
                    num_ties++;
            
            tie_adj_factor = double(num_ties + nn - effective_nn) / double(num_ties);
            
            // adjust weights
            for(size_t k = 0; k < nearest_neighbors.size(); ++k)
                if(distances[curr_pred][nearest_neighbors[k]] == tie_distance)
                    weights[k] *= tie_adj_factor;
        }
        
        // make prediction
        total_weight = accumulate(weights.begin(), weights.end(), 0.0);
        predicted[curr_pred] = 0;
        for(size_t k = 0; k < effective_nn; ++k)
            predicted[curr_pred] += weights[k] * observed[nearest_neighbors[k]];
        predicted[curr_pred] = predicted[curr_pred] / total_weight;
    }
    return;
}

void ForecastMachine::smap_prediction(const size_t start, const size_t end)
{
    size_t curr_pred, effective_nn, E = data_vectors[0].size();
    double avg_distance;
    vec weights;
    std::vector<size_t> nearest_neighbors;
    MatrixXd A, S_inv;
    VectorXd B, S, x;
    double max_s, pred;
    
    for(size_t i = start; i < end; ++i)
    {
        curr_pred = which_pred[i];
        
        // find nearest neighbors
        if(CROSS_VALIDATION)
        {
            nearest_neighbors = find_nearest_neighbors(curr_pred, adjust_lib(curr_pred));
        }
        else
        {
            nearest_neighbors = find_nearest_neighbors(curr_pred, lib_indices);
        }
        effective_nn = nearest_neighbors.size();
        
        weights.assign(effective_nn, 1.0); // default is for theta = 0
        if(theta > 0.0)
        {
            // compute average distance
            avg_distance = 0;
            for(auto& neighbor: nearest_neighbors)
            {
                avg_distance += distances[curr_pred][neighbor];
            }
            avg_distance /= effective_nn;
            
            // compute weights
            for(size_t i = 0; i < effective_nn; ++i)
                weights[i] = exp(-theta * distances[curr_pred][nearest_neighbors[i]] / avg_distance);
        }
        
        // setup matrices for SVD
        A.resize(effective_nn, E+1);
        B.resize(effective_nn);
        
        for(size_t i = 0; i < effective_nn; ++i)
        {
            B(i) = weights[i] * observed[nearest_neighbors[i]];
            
            for(size_t j = 0; j < E; ++j)
                A(i, j) = weights[i] * data_vectors[nearest_neighbors[i]][j];
            A(i, E) = weights[i];
        }
        
        // perform SVD
        JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
        
        // remove singular values close to 0
        S = svd.singularValues();
        S_inv = MatrixXd::Zero(E+1, E+1);
        max_s = S(0) * 1e-5;
        for(size_t j = 0; j <= E; ++j)
        {
            if(S(j) >= max_s)
                S_inv(j, j) = 1/S(j);
        }
        
        // perform back-substitution to solve
        x = svd.matrixV() * S_inv * svd.matrixU().transpose() * B;
        
        pred = 0;
        for(size_t j = 0; j < E; ++j)
            pred += x(j) * data_vectors[curr_pred][j];
        pred += x(E);
        
        // save prediction
        predicted[curr_pred] = pred;
    }
    return;
}

std::vector<bool> ForecastMachine::adjust_lib(const size_t curr_pred)
{
    std::vector<bool> valid_lib_indices = lib_indices;
    
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
    return valid_lib_indices;
}

std::vector<size_t> which_indices_true(const std::vector<bool>& indices)
{
    std::vector<size_t> which;
    int index = 0;
    for(std::vector<bool>::const_iterator iter = indices.begin();
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
	for (auto a_iter = A.begin(), b_iter = B.begin();
         a_iter != A.end(); ++a_iter, ++b_iter)
	{
		dist += fabs(*a_iter - *b_iter);
	}
	return dist;
}

double l2_distance_func(const vec& A, const vec& B)
{
    double dist = 0;
	for (auto a_iter = A.begin(), b_iter = B.begin();
         a_iter != A.end(); ++a_iter, ++b_iter)
	{
		dist += (*a_iter - *b_iter) * (*a_iter - *b_iter);
	}
	return sqrt(dist);
}

std::vector<size_t> sort_indices(const vec& v, std::vector<size_t> idx)
{
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

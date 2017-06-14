#include "forecast_machine.h"

static const double min_weight = 0.000001;
const double ForecastMachine::qnan = std::numeric_limits<double>::quiet_NaN();

ForecastMachine::ForecastMachine():
lib_indices(std::vector<bool>()), pred_indices(std::vector<bool>()),
which_lib(std::vector<size_t>()), which_pred(std::vector<size_t>()),
time(vec()), data_vectors(std::vector<vec>()), smap_coefficients(std::vector<vec>()),
targets(vec()), predicted(vec()), predicted_var(vec()),
const_targets(vec()), const_predicted(vec()),
num_vectors(0), distances(std::vector<vec>()),
CROSS_VALIDATION(false), SUPPRESS_WARNINGS(false), SAVE_SMAP_COEFFICIENTS(false),
pred_mode(SIMPLEX), norm_mode(L2_NORM),
nn(0), exclusion_radius(-1), epsilon(-1), p(0.5),
lib_ranges(std::vector<time_range>()), pred_ranges(std::vector<time_range>())
{
    //num_threads = std::thread::hardware_concurrency();
}

/*
void ForecastMachine::debug_print_vectors()
{
    size_t i = 0;
    for(auto& v: data_vectors)
    {
        std::cerr << i << ": <";
        for(auto& i: v)
        {
            std::cerr << i << " ";
        }
        std::cerr << "> --> " << targets[i] << "\n";
        ++i;
    }
    return;
}

void ForecastMachine::debug_print_lib_and_pred()
{
    size_t i = 0;
    for(auto l: lib_indices)
    {
        std::cerr << i << ": " << "lib = " << l << ", pred = " << pred_indices[i] << "\n";
        ++i;
    }
    return;
}
*/

void ForecastMachine::init_distances()
{
    // select distance function
    switch(norm_mode)
    {
        case L1_NORM:
            //dist_func = &l1_distance_func;
            dist_func = [](const vec& A, const vec& B)
            {
                double dist = 0;
                for (auto a_iter = A.begin(), b_iter = B.begin();
                     a_iter != A.end(); ++a_iter, ++b_iter)
                {
                    dist += fabs(*a_iter - *b_iter);
                }
                return dist;
            };
            break;
        case L2_NORM:
            //dist_func = &l2_distance_func;
            dist_func = [](const vec& A, const vec& B)
            {
                double dist = 0;
                for (auto a_iter = A.begin(), b_iter = B.begin();
                     a_iter != A.end(); ++a_iter, ++b_iter)
                {
                    dist += (*a_iter - *b_iter) * (*a_iter - *b_iter);
                }
                return sqrt(dist);
            };
            break;
        case P_NORM:
            dist_func = [&](const vec& A, const vec& B){
                            double dist = 0;
                            for (auto a_iter = A.begin(), b_iter = B.begin();
                                 a_iter != A.end(); ++a_iter, ++b_iter)
                            {
                                dist += pow(fabs(*a_iter - *b_iter), p);
                            }
                            return pow(dist, 1/p);
                        };
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
    /*
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
    */
    for(auto& curr_pred: which_pred)
    {
        for(auto& curr_lib: which_lib)
        {
            if(std::isnan(distances[curr_pred][curr_lib]))
            {
                distances[curr_pred][curr_lib] = dist_func(data_vectors[curr_pred],
                                                            data_vectors[curr_lib]);
                distances[curr_lib][curr_pred] = distances[curr_pred][curr_lib];
            }
        }
    }
    /*
                                }));

        // set up rows for next calc
        start = end;
        end = start + rows;
    }

    for(auto& tt: workers)
        tt.join();
    */
    return;
}

std::vector<size_t> ForecastMachine::find_nearest_neighbors(const vec& dist)
{
    if(nn < 1)
    {
        return sort_indices(dist, which_lib);
    }
    // else
    std::vector<size_t> neighbors;
    std::vector<size_t> nearest_neighbors;
    double curr_distance;

    if(nn > log(double(which_lib.size())))
    {
        neighbors = sort_indices(dist, which_lib);
        std::vector<size_t>::iterator curr_lib;

        // find nearest neighbors
        for(curr_lib = neighbors.begin(); curr_lib != neighbors.end(); ++curr_lib)
        {
            nearest_neighbors.push_back(*curr_lib);
            if(nearest_neighbors.size() >= nn)
                break;
        }
        if(curr_lib == neighbors.end())
            return nearest_neighbors;

        double tie_distance = dist[nearest_neighbors.back()];

        // check for ties
        for(++curr_lib; curr_lib != neighbors.end(); ++curr_lib)
        {
            if(dist[*curr_lib] > tie_distance) // distance is bigger
                break;
            nearest_neighbors.push_back(*curr_lib); // add to nearest neighbors
        }
    }
    else
    {
        size_t i;
        nearest_neighbors.push_back(which_lib[0]);
        for(auto curr_lib: which_lib)
        {
            curr_distance = dist[curr_lib];
            if(curr_distance <= dist[nearest_neighbors.back()])
            {
                i = nearest_neighbors.size();
                while((i > 0) && (curr_distance < dist[nearest_neighbors[i-1]]))
                {
                    i--;
                }
                nearest_neighbors.insert(nearest_neighbors.begin()+i, curr_lib);

                if((nearest_neighbors.size() > nn) &&
                   (dist[nearest_neighbors[nn-1]] < dist[nearest_neighbors.back()]))
                {
                    nearest_neighbors.pop_back();
                }
            }
        }
    }

    // filter for max_distance
    if(epsilon >= 0)
    {
        for(auto neighbor_iter = nearest_neighbors.begin(); neighbor_iter != nearest_neighbors.end(); ++neighbor_iter)
        {
            if(dist[*neighbor_iter] > epsilon)
            {
                nearest_neighbors.erase(neighbor_iter, nearest_neighbors.end());
                break;
            }
        }
    }

    return nearest_neighbors;
}

void ForecastMachine::forecast()
{
    predicted.assign(num_vectors, qnan); // initialize predictions
    const_predicted.assign(num_vectors, qnan);
    predicted_var.assign(num_vectors, qnan);
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
            std::ostringstream temp;
            temp << "end_of_range = ";
            temp << end_of_range;
            temp << ", but num_vectors = ";
            temp << num_vectors;
            std::string temp_str = temp.str();
            LOG_WARNING(temp_str.c_str());
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
    if(std::isnan(targets[vec_index])) return false;

    // if all is good, then:
    return true;
}

PredStats ForecastMachine::make_stats()
{
    return compute_stats_internal(targets, predicted);
}

PredStats ForecastMachine::make_const_stats()
{
    return compute_stats_internal(targets, const_predicted);
}

void ForecastMachine::LOG_WARNING(const char* warning_text)
{
    if(!SUPPRESS_WARNINGS)
        Rcpp::warning(warning_text);
    return;
}

// *** PRIVATE METHODS FOR INTERNAL USE ONLY *** //

void ForecastMachine::simplex_forecast()
{
    /*
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
    */
    simplex_prediction(0, which_pred.size());
    const_prediction(0, which_pred.size());
    return;
}

void ForecastMachine::smap_forecast()
{
    /*
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
    */
    if(SAVE_SMAP_COEFFICIENTS)
    {
        smap_coefficients.assign(num_vectors, vec(data_vectors[0].size()+1, qnan));
    }
    smap_prediction(0, which_pred.size());
    const_prediction(0, which_pred.size());
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
    std::vector<size_t> temp_lib;
    
    for(size_t k = start; k < end; ++k)
    {
        curr_pred = which_pred[k];
        
        // find nearest neighbors
        if(CROSS_VALIDATION)
        {
            temp_lib = which_lib;
            adjust_lib(curr_pred);
            nearest_neighbors = find_nearest_neighbors(distances[curr_pred]);
            which_lib = temp_lib;
        }
        else
        {
            nearest_neighbors = find_nearest_neighbors(distances[curr_pred]);
        }
        effective_nn = nearest_neighbors.size();
        
        if(effective_nn == 0)
        {
            predicted[curr_pred] = qnan;
            LOG_WARNING("no nearest neighbors found; using NA for forecast");
            continue;
        }
        
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
            predicted[curr_pred] += weights[k] * targets[nearest_neighbors[k]];
        predicted[curr_pred] = predicted[curr_pred] / total_weight;
        
        //compute variance
        predicted_var[curr_pred] = 0;
        for(size_t k = 0; k < effective_nn; ++k)
            predicted_var[curr_pred] += weights[k] * pow(targets[nearest_neighbors[k]] - predicted[curr_pred], 2);
        predicted_var[curr_pred] = predicted_var[curr_pred] / total_weight;
    }
    return;
}

void ForecastMachine::smap_prediction(const size_t start, const size_t end)
{
    size_t curr_pred, effective_nn, E = data_vectors[0].size();
    double avg_distance;
    //    vec weights;
    std::vector<size_t> nearest_neighbors;
    MatrixXd A, S_inv;
    VectorXd B, S, x, weights;
    double max_s, pred;
    std::vector<size_t> temp_lib;
    
    for(size_t k = start; k < end; ++k)
    {
        curr_pred = which_pred[k];
        
        // find nearest neighbors
        if(CROSS_VALIDATION)
        {
            temp_lib = which_lib;
            adjust_lib(curr_pred);
            nearest_neighbors = find_nearest_neighbors(distances[curr_pred]);
            which_lib = temp_lib;
        }
        else
        {
            nearest_neighbors = find_nearest_neighbors(distances[curr_pred]);
        }
        effective_nn = nearest_neighbors.size();
        
        if(effective_nn == 0)
        {
            predicted[curr_pred] = qnan;
            LOG_WARNING("no nearest neighbors found; using NA for forecast");
            continue;
        }
        weights = Eigen::VectorXd::Constant(effective_nn, 1.0);
        //        weights.assign(effective_nn, 1.0); // default is for theta = 0
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
                weights(i) = exp(-theta * distances[curr_pred][nearest_neighbors[i]] / avg_distance);
        }
        
        // setup matrices for SVD
        A.resize(effective_nn, E+1);
        B.resize(effective_nn);
        
        for(size_t i = 0; i < effective_nn; ++i)
        {
            B(i) = weights(i) * targets[nearest_neighbors[i]];
            
            for(size_t j = 0; j < E; ++j)
                A(i, j) = weights(i) * data_vectors[nearest_neighbors[i]][j];
            A(i, E) = weights(i);
        }
        
        // perform SVD
        Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        
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
        if(SAVE_SMAP_COEFFICIENTS)
        {
            for(size_t j = 0; j <= E; ++j)
                smap_coefficients[curr_pred][j] = x(j);
        }
        // save prediction
        predicted[curr_pred] = pred;
        
        // compute variance of prediction
        VectorXd w_resid = B - A * x;
        double total_w = 0;
        for(size_t i = 0; i < effective_nn; ++i)
        {
            total_w += weights(i) * weights(i);
        }
        predicted_var[curr_pred] = w_resid.dot(w_resid) / total_w;
    }
    return;
}

void ForecastMachine::const_prediction(const size_t start, const size_t end)
{
    size_t curr_pred;
    for(size_t k = start; k < end; ++k)
    {
        curr_pred = which_pred[k];
        const_predicted[curr_pred] = const_targets[curr_pred];
    }
    return;
}

void ForecastMachine::adjust_lib(const size_t curr_pred)
{
    // clear out lib indices we don't want from which_lib
    if(exclusion_radius >= 0)
    {
        auto f = [&](const size_t curr_lib) {
            return (curr_lib == curr_pred) || ((time[curr_lib] >= (time[curr_pred] - exclusion_radius)) && (time[curr_lib] <= (time[curr_pred] + exclusion_radius)));
        };
        which_lib.erase(std::remove_if(which_lib.begin(), which_lib.end(), f), which_lib.end());
    }
    else
    {
        which_lib.erase(std::remove(which_lib.begin(), which_lib.end(), curr_pred), which_lib.end());
    }
    return;
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

std::vector<size_t> sort_indices(const vec& v, std::vector<size_t> idx)
{
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

PredStats compute_stats_internal(const vec& obs, const vec& pred)
{
    size_t num_pred = 0;
    double sum_errors = 0;
    double sum_squared_errors = 0;
    double sum_obs = 0;
    double sum_pred = 0;
    double sum_squared_obs = 0;
    double sum_squared_pred = 0;
    double sum_prod = 0;
    size_t same_sign = 0;
    size_t num_vectors = obs.size();
    if(pred.size() < num_vectors)
        num_vectors = pred.size();
    
    for(size_t k = 0; k < num_vectors; ++k)
    {
        if(!std::isnan(obs[k]) && !std::isnan(pred[k]))
        {
            ++ num_pred;
            sum_errors += fabs(obs[k] - pred[k]);
            sum_squared_errors += (obs[k] - pred[k]) * (obs[k] - pred[k]);
            sum_obs += obs[k];
            sum_pred += pred[k];
            sum_squared_obs += obs[k] * obs[k];
            sum_squared_pred += pred[k] * pred[k];
            sum_prod += obs[k] * pred[k];
            if((obs[k] >= 0 && pred[k] >= 0) ||
               (obs[k] <= 0 && pred[k] <= 0))
                ++ same_sign;
        }
    }
    
    PredStats output;
    output.num_pred = num_pred;
    output.rho = (sum_prod * num_pred - sum_obs * sum_pred) /
        sqrt((sum_squared_obs * num_pred - sum_obs * sum_obs) *
            (sum_squared_pred * num_pred - sum_pred * sum_pred));
    output.mae = sum_errors / double(num_pred);
    output.rmse = sqrt(sum_squared_errors / double(num_pred));
    output.perc = double(same_sign) / double(num_pred);
    output.p_val = R::pnorm(atanh(output.rho), 0.0, 1.0 / sqrt(double(output.num_pred-3)), 0.0, 0.0);
    
    return output;
}

// [[Rcpp::export]]
DataFrame compute_stats(std::vector<double> observed, std::vector<double> predicted)
{
    PredStats output = compute_stats_internal(observed, predicted);
    return DataFrame::create( Named("num_pred") = output.num_pred,
                              Named("rho") = output.rho,
                              Named("mae") = output.mae,
                              Named("rmse") = output.rmse,
                              Named("perc") = output.perc,
                              Named("p_val") = output.p_val);
}

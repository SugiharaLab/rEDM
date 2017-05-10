#ifndef DATA_TYPES_H
#define DATA_TYPES_H

// shortcut name for vector<double> for in attractor reconstruction
typedef std::vector<double> vec;
typedef std::pair<size_t, size_t> time_range;

// which prediction method to use
enum PredEnum
{
    SIMPLEX,
    SMAP,
    FAST_LINEAR
};

enum NormEnum
{
    L1_NORM,
    L2_NORM, 
    P_NORM
};

struct PredStats
{
    size_t num_pred;
    double rho;
    double mae;
    double rmse;
    double perc;
    double p_val;
};

#endif

#ifndef DATA_TYPES_H
#define DATA_TYPES_H

using namespace std;

// shortcut name for vector<double> for in attractor reconstruction
typedef vector<double> vec;
typedef pair<size_t, size_t> time_range;

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
    L2_NORM
};

#endif
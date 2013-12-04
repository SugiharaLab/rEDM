#ifndef DATA_TYPES_H
#define DATA_TYPES_H

using namespace std;

// an ordered pair that specifies a range of vectors from 'first' to 'last'
struct time_range 
{
  int first, last;
};

// shortcut name for vector<double> for in attractor reconstruction
typedef vector<double> vec;

// which prediction method to use
enum PredEnum
{
  SIMPLEX,
  SMAP,
  FAST_LINEAR
};

#endif
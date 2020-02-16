#ifndef RCPPEDMCOMMON
#define RCPPEDMCOMMON

#define RCPPTHREAD_OVERRIDE_COUT 1 // std::cout override

#include <Rcpp.h>
#include <R.h>
#include <RcppThread.h>
#include <iostream>
#include "Common.h"

namespace r = Rcpp;

// Forward declarations
DataFrame< double > DFToDataFrame ( Rcpp::DataFrame df );

r::DataFrame DataFrameToDF ( DataFrame< double > dataFrame );

r::DataFrame ReadDataFrame ( std::string path, std::string file );

r::DataFrame PredictNonlinear_rcpp( std::string  pathIn,
                                    std::string  dataFile,
                                    r::DataFrame dataList,
                                    std::string  pathOut,
                                    std::string  predictFile,
                                    std::string  lib,
                                    std::string  pred,
                                    std::string  theta,
                                    int          E,
                                    int          Tp,
                                    int          knn,
                                    int          tau,
                                    std::string  columns,
                                    std::string  target,
                                    bool         embedded,
                                    bool         verbose,
                                    unsigned     numThreads );

r::DataFrame PredictInterval_rcpp( std::string  pathIn,
                                   std::string  dataFile,
                                   r::DataFrame dataList,
                                   std::string  pathOut,
                                   std::string  predictFile,
                                   std::string  lib,
                                   std::string  pred,
                                   int          maxTp,
                                   int          E,
                                   int          tau,
                                   std::string  columns,
                                   std::string  target,
                                   bool         embedded,
                                   bool         verbose,
                                   unsigned     numThreads );

r::List Multiview_rcpp ( std::string  pathIn,
                         std::string  dataFile,
                         r::DataFrame dataList,
                         std::string  pathOut,
                         std::string  predictFile,
                         std::string  lib,
                         std::string  pred,
                         int          E,
                         int          Tp,
                         int          knn,
                         int          tau, 
                         std::string  columns,
                         std::string  target,
                         int          multiview,
                         int          exlcusionRadius,
                         bool         verbose,
                         unsigned int numThreads );

r::DataFrame EmbedDimension_rcpp( std::string  pathIn,
                                  std::string  dataFile,
                                  r::DataFrame dataList,
                                  std::string  pathOut,
                                  std::string  predictFile,
                                  std::string  lib,
                                  std::string  pred,
                                  int          maxE,
                                  int          Tp,
                                  int          tau,
                                  std::string  columns,
                                  std::string  target,
                                  bool         embedded,
                                  bool         verbose,
                                  unsigned     numThreads );

r::DataFrame Embed_rcpp( std::string  path,
                         std::string  dataFile,
                         r::DataFrame df,
                         int          E,
                         int          tau,
                         std::string  columns,
                         bool         verbose );

r::DataFrame MakeBlock_rcpp( r::DataFrame             dataList,
                             int                      E,
                             int                      tau,
                             std::vector<std::string> columnNames,
                             bool                     verbose );

r::List ComputeError_rcpp ( std::vector<double> vec1, 
                            std::vector<double> vec2 );

r::DataFrame CCM_rcpp( std::string  pathIn,
                       std::string  dataFile,
                       r::DataFrame dataList,
                       std::string  pathOut,
                       std::string  predictFile,
                       int          E,
                       int          Tp,
                       int          knn,
                       int          tau, 
                       std::string  columns,
                       std::string  target,
                       std::string  libSizes,
                       int          sample,
                       bool         random,
                       bool         replacement,
                       unsigned     seed, 
                       bool         verbose );

r::DataFrame Simplex_rcpp( std::string  pathIn,
                           std::string  dataFile,
                           r::DataFrame dataList,
                           std::string  pathOut,
                           std::string  predictFile,
                           std::string  lib,
                           std::string  pred, 
                           int          E,
                           int          Tp,
                           int          knn,
                           int          tau, 
                           int          exclusionRadius, 
                           std::string  columns,
                           std::string  target,
                           bool         embedded,
                           bool         const_predict,
                           bool         verbose );

r::List SMap_rcpp( std::string  pathIn, 
                   std::string  dataFile,
                   r::DataFrame dataList,
                   std::string  pathOut,
                   std::string  predictFile,
                   std::string  lib,
                   std::string  pred, 
                   int          E,
                   int          Tp,
                   int          knn,
                   int          tau,
                   double       theta,
                   int          exclusionRadius, 
                   std::string  columns,
                   std::string  target,
                   std::string  smapFile,
                   std::string  jacobians,
                   bool         embedded,
                   bool         const_predict,
                   bool         verbose );
#endif

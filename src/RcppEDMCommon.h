
// R to C++ interface using Rcpp
// Functional flow: R func calls Rcpp func calls C++ func.

#ifndef RCPPEDMCOMMON
#define RCPPEDMCOMMON

#define RCPPTHREAD_OVERRIDE_COUT 1 // std::cout override

#include <Rcpp.h>
#include <R.h>
#include <RcppThread.h>
#include <iostream>
#include "API.h"

namespace r = Rcpp;

// Forward declarations
DataFrame< double > DFToDataFrame ( Rcpp::DataFrame df );

r::DataFrame DataFrameToDF ( DataFrame< double > dataFrame );

r::DataFrame ReadDataFrame ( std::string path, std::string file );

r::List ParamMaptoList( std::map< std::string, std::string > m );

r::List Simplex_rcpp( std::string       pathIn,
                      std::string       dataFile,
                      r::DataFrame      dataList,
                      std::string       pathOut,
                      std::string       predictFile,
                      std::string       lib,
                      std::string       pred,
                      int               E,
                      int               Tp,
                      int               knn,
                      int               tau,
                      int               exclusionRadius,
                      std::string       columns,
                      std::string       target,
                      bool              embedded,
                      // bool           const_predict, // Rcpp 20 arg limit
                      bool              verbose,
                      std::vector<bool> validLib,
                      int               generateSteps,
                      // bool           generateLibrary, // Rcpp 20 arg limit
                      bool              parameterList );

r::List SMap_rcpp( std::string       pathIn,
                   std::string       dataFile,
                   r::DataFrame      dataList,
                   //std::string     pathOut,     // Rcpp 20 arg limit
                   //std::string     predictFile, // Rcpp 20 arg limit
                   std::string       lib,
                   std::string       pred,
                   int               E,
                   int               Tp,
                   int               knn,
                   int               tau,
                   double            theta,
                   int               exclusionRadius,
                   std::string       columns,
                   std::string       target,
                   //std::string     smapCoefFile,  // Rcpp 20 arg limit
                   //std::string     smapSVFile,    // Rcpp 20 arg limit
                   //SVDValues       (*solver)      // Not supported by glmnet
                   bool              embedded,
                   //bool            const_predict, // Rcpp 20 arg limit
                   bool              verbose,
                   std::vector<bool> validLib,
                   bool              ignoreNan,
                   int               generateSteps,
                   //bool            generateLibrary, // Rcpp 20 arg limit
                   bool              parameterList );

r::List CCM_rcpp( std::string  pathIn,
                  std::string  dataFile,
                  r::DataFrame dataList,
                  //std::string  pathOut,     // Rcpp 20 arg limit
                  //std::string  predictFile, // Rcpp 20 arg limit
                  int          E,
                  int          Tp,
                  int          knn,
                  int          tau,
                  int          exclusionRadius,
                  std::string  columns,
                  std::string  target,
                  std::string  libSizes,
                  int          sample,
                  bool         random,
                  // bool      replacement,  // Rcpp 20 arg limit
                  unsigned     seed,
                  bool         embedded,
                  bool         includeData,
                  bool         parameterList,
                  bool         verbose );

r::List Multiview_rcpp ( std::string  pathIn,
                         std::string  dataFile,
                         r::DataFrame dataList,
                         //std::string  pathOut,      // Rcpp 20 arg limit
                         //std::string  predictFile,  // Rcpp 20 arg limit
                         std::string  lib,
                         std::string  pred,
                         int          D,
                         int          E,
                         int          Tp,
                         int          knn,
                         int          tau, 
                         std::string  columns,
                         std::string  target,
                         int          multiview,
                         int          exlcusionRadius,
                         bool         trainLib,
                         bool         excludeTarget,
                         bool         parameterList,
                         bool         verbose,
                         unsigned int numThreads );

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
                                    int          exclusionRadius,
                                    std::string  columns,
                                    std::string  target,
                                    bool         embedded,
                                    bool         verbose,
                                    std::vector<bool> validLib,
                                    bool         ignoreNan,
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
                                   int          exclusionRadius,
                                   std::string  columns,
                                   std::string  target,
                                   bool         embedded,
                                   bool         verbose,
                                   std::vector<bool> validLib,
                                   unsigned     numThreads );

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
                                  int          exclusionRadius,
                                  std::string  columns,
                                  std::string  target,
                                  bool         embedded,
                                  bool         verbose,
                                  std::vector<bool> validLib,
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
                             bool                     deletePartial );

r::List ComputeError_rcpp ( std::vector<double> vec1, 
                            std::vector<double> vec2 );
#endif

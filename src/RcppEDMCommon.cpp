// Expose and map cpp wrapper functions to EDM module via Rcpp
// See RCPP_MODULE() at end of file.
//
// Requirements for a function to be exposed to R via Rcpp modules are:
//   The function takes between 0 and 65 parameters.
//   The type of each input parameter must be manageable by Rcpp::as template.
//   The return type of the function must be either void or any type that can
//   be managed by the Rcpp::wrap template.
//   The function name itself has to be unique in the module. In other words,
//   no two functions with the same name but different signatures are allowed.
//   C++ allows overloading functions. This might be added in future versions
//   of modules.

#include "RcppEDMCommon.h"

//-------------------------------------------------------------------------
// Definitions of formal arguments and default params of the R functions
// that encapsulate the C++ functions in an Rcpp::List.
//-------------------------------------------------------------------------
auto ReadDataFrameArgs = r::List::create( r::_["path"] = "",
                                          r::_["file"] = "" );

auto MakeBlockArgs = r::List::create( 
    r::_["dataFrame"]     = r::DataFrame(),
    r::_["E"]             = 0,
    r::_["tau"]           = -1,
    r::_["columnNames"]   = std::vector<std::string>(),
    r::_["deletePartial"] = false );

auto EmbedArgs = r::List::create( 
    r::_["path"]      = std::string(""),
    r::_["dataFile"]  = std::string(""),
    r::_["dataFrame"] = r::DataFrame(),
    r::_["E"]         = 0,
    r::_["tau"]       = -1,
    r::_["columns"]   = std::string(""),
    r::_["verbose"]   = false );

auto SimplexArgs = r::List::create(
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    r::_["pathOut"]         = std::string("./"),
    r::_["predictFile"]     = std::string(""),
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["E"]               = 0,
    r::_["Tp"]              = 1,
    r::_["knn"]             = 0,
    r::_["tau"]             = -1,
    r::_["exclusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["embedded"]        = false,
    r::_["const_predict"]   = false,
    r::_["verbose"]         = false,
    r::_["validLib"]        = std::vector<bool>(),
    r::_["generateSteps"]   = 0,
    r::_["parameterList"]   = false );

auto SMapArgs = r::List::create(
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    //r::_["pathOut"]       = std::string("./"), // Rcpp 20 arg limit
    //r::_["predictFile"]   = std::string(""),   // Rcpp 20 arg limit
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["E"]               = 0,
    r::_["Tp"]              = 1,
    r::_["knn"]             = 0,
    r::_["tau"]             = -1,
    r::_["theta"]           = 0,
    r::_["exclusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["smapFile"]        = std::string(""),
    //r::_["jacobians"]     = std::string(""), // Rcpp 20 arg limit
    r::_["embedded"]        = false,
    r::_["const_predict"]   = false,
    r::_["verbose"]         = false,
    r::_["validLib"]        = std::vector<bool>(),
    r::_["generateSteps"]   = 0,
    r::_["parameterList"]   = false );

auto MultiviewArgs = r::List::create( 
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    //r::_["pathOut"]       = std::string("./"), // Rcpp 20 arg limit
    //r::_["predictFile"]   = std::string(""),   // Rcpp 20 arg limit
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["D"]               = 0,
    r::_["E"]               = 1,
    r::_["Tp"]              = 1,
    r::_["knn"]             = 0,
    r::_["tau"]             = -1,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["multiview"]       = 0,
    r::_["exlcusionRadius"] = 0,
    r::_["trainLib"]        = true,
    r::_["excludeTarget"]   = false,
    r::_["parameterList"]   = false,
    r::_["verbose"]         = false,
    r::_["numThreads"]      = 4 );

auto CCMArgs = r::List::create( 
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    //r::_["pathOut"]       = std::string("./"), // Rcpp 20 arg limit
    //r::_["predictFile"]   = std::string(""),   // Rcpp 20 arg limit
    r::_["E"]               = 0,
    r::_["Tp"]              = 0,
    r::_["knn"]             = 0,
    r::_["tau"]             = -1,
    r::_["exlcusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["libSizes"]        = std::string(""),
    r::_["sample"]          = 0,
    r::_["random"]          = true,
    r::_["replacement"]     = false,
    r::_["seed"]            = 0,
    r::_["embedded"]        = false,
    r::_["includeData"]     = false,
    r::_["parameterList"]   = false,
    r::_["verbose"]         = false );
    
auto EmbedDimensionArgs     = r::List::create( 
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    r::_["pathOut"]         = std::string("./"),
    r::_["predictFile"]     = std::string(""),
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["maxE"]            = 10,
    r::_["Tp"]              = 1,
    r::_["tau"]             = -1,
    r::_["exclusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["embedded"]        = false,
    r::_["verbose"]         = false,
    r::_["validLib"]        = std::vector<bool>(),
    r::_["numThreads"]      = 4 );

auto PredictIntervalArgs = r::List::create( 
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    r::_["pathOut"]         = std::string("./"),
    r::_["predictFile"]     = std::string(""),
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["maxTp"]           = 10,
    r::_["E"]               = 0,
    r::_["tau"]             = -1,
    r::_["exclusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["embedded"]        = false,
    r::_["verbose"]         = false,
    r::_["validLib"]        = std::vector<bool>(),
    r::_["numThreads"]      = 4 );

auto PredictNonlinearArgs = r::List::create( 
    r::_["pathIn"]          = std::string("./"),
    r::_["dataFile"]        = std::string(""),
    r::_["dataFrame"]       = r::DataFrame(),
    r::_["pathOut"]         = std::string("./"),
    r::_["predictFile"]     = std::string(""),
    r::_["lib"]             = std::string(""),
    r::_["pred"]            = std::string(""),
    r::_["theta"]           = std::string(""),
    r::_["E"]               = 0,
    r::_["Tp"]              = 1,
    r::_["knn"]             = 0,
    r::_["tau"]             = -1,
    r::_["exclusionRadius"] = 0,
    r::_["columns"]         = std::string(""),
    r::_["target"]          = std::string(""),
    r::_["embedded"]        = false,
    r::_["verbose"]         = false,
    r::_["validLib"]        = std::vector<bool>(),
    r::_["numThreads"]      = 4 );

//-------------------------------------------------------------------------
// Export / map the functions
//   First argument:  R function name, see ../R/EDM.R
//   Second argument: pointer to Rcpp interface function
//   Third argument:  arguments of the R function that encapsulates the 
//                    C++ function in a Rcpp::List
//-------------------------------------------------------------------------
RCPP_MODULE(EDMInternal) {
    r::function( "RtoCpp_ComputeError",  &ComputeError_rcpp                  );
    r::function( "RtoCpp_ReadDataFrame", &ReadDataFrame,   ReadDataFrameArgs );
    r::function( "RtoCpp_MakeBlock",     &MakeBlock_rcpp,  MakeBlockArgs     );
    r::function( "RtoCpp_Embed",         &Embed_rcpp,      EmbedArgs         );
    r::function( "RtoCpp_Simplex",       &Simplex_rcpp,    SimplexArgs       );
    r::function( "RtoCpp_SMap",          &SMap_rcpp,       SMapArgs          );
    r::function( "RtoCpp_Multiview",     &Multiview_rcpp,  MultiviewArgs     );
    r::function( "RtoCpp_CCM",           &CCM_rcpp,        CCMArgs           );
    r::function( "RtoCpp_EmbedDimension",   &EmbedDimension_rcpp, 
                                             EmbedDimensionArgs   );
    r::function( "RtoCpp_PredictInterval",  &PredictInterval_rcpp, 
                                             PredictIntervalArgs  );
    r::function( "RtoCpp_PredictNonlinear", &PredictNonlinear_rcpp, 
                                             PredictNonlinearArgs );
}

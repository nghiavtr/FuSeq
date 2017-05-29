/*
Date:01/11/2016
Note:This implementation is adapted from CollapsedEMOptimizer.hpp of Sailfish for our purpose. 
*/

#ifndef EXPORT_FEQ_HPP
#define EXPORT_FEQ_HPP

#include <memory>
#include <unordered_map>
#include <functional>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "ReadExperiment.hpp"
#include "SailfishOpts.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

//for FuSeq
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
//for FuSeq


//class BootstrapWriter;

class ExportFeq {
    public:
        ExportFeq();

        bool writeFeq(ReadExperiment& readExp,
                      SailfishOpts& sopt,
                      const boost::filesystem::path& feqOutputPath = "/" //for FuSeq                      
                      );
             
};

#endif // EXPORT_FEQ_HPP


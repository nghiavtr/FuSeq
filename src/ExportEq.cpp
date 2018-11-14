/*
Date:08 Nov 2018
Note:This implementation is adapted from CollapsedEMOptimizer.cpp of Sailfish for our purpose. 
*/
#include <vector>
#include <unordered_map>
#include <atomic>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#include "concurrentqueue.h"

#include <boost/math/special_functions/digamma.hpp>
//for FuSeq
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "EmpiricalDistribution.hpp"
#include "ReadLibrary.hpp"
#include "RapMapUtils.hpp"
//for FuSeq

// C++ string formatting library
#include "spdlog/details/format.h"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

#include "ExportEq.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishMath.hpp"
#include "ReadExperiment.hpp"
#include "BootstrapWriter.hpp"
#include "MultinomialSampler.hpp"

ExportEq::ExportEq() {}


bool ExportEq::writeEq(ReadExperiment& readExp,
        SailfishOpts& sopt,
        const boost::filesystem::path& feqOutputPath //for FuSeq
        ){ 
    
    std::vector<Transcript>& transcripts = readExp.transcripts();
    Eigen::VectorXd effLens(transcripts.size());

    // Fill in the effective length vector
    double totalLen{0.0};
    for (size_t i = 0; i < transcripts.size(); ++i) {
        effLens(i) = (sopt.noEffectiveLengthCorrection) ?
                        transcripts[i].RefLength : transcripts[i].EffectiveLength;
        if (effLens(i) <= 1.0) { effLens(i) = 1.0; }
        totalLen += effLens(i);
    }

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    
    std::cout << "=========================================================================\n" << std::endl;
    std::cout << "[FuSeq] -- We are export some data here -- \n" << std::endl;
    std::cout << "=========================================================================\n" << std::endl;
        std::cout << "[FuSeq] -- Export fragment information \n" << std::endl; 
    boost::filesystem::path fragfname = feqOutputPath / "fragmentInfo.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> FragOutput(std::fopen(fragfname.c_str(), "w"), std::fclose);

    auto kmerLen = rapmap::utils::my_mer::k();
    auto& rl = readExp.readLibraries().front();
    if (rl.format().type == ReadType::SINGLE_END){
        fmt::print(FragOutput.get(), "readlen\tnumObservedReads\tnumMappedReads\tnumHits\tkmer\n");
        fmt::print(FragOutput.get(), "{}\t{}\t{}\t{}\t{}\n",readExp.getReadLength(),readExp.numObservedFragments(),readExp.numMappedFragments(),readExp.numFragHitsAtomic(),kmerLen);
        std::cout <<"Read length "<<readExp.getReadLength()<< "Observed reads "<<readExp.numObservedFragments() << " Mapped reads "<<readExp.numMappedFragments() <<" Total hits " <<readExp.numFragHitsAtomic()<< " kmer "<<kmerLen<< "\n" << std::endl;
    }else{
        fmt::print(FragOutput.get(), "readlen\tfragLengthMedian\tfragLengthMean\tfragLengthSd\tnumObservedFragments\tnumMappedFragments\tnumHits\tkmer\n");
        fmt::print(FragOutput.get(), "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",readExp.getReadLength(),readExp.fragLengthDist()->median(),readExp.fragLengthDist()->mean(),readExp.fragLengthDist()->sd(),readExp.numObservedFragments(),readExp.numMappedFragments(),readExp.numFragHitsAtomic(),kmerLen);
        std::cout <<"Read length "<<readExp.getReadLength()<< " fragLengthMedian " << readExp.fragLengthDist()->median() << " fragLengthMean " << readExp.fragLengthDist()->mean() << " fragLengthSd " <<  readExp.fragLengthDist()->sd() << "Observed Fragments "<<readExp.numObservedFragments() << " Mapped Fragments "<<readExp.numMappedFragments() <<" Total hits " <<readExp.numFragHitsAtomic()<< " kmer "<<kmerLen<< "\n" << std::endl;
    }
    std::cout << "[FuSeq] -- Extracting equivalence class \n" << std::endl;
    boost::filesystem::path sequgiofname = feqOutputPath / "eqClass.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> SequgioOutput(std::fopen(sequgiofname.c_str(), "w"), std::fclose);
    fmt::print(SequgioOutput.get(), "Transcript\tWeight\tCount\tEffLength\tRefLength\teqClass\n");
    // we can export information of transcript mapped read counts before bias correction and optimization
    uint32_t eqClassID=0;
    for (auto& kv : eqVec) {
        eqClassID=eqClassID+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;          
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            fmt::print(SequgioOutput.get(), "{}\t{}\t{}\t{}\t{}\t{}\n",transcripts[t].RefName,v.weights[i],v.count,effLens(t),transcripts[t].RefLength ,eqClassID);
    
        }
    }
    return true;
}


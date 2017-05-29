/*
Date:01/11/2016
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

#include "ExportFeq.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishMath.hpp"
#include "ReadExperiment.hpp"
#include "BootstrapWriter.hpp"
#include "MultinomialSampler.hpp"

ExportFeq::ExportFeq() {}


bool ExportFeq::writeFeq(ReadExperiment& readExp,
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
    fmt::print(FragOutput.get(), "readlen\tfragLengthMedian\tfragLengthMean\tfragLengthSd\tnumObservedFragments\tnumMappedFragments\tnumHits\tkmer\n");
    fmt::print(FragOutput.get(), "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",readExp.getReadLength(),readExp.fragLengthDist()->median(),readExp.fragLengthDist()->mean(),readExp.fragLengthDist()->sd(),readExp.numObservedFragments(),readExp.numMappedFragments(),readExp.numFragHitsAtomic(),kmerLen);
    std::cout <<"Read length "<<readExp.getReadLength()<< " fragLengthMedian " << readExp.fragLengthDist()->median() << " fragLengthMean " << readExp.fragLengthDist()->mean() << " fragLengthSd " <<  readExp.fragLengthDist()->sd() << "Observed Fragments "<<readExp.numObservedFragments() << " Mapped Fragments "<<readExp.numMappedFragments() <<" Total hits " <<readExp.numFragHitsAtomic()<< " kmer "<<kmerLen<< "\n" << std::endl;
       
    std::cout << "[FuSeq] -- Extracting equivalence class \n" << std::endl;
    boost::filesystem::path sequgiofname = feqOutputPath / "rawCount.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> SequgioOutput(std::fopen(sequgiofname.c_str(), "w"), std::fclose);
    fmt::print(SequgioOutput.get(), "Transcript\tWeight\tCount\teffLens\tRefLength\tEffectiveLength\teqClass\n");
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
            fmt::print(SequgioOutput.get(), "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",transcripts[t].RefName,v.weights[i],v.count ,effLens(t),transcripts[t].RefLength ,transcripts[t].EffectiveLength, eqClassID);
    
        }
    }

    std::cout << "[FuSeq] -- Extracting RR fusion equivalence classes \n" << std::endl;
    std::vector<std::pair<const TranscriptGroup, TGValue>>& feqVecRR = readExp.feqClassBuilderRR().eqVec();

    boost::filesystem::path feqfnameRR = feqOutputPath / "feq_RR.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> feqOutputRR(std::fopen(feqfnameRR.c_str(), "w"), std::fclose);
    fmt::print(feqOutputRR.get(),"Transcript\tCount\tRead\tFeq\n");

    uint32_t feqClassID_RR=0;
    for (auto& kv : feqVecRR) {
        feqClassID_RR=feqClassID_RR+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;
        uint32_t readType=1;
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            if (t==0){                
                readType=2;
            }else{
                fmt::print(feqOutputRR.get(), "{}\t{}\t{}\t{}\n",transcripts[t-1].RefName,v.count, readType, feqClassID_RR);
            }
         }
    }

    std::cout << "[FuSeq] -- Extracting RF fusion equivalence classes \n" << std::endl;
    std::vector<std::pair<const TranscriptGroup, TGValue>>& feqVecRF = readExp.feqClassBuilderRF().eqVec();

    boost::filesystem::path feqfnameRF = feqOutputPath / "feq_RF.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> feqOutputRF(std::fopen(feqfnameRF.c_str(), "w"), std::fclose);
    fmt::print(feqOutputRF.get(),"Transcript\tCount\tRead\tFeq\n");

    uint32_t feqClassID_RF=0;
    for (auto& kv : feqVecRF) {
        feqClassID_RF=feqClassID_RF+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;
        uint32_t readType=1;
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            if (t==0){                
                readType=2;
            }else{
                fmt::print(feqOutputRF.get(), "{}\t{}\t{}\t{}\n",transcripts[t-1].RefName,v.count, readType, feqClassID_RF);
            }
         }
    }

        std::cout << "[FuSeq] -- Extracting FF fusion equivalence classes \n" << std::endl;
    std::vector<std::pair<const TranscriptGroup, TGValue>>& feqVecFF = readExp.feqClassBuilderFF().eqVec();

    boost::filesystem::path feqfnameFF = feqOutputPath / "feq_FF.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> feqOutputFF(std::fopen(feqfnameFF.c_str(), "w"), std::fclose);
    fmt::print(feqOutputFF.get(),"Transcript\tCount\tRead\tFeq\n");

    uint32_t feqClassID_FF=0;
    for (auto& kv : feqVecFF) {
        feqClassID_FF=feqClassID_FF+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;
        uint32_t readType=1;
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            if (t==0){                
                readType=2;
            }else{
                fmt::print(feqOutputFF.get(), "{}\t{}\t{}\t{}\n",transcripts[t-1].RefName,v.count, readType, feqClassID_FF);
            }
         }
    }

    std::cout << "[FuSeq] -- Extracting FR fusion equivalence classes \n" << std::endl;
    std::vector<std::pair<const TranscriptGroup, TGValue>>& feqVecFR = readExp.feqClassBuilderFR().eqVec();

    boost::filesystem::path feqfnameFR = feqOutputPath / "feq_FR.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> feqOutputFR(std::fopen(feqfnameFR.c_str(), "w"), std::fclose);
    fmt::print(feqOutputFR.get(),"Transcript\tCount\tRead\tFeq\n");

    uint32_t feqClassID_FR=0;
    for (auto& kv : feqVecFR) {
        feqClassID_FR=feqClassID_FR+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;
        uint32_t readType=1;
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            if (t==0){                
                readType=2;
            }else{
                fmt::print(feqOutputFR.get(), "{}\t{}\t{}\t{}\n",transcripts[t-1].RefName,v.count, readType, feqClassID_FR);
            }
         }
    }


    std::cout << "[FuSeq] -- Extracting UN fusion equivalence classes \n" << std::endl;
    std::vector<std::pair<const TranscriptGroup, TGValue>>& feqVecUN = readExp.feqClassBuilderUN().eqVec();

    boost::filesystem::path feqfnameUN = feqOutputPath / "feq_UN.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE *)> feqOutputUN(std::fopen(feqfnameUN.c_str(), "w"), std::fclose);
    fmt::print(feqOutputUN.get(),"Transcript\tCount\tRead\tFeq\n");

    uint32_t feqClassID_UN=0;
    for (auto& kv : feqVecUN) {
        feqClassID_UN=feqClassID_UN+1;
        auto& tg = kv.first;
        // The size of the label
        size_t classSize = tg.txps.size();
        // The weights of the label
        TGValue& v = kv.second;
        uint32_t readType=1;
        for (size_t i = 0; i < classSize; ++i) {
            auto& t = tg.txps[i];
            if (t==0){                
                readType=2;
            }else{
                fmt::print(feqOutputUN.get(), "{}\t{}\t{}\t{}\n",transcripts[t-1].RefName,v.count, readType, feqClassID_UN);
            }
         }
    }

    return true;
}

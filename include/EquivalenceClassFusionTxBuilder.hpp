/*
Date:01/11/2016
Note:This implementation is adapted from EquivalenceClassBuilder.hpp of Sailfish for our purpose. 
*/

#ifndef EQUIVALENCE_CLASS_FUSION_TX_BUILDER_HPP
#define EQUIVALENCE_CLASS_FUSION_TX_BUILDER_HPP

#include <unordered_map>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>

// Logger includes
#include "spdlog/spdlog.h"

#include "cuckoohash_map.hh"
#include "concurrentqueue.h"
#include "TranscriptGroup.hpp"


struct TGFusionTxValue {
    TGFusionTxValue(const TGFusionTxValue& o) {
        weights = o.weights;
        count.store(o.count.load());

        minLeftReadPos.store(o.minLeftReadPos.load());
        maxLeftReadPos.store(o.maxLeftReadPos.load());

        minRightReadPos.store(o.minRightReadPos.load());
        maxRightReadPos.store(o.maxRightReadPos.load());

    }

    TGFusionTxValue(std::vector<double>& weightIn, uint64_t countIn, int32_t minLeftReadPosIn, int32_t maxLeftReadPosIn, int32_t minRightReadPosIn, int32_t maxRightReadPosIn) :
        weights(weightIn) { count.store(countIn); minLeftReadPos.store(minLeftReadPosIn);maxLeftReadPos.store(maxLeftReadPosIn); minRightReadPos.store(minRightReadPosIn);maxRightReadPos.store(maxRightReadPosIn);}

    // const is a lie
    void normalizeAux() const {
        double sumOfAux{0.0};
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux += weights[i];
        }
        double norm = 1.0 / sumOfAux;
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= norm;
        }
        /* LOG SPACE
        double sumOfAux = salmon::math::LOG_0;
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux = salmon::math::logAdd(sumOfAux, weights[i]);
        }
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] = std::exp(weights[i] - sumOfAux);
        }
        */
    }

    // forget synchronizing this for the time being
    mutable std::vector<double> weights;
    std::atomic<uint64_t> count{0};
    std::atomic<int32_t> minLeftReadPos{2147483646};
    std::atomic<int32_t> maxLeftReadPos{0};

    std::atomic<int32_t> minRightReadPos{2147483646};
    std::atomic<int32_t> maxRightReadPos{0};

};

class EquivalenceClassFusionTxBuilder {
    public:
        EquivalenceClassFusionTxBuilder(std::shared_ptr<spdlog::logger> loggerIn) :
		logger_(loggerIn) {
            countMap_.reserve(1000000);
        }

        ~EquivalenceClassFusionTxBuilder() {}

        void start() { active_ = true; }

        bool finish() {
            active_ = false;
            size_t totalCount{0};
            auto lt = countMap_.lock_table();
            for (auto& kv : lt) {
            //for (auto kv = countMap_.begin(); !kv.is_end(); ++kv) {
                kv.second.normalizeAux();
                totalCount += kv.second.count;
                countVec_.push_back(kv);
            }

    	    logger_->info("Computed {} rich equivalence classes "
			  "for further processing", countVec_.size());
            logger_->info("Counted {} total reads in the equivalence classes ",
                    totalCount);
            return true;
        }

        inline void insertGroup(TranscriptGroup g, uint32_t count, int32_t minLeftReadPos, int32_t maxLeftReadPos, int32_t minRightReadPos, int32_t maxRightReadPos) {
            std::vector<double> weights(g.txps.size(), 0.0);
            //auto updatefn = [count](TGFusionTxValue& x) { x.count = count; };
            TGFusionTxValue v(weights, count, minLeftReadPos, maxLeftReadPos, minRightReadPos, maxRightReadPos);
            countVec_.push_back(std::make_pair(g, v));
            //countMap_.upsert(g, updatefn, v);
        }

        inline void addGroup(TranscriptGroup&& g,
                             std::vector<double>& weights, int32_t leftReadPos, int32_t rightReadPos) {

            auto upfn = [&weights, &leftReadPos, &rightReadPos](TGFusionTxValue& x) -> void {
                // update the count
                x.count++;
                // update the weights
                for (size_t i = 0; i < x.weights.size(); ++i) {
                    // Possibly atomicized in the future
                    weights[i] += x.weights[i];
                    /* LOG SPACE
                    x.weights[i] =
                        salmon::math::logAdd(x.weights[i], weights[i]);
                    */
                }
                //update minPos and maxPos of the left
                int32_t old_minLeft = x.minLeftReadPos.load();
                while(old_minLeft > leftReadPos &&
                    x.minLeftReadPos.compare_exchange_weak(old_minLeft, leftReadPos));

                int32_t old_maxLeft = x.maxLeftReadPos.load();                
                while(old_maxLeft < leftReadPos &&
                    !x.maxLeftReadPos.compare_exchange_weak(old_maxLeft, leftReadPos));


                //update minPos and maxPos of the right
                int32_t old_minRight = x.minRightReadPos.load();
                while(old_minRight > rightReadPos &&
                    x.minRightReadPos.compare_exchange_weak(old_minRight, rightReadPos));

                int32_t old_maxRight = x.maxRightReadPos.load();                
                while(old_maxRight < rightReadPos &&
                    !x.maxRightReadPos.compare_exchange_weak(old_maxRight, rightReadPos));    

           
            };
            TGFusionTxValue v(weights, 1, leftReadPos, leftReadPos+1, rightReadPos, rightReadPos+1);
            countMap_.upsert(g, upfn, v);
        }

        std::vector<std::pair<const TranscriptGroup, TGFusionTxValue>>& eqVec() {
            return countVec_;
        }

    private:
        std::atomic<bool> active_;
	    cuckoohash_map<TranscriptGroup, TGFusionTxValue, TranscriptGroupHasher> countMap_;
        std::vector<std::pair<const TranscriptGroup, TGFusionTxValue>> countVec_;
    	std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_FUSION_BUILDER_HPP

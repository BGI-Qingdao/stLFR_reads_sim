#ifndef __APPCOMMON_RANDOM_H__
#define __APPCOMMON_RANDOM_H__

#include <random>
#include <map>

namespace BGIQD {
    namespace Random {

        int RandomStartPosByLength(int len )
        {
            std::default_random_engine generator;
            std::uniform_int_distribution<int> distribution(0,len-1);
            return distribution(generator);
        }

        struct DiscreteRandomWithBin
        {
            struct ConfIni
            {
                int start ;
                int bin ;
                int weight ;
            };

            std::vector<ConfIni> keybin;
            std::default_random_engine generator;
            std::discrete_distribution<int> distribution;

            void InitDistribution()
            {
                std::vector<int> weight ;
                for( const auto & x : keybin )
                    weight.push_back(x.weight);
                distribution = std::discrete_distribution<int>
                    (weight.begin() ,weight.end()) ;
            }
            int operator() ()
            {
                int index = distribution(generator);
                int step = RandomStartPosByLength(keybin.at(index).bin);
                return keybin.at(index).start + step ;
            }
        };
    }
}

#endif

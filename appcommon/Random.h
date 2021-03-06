#ifndef __APPCOMMON_RANDOM_H__
#define __APPCOMMON_RANDOM_H__

#include <random>
#include <map>
#include <cassert>
#include <chrono>

namespace BGIQD {
    namespace Random {

        //static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
        static std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());

        char RandomInsert() 
        {
            std::uniform_int_distribution<int> distribution(0,3);
            int to_index = distribution(generator);
            const static char table[] = "AGCT";
            return table[to_index];
        }

        char RandomSubstitute(char x)
        {
            std::uniform_int_distribution<int> distribution(0,2);
            int to_index = distribution(generator);
            const static char tableA[] = "GCT";
            const static char tableG[] = "ACT";
            const static char tableC[] = "GAT";
            const static char tableT[] = "GCA";
            switch(x)
            {
                case 'A':
                case 'a':
                    return tableA[to_index];
                case 'G':
                case 'g':
                    return tableG[to_index];
                case 'C':
                case 'c':
                    return tableC[to_index];
                case 'T':
                case 't':
                    return tableT[to_index];
            }
            assert(0);
            return 'N';
        }

        bool RandomByProbability(float prob )
        {
            int accuracy = 1000000 ;
            int hit = prob * accuracy ;
            return generator() % accuracy < hit ;
            //std::uniform_int_distribution<int> distribution(0,accuracy-1);
            //return distribution(generator) < hit ;

        }

        int RandomChoose( float p1 ,float p2 , float p3 )
        {
            int accuracy = 1000000 ;
            std::vector<int> hits;
            hits.push_back( p1 * accuracy) ;
            hits.push_back( p2 * accuracy) ;
            hits.push_back( p3 * accuracy) ;
            return  std::discrete_distribution<int>
                (hits.begin() ,hits.end())(generator) ;
        }

        long long RandomStartPosByLength(long  long len )
        {
            return generator() % len ;
            //std::uniform_int_distribution<long long > distribution(0,len-1);
            //return distribution(generator);
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

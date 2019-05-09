#ifndef __APPCOMMON_MUTTAION_H__
#define __APPCOMMON_MUTTAION_H__

#include "biocommon/align_common/align_result.h"
#include "biocommon/seq/seq.h"
#include "appcommon/Random.h"
#include <ctype.h>
#include <cassert>

namespace BGIQD {
    namespace Random{

        struct MutationResult
        {
            std::string seq;
            std::vector<BGIQD::ALIGN_COMMON::MatchInfo> details ;
        };

        struct MutationEngine
        {
            float mutation_rate ;
            float insert_percent ;
            float delete_percent ;
            float substitute_percent ;

            bool Valid() const {
                return mutation_rate >= 0.0f
                    && insert_percent >= 0.0f 
                    && delete_percent >= 0.0f 
                    && substitute_percent >= 0.0f 
                    && ( insert_percent + delete_percent + substitute_percent ) == 1.0f ;
            }

            MutationResult operator()
                (const std::string & s) const
                {
                    MutationResult ret ;
                    auto & ret_seq = ret.seq ;
                    BGIQD::ALIGN_COMMON::CIGAR last_cigar = BGIQD::ALIGN_COMMON::CIGAR::M;
                    int last_cigar_num = 0 ;

                    auto add_cigar = [&]( const BGIQD::ALIGN_COMMON::CIGAR & cig)
                    {
                        if( cig == last_cigar )
                            last_cigar_num ++ ;
                        else
                        {
                            if( last_cigar_num > 0 )
                            {
                                BGIQD::ALIGN_COMMON::MatchInfo tmp;
                                tmp.type = last_cigar ;
                                tmp.len = last_cigar_num ;
                                ret.details.push_back(tmp);
                            }
                            last_cigar = cig ;
                            last_cigar_num = 1;
                        }
                    };

                    for( const auto c : s )
                    {
                        assert( c == 'A' || c=='a'
                                || c== 'G' || c=='g'
                                || c== 'C' || c=='c'
                                || c== 'T' || c=='t' );
                        if( ! RandomByProbability(mutation_rate) )
                        {
                            ret_seq.push_back(toupper(c));
                            add_cigar(BGIQD::ALIGN_COMMON::CIGAR::EQUAL);
                        }
                        else
                        {
                            int type = RandomChoose( substitute_percent
                                    ,insert_percent ,delete_percent );
                            if( type == 0 )
                            {//substitute
                                ret_seq.push_back(
                                        RandomSubstitute(c));
                                add_cigar(BGIQD::ALIGN_COMMON::CIGAR::X);
                            }
                            else if ( type == 1 )
                            {// insert
                                ret_seq.push_back(
                                        RandomInsert());
                                add_cigar(BGIQD::ALIGN_COMMON::CIGAR::I);
                                ret_seq.push_back(toupper(c));
                                add_cigar(BGIQD::ALIGN_COMMON::CIGAR::EQUAL);
                            }
                            else
                            {//delete
                                add_cigar(BGIQD::ALIGN_COMMON::CIGAR::D);
                            }
                        }
                    }
                    // deal last cigar
                    if( last_cigar_num > 0 )
                    {
                        BGIQD::ALIGN_COMMON::MatchInfo tmp;
                        tmp.type = last_cigar ;
                        tmp.len = last_cigar_num ;
                        ret.details.push_back(tmp);
                    }
                    return ret ;
                }
        };
    }
}

#endif 

#ifndef __APPCOMMON_REF_H__
#define __APPCOMMON_REF_H__

#include "biocommon/fasta/fasta.h"

#include <vector>
#include <map>
#include <tuple>
#include <algorithm>


namespace BGIQD {
    namespace stLFRSim {

        typedef BGIQD::FASTA::NormalHead RefHead;

        typedef BGIQD::FASTA::Fasta<RefHead> RefFa;

        struct RefChromesome
        {
            int length ;
            RefFa fa ;
            //      [n_start , n_end ] by 0base 
            std::vector< std::tuple<int , int > > n_area ;

            void Init()
            {
                length = fa.seq.atcgs.length();
                int n_before = -1 ;
                bool is_n = false ;
                int i = 0 ;
                for( const auto & c : fa.seq.atcgs )
                {
                    switch(c)
                    {
                        case 'A':
                        case 'a':
                        case 'G':
                        case 'g':
                        case 'T':
                        case 't':
                        case 'C':
                        case 'c':
                            is_n = false ;
                            break;
                        default :
                            is_n = true ;
                    }
                    if( n_before < 0 && is_n )
                        n_before = i ;
                    else if ( !is_n && n_before > 0 )
                    {
                        n_area.push_back( 
                                std::make_tuple( n_before, i)) ;
                        n_before = -1 ;
                    }
                    i ++ ;
                }
                if ( n_before > 0 )
                {
                    n_area.push_back( 
                            std::make_tuple( n_before, i)) ;
                }
            }
            bool IsValidArea( int start , int length )
            {
                auto itr = std::lower_bound( n_area.begin() 
                        , n_area.end() 
                        , std::make_tuple( start +length -1, 100000000L) );
                if( itr != n_area.end() )
                {
                    int tmp_start ; int tmp_end ;
                    std::tie(tmp_start,tmp_end) = *itr ;
                    if( start <= tmp_end )
                        return false ;
                }
                return true ;
            }
        };

        struct Ref
        {
            std::vector<RefChromesome> refs;
            //      [n_start , n_end ] by 0base 
            std::vector<std::tuple<long , long>>  chromesome_area ;
            void InitAreas()
            {
                long curr_pos = 0 ;
                for( auto & ref : refs )
                {
                    ref.Init() ;
                    chromesome_area.push_back(
                            std::make_tuple(curr_pos ,curr_pos + ref.length -1)) ;
                    curr_pos += ref.length ;
                }
            }
            bool IsValidArea( long start , int length ) const
            {
                auto itr = std::lower_bound( 
                        chromesome_area.begin() 
                        ,chromesome_area.end()
                        ,std::make_tuple(start,100000000L) );
                if ( itr == chromesome_area.end() )
                    return false ;
                long tmp_start ; long tmp_end ;
                std::tie( tmp_start , tmp_end ) 
                    = *itr ;
                if( start + length -1 > tmp_end )
                    return false ;
                return true ;
            }
        };
    }
}
#endif

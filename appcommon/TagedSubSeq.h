#ifndef __APPCOMMON_TAGEDSUBSEQ_H__
#define __APPCOMMON_TAGEDSUBSEQ_H__

#include "appcommon/Ref.h"
#include "biocommon/seq/seq.h"
#include "biocommon/seq/tool_func.h"
namespace BGIQD {
    namespace stLFRSim {

        template<class T>
            struct SubItem
            {
                typedef T BaseType;
                const BaseType & ref ;
                const long long start_pos ;
                const int length ;
                SubItem( const BaseType & r , long long s , int l )
                    : ref(r) 
                      , start_pos(s)
                      , length(l)
                {}

                SubItem( const SubItem& o) : SubItem(o.ref,o.start_pos,o.length) {}
            };


        typedef SubItem<RefChromesome> LongRead;
        typedef SubItem<LongRead>      InsertFragment;

        struct PE
        {
            std::string read1;
            std::string read2;
        };

        PE GetPE(const InsertFragment & IF 
                , int r1_len
                , int r2_len 
                )
        {
            PE ret;
            int r1_start = IF.start_pos 
                + IF.ref.start_pos ;
            int r2_start = r1_start + IF.length -r2_len ;
            const std::string & ref = IF.ref.ref.fa.seq.atcgs;
            ret.read1 =ref.substr(r1_start,r1_len);
            ret.read2 =BGIQD::SEQ::seqCompleteReverse(
                    ref.substr(r2_start,r2_len));

            if( rand() % 2 == 0 )
            {
                std::swap(ret.read1,ret.read2);
            }
            return ret ;
        }
    }
}
#endif // __APPCOMMON_TAGEDSUBSEQ_H__

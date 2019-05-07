#ifndef __APPCOMMON_TAGEDSUBSEQ_H__
#define __APPCOMMON_TAGEDSUBSEQ_H__

#include "appcommon/Ref.h"
namespace BGIQD {
    namespace stLFRSim {

        struct SubSeq
        {
            const RefChromesome & ref ;
            const int start_pos ;
            const int length ;
            SubSeq( const RefChromesome & r , int s , int l )
                : ref(r) 
                , start_pos(s)
                , length(l)
            {}

            SubSeq( const SubSeq & o) : SubSeq(o.ref,o.start_pos,o.length) {}
        };

        struct LongRead : public  SubSeq
        {
            ;
        };

        struct InsertFragment : public SubSeq
        {
            int lr_start_pos ;
        };

    }
}
#endif // __APPCOMMON_TAGEDSUBSEQ_H__

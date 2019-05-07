#ifndef __APPCOMMON_REF_H__
#define __APPCOMMON_REF_H__

#include "biocommon/fasta/fasta.h"

namespace BGIQD {
    namespace stLFRSim {

        typedef BGIQD::FASTA::NormalHead RefHead;

        typedef BGIQD::FASTA::Fasta<RefHead> RefFa;

        struct RefChromesome
        {
            int id ;

            RefFa fa ;

        };
    }
}
#endif

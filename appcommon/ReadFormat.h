#ifndef __APPCOMMON_READFORMAT_H__
#define __APPCOMMON_READFORMAT_H__

#include "appcommon/TagedSubSeq.h"
#include "appcommon/Mutation.h"
#include "biocommon/seq/seq.h"
#include "biocommon/align_common/align_result.h"
#include <iostream>


namespace BGIQD {
    namespace stLFRSim {

        void FormatPrint(
                std::ostream & ost
             ,  long long read_index
             ,  int barcode_id
             ,  int obarcode_id
             ,  int read_id 
             ,  const BGIQD::stLFRSim::InsertFragment & IS
             ,  const BGIQD::Random::MutationResult & read
             )
        {
            ost<<"@stlfrsim_"<<read_index
                <<"#barcode_"<<barcode_id
                <<'/'<<read_id
                <<'\t'<<barcode_id
                <<'\t'<<obarcode_id
                <<'\t'<<IS.ref.ref.fa.head.Id
                <<'\t'<<IS.ref.ref.length
                <<'\t'<<IS.ref.start_pos
                <<'\t'<<IS.ref.length
                <<'\t'<<IS.start_pos
                <<'\t'<<IS.length
                <<'\t'<<BGIQD::ALIGN_COMMON::MatchInfos2CIGAR(read.details)
                <<'\n';
            ost<<read.seq<<'\n';
            ost<<"+\n";
            ost<<std::string(read.seq.size(),'F')<<'\n';
        }
    }
}
#endif 

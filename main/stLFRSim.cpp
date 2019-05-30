#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"

#include "appcommon/TagedSubSeq.h"
#include "appcommon/Pool.h"
#include "appcommon/Random.h"
#include "appcommon/Ref.h"
#include "appcommon/Mutation.h"
#include "appcommon/ReadFormat.h"

#include <string>
#include <sstream>
struct AppConfig
{
    std::string lr_length_file;
    BGIQD::Random::DiscreteRandomWithBin lr_length_dis;
    std::string pe_num_file;
    BGIQD::Random::DiscreteRandomWithBin pe_num_dis;
    std::string pe_length_file;
    BGIQD::Random::DiscreteRandomWithBin pe_length_dis;

    BGIQD::stLFRSim::Ref the_ref ;
    std::string ref_name ;

    std::string o_prefix;
    float max_slr_cov ;
    int read_len ;

    std::string r1_fq ;
    std::string r2_fq ;
    std::ostream * or1;
    std::ostream * or2;

    typedef BGIQD::stLFRSim::InsertFragment IS ;
    typedef BGIQD::stLFRSim::Pool<IS> ReadPairPool;
    ReadPairPool buffer;

    BGIQD::Random::DiscreteRandomWithBin LoadDistributionFromFile( const std::string & file )
    {
        BGIQD::Random::DiscreteRandomWithBin ret ;
        auto fin = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
        if( fin == NULL )
            FATAL(" open distribution file to read failed !!! ");
        auto each_line = [&ret] ( const std::string & line )
        {
            std::istringstream ist(line);
            int s ; int b ; int w ;
            ist >> s>>b>>w;
            ret.keybin.push_back( { s ,b , w } );
            return ret ;
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*fin,each_line);
        return ret ;
    }

    void LoadDistribution()
    {
        lr_length_dis = LoadDistributionFromFile(lr_length_file);
        pe_length_dis = LoadDistributionFromFile(pe_length_file);
        pe_num_dis    = LoadDistributionFromFile(pe_num_file);

        lr_length_dis.InitDistribution();
        pe_length_dis.InitDistribution();
        pe_num_dis.InitDistribution();
    }

    void LoadReference()
    {
        typedef BGIQD::stLFRSim::RefFa Fa;
        typedef BGIQD::FASTA::FastaReader<Fa>  Reader;
        auto fin = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(ref_name);
        if( fin == NULL )
            FATAL(" open ref file to read failed !!! ");
        Fa tmp ;
        while( Reader::LoadNextFasta(*fin , tmp) )
        {
            the_ref.refs.emplace_back(tmp);
        }
        the_ref.Init();
    }


    int RandomLRLengthByDistribution()
    {
        return lr_length_dis();
    }

    bool ValidLR(long long start , int length ) const 
    {
        return the_ref.IsValidArea(start, length );
    }

    BGIQD::stLFRSim::LongRead  GetLR(long long start , int length ) const 
    {
        long long chromesome_start ;
        const auto & chromesome = the_ref.GetChromesome(start, chromesome_start );
        return BGIQD::stLFRSim::LongRead(chromesome ,start-chromesome_start , length );
    }

    int RandomPENumByDistribution() 
    {
        return pe_num_dis();
    }
    bool ValidPENum( int pe_num , int lr_length ) const 
    {
        return float( pe_num * 2 * read_len  ) / float ( lr_length ) < max_slr_cov ;
    }


    int RandomPELengthByDistribution()
    {
        return pe_length_dis();
    }

    bool ValidInsertFragment( BGIQD::stLFRSim::LongRead & lr ,
            int start , int len ) const 
    {
        return lr.ref.IsValidArea(start+lr.start_pos,len)
            && start + len < lr.length ;
    }

    void Init()
    {
        r1_fq = o_prefix+"1.fq";
        r2_fq = o_prefix+"2.fq";
        or1 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(r1_fq);
        or2 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(r2_fq);
        if( NULL ==  or1 )
            FATAL("failed to open o_prefix.1.fq to write !!!");
        if( NULL == or2 )
            FATAL("failed to open o_prefix.2.fq to write !!!");
    }

    long long RefLen() const 
    {
        return the_ref.length ;
    }

    BGIQD::Random::MutationEngine the_mut;


    void PrintReadsFromBuff() 
    {
        static int barcode_id = 1 ;
        static long long read_id = 1; 
        while( buffer.Size() > 0)
        {
            auto IF = buffer.Top();
            auto basic_pe = BGIQD::stLFRSim::GetPE
                ( *IF , read_len,read_len);
            auto r1 = the_mut(basic_pe.read1);
            auto r2 = the_mut(basic_pe.read2);
            FormatPrint(*or1,read_id,barcode_id,1,*IF,r1);
            FormatPrint(*or2,read_id,barcode_id,2,*IF,r2);
            buffer.Pop() ;
            read_id++ ;
        }
        barcode_id ++ ;
    }

    void ClearBuff() 
    {
        buffer.Clear();
    }

    void  AddInsertFragment2Buff( 
            BGIQD::stLFRSim::LongRead & lr , int start , int len ) 
    {
        IS * next = buffer.Push();
        new (next) IS(lr,start,len);
    }

}config;

int main(int argc , char ** argv  )
{
    ////////////////////////////////////////////////////////
    //
    //STEP 1 : Processing input parameters.
    //
    ////////////////////////////////////////////////////////
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , ref
            , "reference fasta file" );
    DEFINE_ARG_REQUIRED(std::string , o_prefix
            , "output file prefix . print into o_prefix.1.fq && o_prefix.2.fq" );
        DEFINE_ARG_REQUIRED(std::string , lr_length_distribution 
                , "distribution file of long read length" );
    DEFINE_ARG_REQUIRED(std::string , pe_num_distribution 
            , "distribution file of number of read-pair in 1 long read" );
    DEFINE_ARG_REQUIRED(std::string , if_lenth_distribution 
            , "distribution file of insert fragment length" );
    DEFINE_ARG_REQUIRED(long long, readpair_num
            , "total number of final generated read-pairs" );

    DEFINE_ARG_OPTIONAL(float , mutation_rate, "mutation rate" , "0.005" );
    DEFINE_ARG_OPTIONAL(float , insert_percent, "insert percent" , "0.005" );
    DEFINE_ARG_OPTIONAL(float , delete_percent, "delete percent " , "0.005" );
    DEFINE_ARG_OPTIONAL(float , substitute_percent, "substitute percent " , "0.99" );
    DEFINE_ARG_OPTIONAL(float , max_slr_cov, "max single long read cov" , "0.5" );
    DEFINE_ARG_OPTIONAL(int   , read_len ,     "read length" , "100" );
    END_PARSE_ARGS ;
    
    config.ref_name = ref.to_string() ;
    config.o_prefix = o_prefix.to_string() ;
    config.lr_length_file = lr_length_distribution.to_string();
    config.pe_num_file = pe_num_distribution.to_string() ;
    config.pe_length_file = if_lenth_distribution.to_string() ;
    config.the_mut.mutation_rate = mutation_rate.to_float();
    config.the_mut.insert_percent = insert_percent.to_float() ;
    config.the_mut.delete_percent = delete_percent.to_float() ;
    config.the_mut.substitute_percent = substitute_percent.to_float() ;
    config.read_len = read_len.to_int() ;
    config.max_slr_cov = max_slr_cov.to_float() ;
    config.Init();

    ////////////////////////////////////////////////////////
    //
    //STEP 2 : Loading distributions.
    //
    ////////////////////////////////////////////////////////
    config.LoadDistribution();
    ////////////////////////////////////////////////////////
    //
    //step 3 : Loading reference.
    //
    ////////////////////////////////////////////////////////
    config.LoadReference();

    ////////////////////////////////////////////////////////
    //
    //step 4 : Produce simulation data
    //
    ////////////////////////////////////////////////////////

    long long succ = 0;
    long long fail = 0;
    long long  R = 0 ;
    while( R <  readpair_num.to_long() )
    {
        // STEP 4.1 . Genarating a long read 
        long long lr_start  =  BGIQD::Random::
            RandomStartPosByLength(config.RefLen());
        int  lr_length = config.RandomLRLengthByDistribution() ;
        if( ! config.ValidLR( lr_start , lr_length ) )
            continue ;
        BGIQD::stLFRSim::LongRead lr
            = config.GetLR( lr_start , lr_length );
        // STEP 4.2 . Random a read pair number.
        int pe_num = config.RandomPENumByDistribution() ;
        if( ! config.ValidPENum( pe_num , lr_length ) )
            continue ;
        // STEP 4.3 . Cyclic generate reads until reads enough or too much failure happened.
        int j = 0 ; int k = 0 ;
        while( j < pe_num && j+k < 2*pe_num )
        {
            int pe_start = BGIQD::Random::
                RandomStartPosByLength(lr_length);
            int pe_length = config.RandomPELengthByDistribution() ;
            if ( ! config.ValidInsertFragment( lr, pe_start, pe_length ) )
            {
                k ++ ;
                continue ;
            }
            else
            {
                j++ ;
                config.AddInsertFragment2Buff(lr, pe_start, pe_length );
            }
        }
        // STEP 4.4 . Output and Log and Clean
        if( j >= pe_num ) // succ
        {
            R += j ;
            config.PrintReadsFromBuff() ;
            succ ++ ;
        }
        else
            fail ++ ;
        config.ClearBuff();
    }

    std::cerr<<" Total succ long read : "<<succ<<'\n';
    std::cerr<<" Total fail long read : "<<fail<<'\n';
    std::cerr<<" Total read pair num  : "<<R<<'\n';
    return 0 ;
}

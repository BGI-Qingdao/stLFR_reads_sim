#include "common/args/argsparser.h"

#include "appcommon/TagedSubSeq.h"
#include "appcommon/Pool.h"

#include <string>

struct AppConfig
{

    void ParseDistribution()
    {

    }

    void LoadReference()
    {

    }

    int ReadNum() const 
    {

    }

    long RandomLRStart() const 
    {

    }
    int RandomLRLengthByDistribution() const 
    {

    }

    bool ValidLR(long start , int length ) const 
    {

    }

    BGIQD::stLFRSim::LongRead  GetLR(long start , int length ) const 
    {

    }

    int RandomPENumByDistribution() const 
    {

    }
    bool ValidPENum( int pe_num , int lr_length ) const 
    {

    }

    int RandomPEStart( int lr_length ) const 
    {

    }

    int RandomPELengthByDistribution() const 
    {

    }

    bool ValidInsertFragment( BGIQD::stLFRSim::LongRead & lr ,
            int start , int len ) const 
    {

    }

    void Init() 
    {

    }

    typedef BGIQD::stLFRSim::Pool<BGIQD::stLFRSim::InsertFragment> ReadPairPool;
    void  AddInsertFragment2Buff( 
            BGIQD::stLFRSim::LongRead & lr , int start , int len ) 
    {

    }
}config;



int main(int argc , char ** argv  )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , ref
            , "reference fasta file" );
        DEFINE_ARG_REQUIRED(std::string , lr_length_distribution 
                , "distribution file of long read length" );
    DEFINE_ARG_REQUIRED(std::string , pe_num_distribution 
            , "distribution file of number of read-pair in 1 long read" );
    DEFINE_ARG_REQUIRED(std::string , is_lenth_distribution 
            , "distribution file of insert size length" );
    DEFINE_ARG_REQUIRED(long , readpair_num
            , "total number of final generated read-pairs" );

    DEFINE_ARG_OPTIONAL(float , mutation_rate, "mutation rate" , "0.005" );
    DEFINE_ARG_OPTIONAL(float , insert_percent, "insert percent" , "0.005" );
    DEFINE_ARG_OPTIONAL(float , delete_percent, "delete percent " , "0.005" );
    DEFINE_ARG_OPTIONAL(float , substitute_percent, "substitute percent " , "0.99" );
    DEFINE_ARG_OPTIONAL(float , max_slr_cov, "max single long read cov" , "0.5" );
    END_PARSE_ARGS ;

    config.Init();
    config.ParseDistribution();
    config.LoadReference();

    long R = 0 ;
    while( R < config.ReadNum() )
    {
        long lr_start  = config.RandomLRStart() ;
        int  lr_length = config.RandomLRLengthByDistribution() ;
        if( ! config.ValidLR( lr_start , lr_length ) )
            continue ;
        int pe_num = config.RandomPENumByDistribution() ;
        if( ! config.ValidPENum( pe_num , lr_length ) )
            continue ;
        BGIQD::stLFRSim::LongRead lr
            = config.GetLR( lr_start , lr_length );
        int j = 0 ; int k = 0 ;
        while( j < pe_num || j+k < 2*pe_num )
        {
            int pe_start = config.RandomPEStart(lr_length);
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
        if( j >= pe_num ) // succ
        {
            config.PrintReadsFromBuff() ;
        }
        else
        {

        }
        config.ClearBuff();
    }

    return 0 ;
}


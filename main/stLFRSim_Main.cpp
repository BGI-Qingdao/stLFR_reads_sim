#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"

#include "appcommon/TagedSubSeq.h"
#include "appcommon/Pool.h"
#include "appcommon/Random.h"
#include "appcommon/Ref.h"

#include <string>
#include <sstream>
struct AppConfig
{

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

    std::string lr_length_file;
    BGIQD::Random::DiscreteRandomWithBin lr_length_dis;
    std::string pe_num_file;
    BGIQD::Random::DiscreteRandomWithBin pe_num_dis;
    std::string pe_length_file;
    BGIQD::Random::DiscreteRandomWithBin pe_length_dis;

    void ParseDistribution()
    {
        lr_length_dis = LoadDistributionFromFile(lr_length_file);
        pe_length_dis = LoadDistributionFromFile(pe_length_file);
        pe_num_dis    = LoadDistributionFromFile(pe_num_file);
    }

    
    void LoadReference()
    {

    }


    int RandomLRLengthByDistribution()
    {
        return lr_length_dis();
    }

    bool ValidLR(long start , int length ) const 
    {
        
    }

    BGIQD::stLFRSim::LongRead  GetLR(long start , int length ) const 
    {

    }

    int RandomPENumByDistribution() 
    {
        return pe_num_dis();
    }
    bool ValidPENum( int pe_num , int lr_length ) const 
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

    long RefLen() const 
    {
        return 1 ;
    }
    void PrintReadsFromBuff() 
    {

    }
    void ClearBuff() 
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

    config.ParseDistribution();
    config.LoadReference();

    config.Init();

    long R = 0 ;
    while( R <  readpair_num.to_long() )
    {
        long lr_start  =  BGIQD::Random::
            RandomStartPosByLength(config.RefLen());
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


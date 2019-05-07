#ifndef __BIOCOMMON_FASTA_FASTA_H__
#define __BIOCOMMON_FASTA_FASTA_H__

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "biocommon/seq/seq.h"
#include "common/flags/flags.h"

namespace BGIQD {

    namespace FASTA {

        struct NormalHead
        {
            std::string head;

            void Init( const std::string & line )
            {
                head = line ;
            }
            std::string Head() const { return head ; }

            void Reset() { head.clear() ;} 
        };

        struct Id_Desc_Head
        {
            std::string Id;

            std::string Desc;

            void Init( const std::string & line ) ;

            std::string Head() const { return ">" + Id +Desc ; }

            void Reset() { Id.clear() ; Desc.clear() ;} 
        };

        template<class T>
            struct Fasta
            {
                typedef T Header;

                FLAGS_INT ;

                ADD_A_FLAG(1,UnSet);
                ADD_A_FLAG(2,Set_head);
                ADD_A_FLAG(3,Set_seq);

                void Reset() 
                {
                    head.Reset();
                    seq.Reset();
                    flags = 0 ;
                    Set_UnSet();
                }

                Header head;
                BGIQD::SEQ::seq seq;
                void AddHead(const std::string & line)
                {
                    Clean_UnSet();
                    Set_Set_head();
                    head.Init(line);
                }

                void AddSeq(const std::string & line )
                {
                    if( Is_UnSet() || ! Is_Set_head() )
                    {
                        assert(0);
                    }
                    Set_Set_seq();
                    seq.AddPartSeq(line);
                }
                bool Is_Setted() const { 
                    return  (!Is_UnSet() )
                        && Is_Set_head() 
                        && Is_Set_seq() ;
                }
            };

        template<class T > 
            struct FastaReader
            {
                typedef T Fasta;

                static bool IsHead(const std::string & line)
                {
                    return ( ! line.empty()) && line[0] == '>' ;
                }

                static void LoadAllFasta( std::istream & ist , std::vector<Fasta> & buffer )
                {
                    std::string line ;
                    Fasta fa;
                    fa.Reset();
                    while ( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            if(fa.Is_Setted())
                            {
                                buffer.push_back(fa);
                            }
                            fa.Reset();
                            fa.AddHead(line);
                        }
                        else
                        {
                            fa.AddSeq(line);
                        }
                    }
                    if(fa.Is_Setted())
                    {
                        buffer.push_back(fa);
                    }
                    fa.Reset();
                }

                static bool LoadNextFasta(std::istream & ist , Fasta & fa)
                {
                    std::string line ;
                    fa.Reset();
                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            fa.AddHead(line);
                            break ;
                        }
                    }
                    if( ist.eof() || ! fa.Is_Set_head() )
                        return false ;

                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            break ;
                        }
                        else
                        {
                            fa.AddSeq(line);
                        }
                    }
                    // Put the head line back into istream
                    if ( ! ist.eof() )
                    {
                        ist.rdbuf()->sputbackc('\n');
                        for( auto  i = line.rbegin() ; i!= line.rend() ; i++ )
                        {
                            ist.rdbuf()->sputbackc(*i);
                        }
                    }

                    return  fa.Is_Setted();
                }
            };
    }
}

#endif //__BIOCOMMON_FASTA_FASTA_H__
.PHONY: all clean bin

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11\
				-I../\
				-lz\
				-lpthread\

DEUBG_CXX  =	${CXXFLAGS} -g
RELEASE_CXX=	${CXXFLAGS}

source_cpp =	../common/files/file_reader.cpp \
				../biocommon/fasta/fasta.cpp\
				../biocommon/fastq/fastq.cpp\
				../biocommon/seq/seq.cpp\
                ../biocommon/align_common/align_result.cpp\
		   		../common/files/file_writer.cpp \
		   		../common/files/gzstream.cpp \
				../common/log/log.cpp\
				../common/log/logfilter.cpp\
				../common/time/timetools.cpp\
				../common/string/stringtools.cpp\
				../common/args/argsparser.cpp\

source_o		= ${source_cpp:%.cpp=%.o}

.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs = stLFRSim

all :  ${jobs}




stLFRSim_cpp = stLFRSim.cpp
stLFRSim_o = stLFRSim.o
stLFRSim : clean ${stLFRSim_o} ${source_o} ../bin
	${CXX} ${stLFRSim_o} ${source_o} ${DEUBG_CXX} -o stLFRSim
	mv $@ ../bin/


jobs_o=  ${stLFRSim_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}


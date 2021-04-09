/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date:   Mar, 2021
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <zlib.h>
#include "common.h"
#include "util.h"
using namespace std;

void inline CSEREAD_resize( CSEREAD * cr, int n ) {
	cr->seq[ n] = 0;
	cr->qual[n] = 0;
	cr->size    = n;
}

void find_seed( vector<unsigned int> &seed, CSEREAD *read, const ktrim_param & kp ) {
	seed.clear();
	register char *poffset  = read->seq;
	register char *indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index1 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		indexloc ++;
	}
	poffset  = read->seq + OFFSET_INDEX3;
	indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index3 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		indexloc ++;
	}
	sort( seed.begin(), seed.end() );
}

void workingThread_SE_C( unsigned int tn, unsigned int start, unsigned int end, CSEREAD *workingReads,
							ktrim_stat * kstat, writeBuffer * writebuffer, const ktrim_param & kp ) {

//	fprintf( stderr, "=== working thread %d: %d - %d\n", tn, start, end ), "\n";

	writebuffer->b1stored[tn] = 0;

	register int i, j;
	register unsigned int last_seed;
	vector<unsigned int> seed;
	vector<unsigned int> :: iterator it;
	const char *p, *q;

	register CSEREAD * wkr = workingReads + start;
	for( unsigned int ii=start; ii!=end; ++ii, ++wkr ) {
//		fprintf( stderr, "working: %d, %s\n", ii, wkr->id );
		// quality control
		p = wkr->qual;
		j = wkr->size;

		// update in v1.2: support window check
		i = get_quality_trim_cycle_se( p, j, kp );

		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != j ) {  // quality-trim occurs
			CSEREAD_resize( wkr, i);
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		find_seed( seed, wkr, kp );

		last_seed = impossible_seed;	// a position which cannot be in seed
		for( it=seed.begin(); it!=seed.end(); ++it ) {
			if( *it != last_seed ) {
//				fprintf( stderr, " check seed: %d\n", *it );
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
				if( check_mismatch_dynamic_SE_C( wkr, *it, kp ) )
					break;
		
				last_seed = *it;
			}
		}
		if( it != seed.end() ) {	// adapter found
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				CSEREAD_resize( wkr, *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];

				if( *it <= DIMER_INSERT )
				  ++ kstat->dimer[tn];

				continue;
			}
		} else {	// seed not found, now check the tail 2, if perfect match, drop these 2; Single-end reads do not check tail 1
			i = wkr->size - 2;
			p = wkr->seq;
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] ) {
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				CSEREAD_resize( wkr, i );
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
												"%s%s\n+\n%s\n", wkr->id, wkr->seq, wkr->qual);

	}
}

int process_multi_thread_SE_C( const ktrim_param &kp ) {
	// IO speed-up
//	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	// in this version, two data containers are used and auto-swapped for working and loading data
	CSEREAD *readA = new CSEREAD[ READS_PER_BATCH ];
	CSEREAD *readB = new CSEREAD[ READS_PER_BATCH ];

	register char *readA_data = new char[ MEM_SE_READSET ];
	register char *readB_data = new char[ MEM_SE_READSET ];

	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		readA[i].id   = readA_data + j;
		readB[i].id   = readB_data + j;
		j += MAX_READ_ID;
		readA[i].seq  = readA_data + j;
		readB[i].seq  = readB_data + j;
		j += MAX_READ_CYCLE;
		readA[i].qual = readA_data + j;
		readB[i].qual = readB_data + j;
		j += MAX_READ_CYCLE;
	}

	CSEREAD *workingReads, *loadingReads, *swapReads;

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ kp.thread ];
	kstat.real_adapter = new unsigned int [ kp.thread ];
	kstat.tail_adapter = new unsigned int [ kp.thread ];
	kstat.dimer	       = new unsigned int [ kp.thread ];

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ kp.thread ];
	writebuffer.b1stored = new unsigned int	[ kp.thread ];

	for(unsigned int i=0; i!=kp.thread; ++i) {
		writebuffer.buffer1[i]  = new char[ BUFFER_SIZE_PER_BATCH_READ ];
		writebuffer.b1stored[i] = 0;

		kstat.dropped[i] = 0;
		kstat.real_adapter[i] = 0;
		kstat.tail_adapter[i] = 0;
		kstat.dimer[i] = 0;
	}

	// deal with multiple input files
	vector<string> R1s;
	extractFileNames( kp.FASTQU, R1s );
	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " single-end fastq files will be loaded.\033[0m\n";

	FILE *fout1;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		return 103;
	}

	register unsigned int line = 0;
	unsigned int threadCNT = kp.thread - 1;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq;
		gzFile gfp;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp = gzopen( p, "r" );
			if( gfp == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		} else {
			fq = fopen( p, "rt" );
			if( fq == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		}
		// initialization
		// get first batch of fastq reads
		unsigned int loaded;
		bool metEOF;
		if( file_is_gz ) {
			loaded = load_batch_data_SE_GZ( gfp, readA, READS_PER_BATCH );
			metEOF = gzeof( gfp );
		} else {
			loaded = load_batch_data_SE_C( fq, readA, READS_PER_BATCH );
			metEOF = feof( fq );
		}
		if( loaded == 0 ) break;
//		fprintf( stderr, "Loaded %d, metEOF=%d\n", loaded, metEOF );

		loadingReads = readB;
		workingReads = readA;
		bool nextBatch = true;
		unsigned int threadLoaded;
		unsigned int NumWkThreads;
		while( nextBatch ) {
			// start parallalization
			omp_set_num_threads( kp.thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				// if EOF is met, then all threads are used for analysis
				// otherwise 1 thread will do data loading
				if( metEOF ) {
					NumWkThreads = kp.thread;
					unsigned int start = loaded * tn / kp.thread;
					unsigned int end   = loaded * (tn+1) / kp.thread;
					workingThread_SE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					nextBatch = false;
				} else {	// use 1 thread to load file, others for trimming
					NumWkThreads = threadCNT;
					if( tn == threadCNT ) {
						if( file_is_gz ) {
							threadLoaded = load_batch_data_SE_GZ( gfp, loadingReads, READS_PER_BATCH );
							metEOF = gzeof( gfp );
						} else {
							threadLoaded = load_batch_data_SE_C( fq, loadingReads, READS_PER_BATCH );
							metEOF = feof( fq );
						}
						nextBatch = (threadLoaded!=0);
//						fprintf( stderr, "Loaded %d, metEOF=%d\n", threadLoaded, metEOF );
						//cerr << "Loading thread: " << threadLoaded << ", " << metEOF << ", " << nextBatch << '\n';
					} else {
						unsigned int start = loaded * tn / threadCNT;
						unsigned int end   = loaded * (tn+1) / threadCNT;
						workingThread_SE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					}
				}
			} // parallel body
			// swap workingReads and loadingReads for next loop
			swapReads	= loadingReads;
			loadingReads = workingReads;
			workingReads = swapReads;
			// write output and update fastq statistics
			for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
				fwrite( writebuffer.buffer1[ii], sizeof(char), writebuffer.b1stored[ii], fout1 );
			}
			line += loaded;
			loaded = threadLoaded;
			//cerr << '\r' << line << " reads loaded";
		}

		if( file_is_gz ) {
			gzclose( gfp );
		} else {
			fclose( fq );
		}
	}

	fclose( fout1 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	int dropped_all=0, real_all=0, tail_all=0, dimer_all=0;
	for( unsigned int i=0; i!=kp.thread; ++i ) {
		dropped_all += kstat.dropped[i];
		real_all  += kstat.real_adapter[i];
		tail_all  += kstat.tail_adapter[i];
		dimer_all += kstat.dimer[i];
	}
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}
	fout << "Total\t"	 << line		<< '\n'
		 << "Dropped\t"  << dropped_all << '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t"  << tail_all	<< '\n'
		 << "Dimer\t"    << dimer_all	<< '\n';
	fout.close();

	//free memory
	for(unsigned int i=0; i!=kp.thread; ++i) {
		delete writebuffer.buffer1[i];
	}
	delete [] writebuffer.buffer1;
	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;

	delete [] readA;
	delete [] readB;
	delete [] readA_data;
	delete [] readB_data;

	return 0;
}

int process_single_thread_SE_C( const ktrim_param &kp ) {
//	fprintf( stderr, "SINGLE END SINGLE_THREAD MODE ON\n" );
	// IO speed-up
//	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	CSEREAD *read = new CSEREAD[ READS_PER_BATCH_ST ];
	register char *read_data = new char[ MEM_SE_READSET ];

	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		read[i].id   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual = read_data + j;
		j += MAX_READ_CYCLE;
	}

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ 1 ];
	kstat.real_adapter = new unsigned int [ 1 ];
	kstat.tail_adapter = new unsigned int [ 1 ];
	kstat.dimer        = new unsigned int [ 1 ];
	kstat.dropped[0] = 0;
	kstat.real_adapter[0] = 0;
	kstat.tail_adapter[0] = 0;
	kstat.dimer[0] = 0;

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ 1 ];
	writebuffer.b1stored = new unsigned int	[ 1 ];
	writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

	// deal with multiple input files
	vector<string> R1s;
	extractFileNames( kp.FASTQU, R1s );

	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " single-end fastq files will be loaded.\033[0m\n";

	FILE *fout1;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		return 103;
	}

	//ifstream fq1;
	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		//fq1.open( R1s[fileCnt].c_str() );
		bool file_is_gz = false;
		FILE *fq;
		gzFile gfp;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
//			fprintf( stderr, "GZ file!\n" );
			file_is_gz = true;
			gfp = gzopen( p, "r" );
			if( gfp == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		} else {
			fq = fopen( p, "rt" );
			if( fq == NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				return 104;
			}
		}

		register unsigned int last_seed;
		vector<unsigned int> seed;
		vector<unsigned int> :: iterator it;

		while( true ) {
			// get fastq reads
			//unsigned int loaded = load_batch_data_SE( fq1, read, READS_PER_BATCH_ST );
			unsigned int loaded;
			if( file_is_gz ) {
				loaded = load_batch_data_SE_GZ( gfp, read, READS_PER_BATCH_ST );
//				fprintf( stderr, "Loaded=%d\n", loaded );
			} else {
				loaded = load_batch_data_SE_C( fq, read, READS_PER_BATCH_ST );
			}

			if( loaded == 0 )
				break;
		
//			fprintf( stderr, "Work\n" );
			workingThread_SE_C( 0, 0, loaded, read, &kstat, &writebuffer, kp );

			// write output and update fastq statistics
//			fprintf( stderr, "Output\n" );
			fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], fout1 );

			line += loaded;
			//cerr << '\r' << line << " reads loaded";

			//if( fq1.eof() ) break;
//			fprintf( stderr, "Check\n" );
			if( file_is_gz ) {
				if( gzeof( gfp ) ) break;
			} else {
				if( feof( fq ) ) break;
			}
		}
		//fq1.close();
		if( file_is_gz ) {
			gzclose( gfp );
		} else {
			fclose( fq );
		}
	}
	fclose( fout1 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}

	fout << "Total\t"    << line					<< '\n'
		 << "Dropped\t"  << kstat.dropped[0]		<< '\n'
		 << "Aadaptor\t" << kstat.real_adapter[0]	<< '\n'
		 << "TailHit\t"  << kstat.tail_adapter[0]	<< '\n'
		 << "Dimer\t"    << kstat.dimer[0]			<< '\n';
	fout.close();

	//free memory
//	delete buffer1;
//	delete buffer2;
//	delete fq1buffer;
	delete [] read;
	delete [] read_data;

	return 0;
}


/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: May 2023
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include "common.h"
using namespace std;

void inline CPEREAD_resize( CPEREAD * read, int n ) {
	read->size = n;
	read->seq1[  n ] = 0;
	read->qual1[ n ] = 0;
	read->seq2[  n ] = 0;
	read->qual2[ n ] = 0;
}

bool inline is_revcomp( const char a, const char b ) {
	//TODO: consider how to deal with N, call it positive or negative???
	switch( a ) {
		case 'A': return b=='T';
		case 'C': return b=='G';
		case 'G': return b=='C';
		case 'T': return b=='A';
		default : return false;
	}
	return true;
}

void init_kstat_wrbuffer( ktrim_stat &kstat, writeBuffer &writebuffer, unsigned int nthread, bool write2stdout ) {
	kstat.dropped	   = new unsigned int [ nthread ];
	kstat.real_adapter = new unsigned int [ nthread ];
	kstat.tail_adapter = new unsigned int [ nthread ];
	kstat.dimer        = new unsigned int [ nthread ];
	kstat.pass         = new unsigned int [ nthread ];

	// buffer for storing the modified reads per thread
	writebuffer.buffer1  = new char * [ nthread ];
	writebuffer.buffer2  = new char * [ nthread ];
	writebuffer.b1stored = new unsigned int	[ nthread ];
	writebuffer.b2stored = new unsigned int [ nthread ];

	for(unsigned int i=0; i!=nthread; ++i) {
		if( write2stdout ) {	// only use buffer 1
			writebuffer.buffer1[i]  = new char[ BUFFER_SIZE_PER_BATCH_READ << 1 ];
			writebuffer.b1stored[i] = 0;
		} else {
			writebuffer.buffer1[i]  = new char[ BUFFER_SIZE_PER_BATCH_READ ];
			writebuffer.buffer2[i]  = new char[ BUFFER_SIZE_PER_BATCH_READ ];
			writebuffer.b1stored[i] = 0;
			writebuffer.b2stored[i] = 0;
		}

		kstat.dropped[i] = 0;
		kstat.real_adapter[i] = 0;
		kstat.tail_adapter[i] = 0;
		kstat.dimer[i] = 0;
		kstat.pass[i]  = 0;
	}
}

void find_seed_pe( vector<unsigned int> &seed, const CPEREAD *read, const ktrim_param & kp ) {
	seed.clear();
	register const char *poffset  = read->seq1;
	register const char *indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index1 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		++ indexloc;
	}
	poffset  = read->seq2;
	indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index2 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		++ indexloc;
	}
	sort( seed.begin(), seed.end() );
}

// this function is slower than C++ version
void workingThread_PE_C( unsigned int tn, unsigned int start, unsigned int end, CPEREAD *workingReads,
					ktrim_stat *kstat, writeBuffer *writebuffer, const ktrim_param &kp ) {

	writebuffer->b1stored[tn] = 0;
	writebuffer->b2stored[tn] = 0;

//	vector<unsigned int> seed;
//	vector<unsigned int> :: iterator it, end_of_seed;
	register int *seed = new int[ MAX_SEED_NUM ];
	register int hit_seed;
	register int *it, *end_of_seed;

	register CPEREAD *wkr = workingReads + start;
	for( unsigned int iii=end-start; iii; --iii, ++wkr ) {
		// read size handling
		if( wkr->size > wkr->size2 )
			wkr->size = wkr->size2;

		// remove the tail '\n'
		// in fact, it is not essential to do this step, because '\n' has a very low ascii value (10)
		// therefore it will be quality-trimmed
		CPEREAD_resize( wkr, wkr->size - 1 );

		// quality control
		register int i = get_quality_trim_cycle_pe( wkr, kp );
		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != wkr->size ) {  // quality-trim occurs
			CPEREAD_resize( wkr, i );
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		//find_seed_pe( seed, wkr, kp );
		//TODO: I donot need to find all the seeds, I can find-check, then next
//		seed.clear();
		hit_seed = 0;
		register const char *poffset  = wkr->seq1;
		register const char *indexloc = poffset;
		while( true ) {
			indexloc = strstr( indexloc, kp.adapter_index1 );
			if( indexloc == NULL )
				break;
			//seed.push_back( indexloc - poffset );
			seed[ hit_seed++ ] = indexloc - poffset;
			++ indexloc;
		}
		poffset  = wkr->seq2;
		indexloc = poffset;
		while( true ) {
			indexloc = strstr( indexloc, kp.adapter_index2 );
			if( indexloc == NULL )
				break;
			//seed.push_back( indexloc - poffset );
			seed[ hit_seed++ ] = indexloc - poffset;
			++ indexloc;
		}
		//sort( seed.begin(), seed.end() );
		end_of_seed = seed + hit_seed;
		if( hit_seed != 0 )
			sort( seed, seed + hit_seed );

		register bool no_valid_adapter = true;
		register unsigned int last_seed = impossible_seed;	// a position which cannot be a seed
		//end_of_seed = seed.end();
		//for( it=seed.begin(); it!=end_of_seed; ++it ) {
		for( it=seed; it!=end_of_seed; ++it ) {
			if( *it != last_seed ) {
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
				if( check_mismatch_dynamic_PE_C( wkr, *it, kp ) )
					break;

				last_seed = *it;
			}
		}
		if( it != end_of_seed ) {	// adapter found
			no_valid_adapter = false;
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				CPEREAD_resize( wkr, *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];

				if( *it <= DIMER_INSERT )
					++ kstat->dimer[tn];
				continue;
			}
		} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
			// note: I will NOT consider tail hits as adapter_found, as '-w' should be used when insertDNA is very short
			i = wkr->size - 2;
			register const char *p = wkr->seq1;
			register const char *q = wkr->seq2;
			//Note: 1 mismatch is allowed in tail-checking
			register unsigned int mismatches = 0;
			if( p[i]   != kp.adapter_r1[0] ) mismatches ++;
			if( p[i+1] != kp.adapter_r1[1] ) mismatches ++;
			if( q[i]   != kp.adapter_r2[0] ) mismatches ++;
			if( q[i+1] != kp.adapter_r2[1] ) mismatches ++;
			if( ! is_revcomp(p[5], q[i-6]) ) mismatches ++;
			if( ! is_revcomp(q[5], p[i-6]) ) mismatches ++;
			if( mismatches <= 1 ) {	// tail is good
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				CPEREAD_resize( wkr, i );
			} else {	// tail 2 is not good, check tail 1
				++ i;
				mismatches = 0;
				if( p[i] != kp.adapter_r1[0] ) mismatches ++;
				if( q[i] != kp.adapter_r2[0] ) mismatches ++;
				if( ! is_revcomp(p[5], q[i-6]) ) mismatches ++;
				if( ! is_revcomp(q[5], p[i-6]) ) mismatches ++;
				if( ! is_revcomp(p[6], q[i-7]) ) mismatches ++;
				if( ! is_revcomp(q[6], p[i-7]) ) mismatches ++;

				if( mismatches <= 1 ) {
					++ kstat->tail_adapter[tn];
					if( i < kp.min_length ) {
						++ kstat->dropped[tn];
						continue;
					}
					CPEREAD_resize( wkr, i );
				}
			}
		}

		if( kp.outputReadWithAdaptorOnly && no_valid_adapter )
			continue;

		++ kstat->pass[tn];
		if( kp.write2stdout ) {
			writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
												"%s%s\n+\n%s\n%s%s\n+\n%s\n",
												wkr->id1, wkr->seq1, wkr->qual1,
												wkr->id2, wkr->seq2, wkr->qual2);
		} else {
			writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
												"%s%s\n+\n%s\n", wkr->id1, wkr->seq1, wkr->qual1 );
			writebuffer->b2stored[tn] += sprintf( writebuffer->buffer2[tn]+writebuffer->b2stored[tn],
												"%s%s\n+\n%s\n", wkr->id2, wkr->seq2, wkr->qual2 );
		}
	}

	//wait for my turn to output
	while( true ) {
		if( tn == write_thread ) {
			if( kp.write2stdout ) {
				fwrite( writebuffer->buffer1[tn], sizeof(char), writebuffer->b1stored[tn], stdout );
			} else {
				fwrite( writebuffer->buffer1[tn], sizeof(char), writebuffer->b1stored[tn], kp.fout1 );
				fwrite( writebuffer->buffer2[tn], sizeof(char), writebuffer->b2stored[tn], kp.fout2 );
			}
			++ write_thread;
			break;
		} else {
			this_thread::sleep_for( waiting_time_for_writing );
		}
	}
	
	delete [] seed;
}

int process_multi_thread_PE_C( const ktrim_param &kp ) {
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	// in this version, two data containers are used and auto-swapped for working and loading data
	CPEREAD *readA, *readB;
	register char *readA_data, *readB_data;
	CPEREAD *workingReads, *loadingReads, *swapReads;
	ktrim_stat kstat;
	writeBuffer writebuffer;
//	vector<string> R1s, kp.R2s;
	unsigned int totalFiles;

// now I use 2 threads for init
// prepare memory and file
	omp_set_num_threads( 2 );
	#pragma omp parallel
	{
		unsigned int tn = omp_get_thread_num();

		if( tn == 0 ) {
			readA = new CPEREAD[ READS_PER_BATCH ];
			readB = new CPEREAD[ READS_PER_BATCH ];
			readA_data = new char[ MEM_PE_READSET ];
			readB_data = new char[ MEM_PE_READSET ];

			for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
				readA[i].id1   = readA_data + j;
				readB[i].id1   = readB_data + j;
				j += MAX_READ_ID;
				readA[i].seq1  = readA_data + j;
				readB[i].seq1  = readB_data + j;
				j += MAX_READ_CYCLE;
				readA[i].qual1 = readA_data + j;
				readB[i].qual1 = readB_data + j;
				j += MAX_READ_CYCLE;

				readA[i].id2   = readA_data + j;
				readB[i].id2   = readB_data + j;
				j += MAX_READ_ID;
				readA[i].seq2  = readA_data + j;
				readB[i].seq2  = readB_data + j;
				j += MAX_READ_CYCLE;
				readA[i].qual2 = readA_data + j;
				readB[i].qual2 = readB_data + j;
				j += MAX_READ_CYCLE;
			}
		} else {
			init_kstat_wrbuffer( kstat, writebuffer, kp.thread, kp.write2stdout );

			// deal with multiple input files
			totalFiles = kp.R1s.size();
			//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";
		}
	}

// start analysis
	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = kp.R1s[fileCnt].size() - 3;
		register const char * p = kp.R1s[fileCnt].c_str();
		register const char * q = kp.R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		}
//		fprintf( stderr, "Loading files:\nRead1: %s\nRead2: %s\n", p, q );

		// initialization
		// get first batch of fastq reads
		unsigned int loaded;
		bool metEOF;
		omp_set_num_threads( 2 );
		#pragma omp parallel
		{
			unsigned int tn = omp_get_thread_num();
			if( tn == 0 ) {
				if( file_is_gz ) {
					loaded = load_batch_data_PE_GZ( gfp1, readA, READS_PER_BATCH, true );
					metEOF = gzeof( gfp1 );
				} else {
					loaded = load_batch_data_PE_C( fq1, readA, READS_PER_BATCH, true );
					metEOF = feof( fq1 );
				}
			} else {
				if( file_is_gz ) {
					loaded = load_batch_data_PE_GZ( gfp2, readA, READS_PER_BATCH, false );
				} else {
					loaded = load_batch_data_PE_C( fq2, readA, READS_PER_BATCH, false );
				}
			}
		}
		if( loaded == 0 ) break;

		loadingReads = readB;
		workingReads = readA;
		bool nextBatch = true;
		unsigned int threadLoaded=0, threadLoaded2=0;
		unsigned int NumWkThreads=0;
		while( nextBatch ) {
//			cerr << "Working on " << loaded << " reads\n";
			// start parallalization
			write_thread = 0;
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
					workingThread_PE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					nextBatch = false;
				} else {	// use 2 thread to load files, others for trimming
					NumWkThreads = kp.thread - 2;
					if( tn == kp.thread - 1 ) {
						if( file_is_gz ) {
							threadLoaded = load_batch_data_PE_GZ( gfp1, loadingReads, READS_PER_BATCH, true );
							metEOF = gzeof( gfp1 );
						} else {
							threadLoaded = load_batch_data_PE_C( fq1, loadingReads, READS_PER_BATCH, true );
							metEOF = feof( fq1 );
						}
//		cerr << "R1 loaded " << threadLoaded << ", pos=" << gztell(gfp2) << ", EOF=" << gzeof( gfp1 ) << "\n";
						nextBatch = (threadLoaded!=0);
				//cerr << "Loading thread: " << threadLoaded << ", " << metEOF << ", " << nextBatch << '\n';
					} else if ( tn == kp.thread - 2 ) {
						if( file_is_gz ) {
							threadLoaded2 = load_batch_data_PE_GZ( gfp2, loadingReads, READS_PER_BATCH, false );
						} else {
							threadLoaded2 = load_batch_data_PE_C( fq2, loadingReads, READS_PER_BATCH, false );
						}
//		cerr << "R2 loaded " << threadLoaded2 << ", pos=" << gztell(gfp2) << ", EOF=" << gzeof( gfp2 )<< "\n";
					} else {
						unsigned int start = loaded * tn / NumWkThreads;
						unsigned int end   = loaded * (tn+1) / NumWkThreads;
						workingThread_PE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					}
				}
			} // parallel body
			// write output and update fastq statistics
			// I cannot write fastq in each thread because it may cause unpaired reads
			/*if( ! kp.write2stdout ) {
				omp_set_num_threads( 2 );
				#pragma omp parallel
				{
					unsigned int tn = omp_get_thread_num();
					if( tn == 0 ) {
						for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
							fwrite( writebuffer.buffer1[ii], sizeof(char), writebuffer.b1stored[ii], kp.fout1 );
						}
					} else {
						for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
							fwrite( writebuffer.buffer2[ii], sizeof(char), writebuffer.b2stored[ii], kp.fout2 );
						}
					}
				}
			}*/
			// check whether the read-loading is correct
			if( threadLoaded != threadLoaded2 ) {
				cerr << "ERROR: unequal read number for read1 (" << threadLoaded << ") and read2 (" << threadLoaded2 << ")!\n";
				return 1;
			}

			line += loaded;
			loaded = threadLoaded;
			// swap workingReads and loadingReads for next loop
			swapReads	 = loadingReads;
			loadingReads = workingReads;
			workingReads = swapReads;

//			cerr << '\r' << line << " reads loaded\n";
//			cerr << line << " reads loaded, metEOF=" << metEOF << ", next=" << nextBatch << "\n";
		}//process 1 file

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fclose( fq1 );
			fclose( fq2 );
		}
//		cerr << '\n';
	} // all input files are loaded
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	int dropped_all=0, real_all=0, tail_all=0, dimer_all=0, pass_all=0;
	for( unsigned int i=0; i!=kp.thread; ++i ) {
		dropped_all += kstat.dropped[i];
		real_all  += kstat.real_adapter[i];
		tail_all  += kstat.tail_adapter[i];
		dimer_all += kstat.dimer[i];
		pass_all  += kstat.pass[i];
	}
	fprintf( kp.flog, "Total\t%u\nDropped\t%u\nAadaptor\t%u\nTailHit\t%u\nDimer\t%u\nPass\t%u\n",
				line, dropped_all, real_all, tail_all, dimer_all, pass_all );

	//free memory
	for(unsigned int i=0; i!=kp.thread; ++i) {
		delete writebuffer.buffer1[i];
		if( ! kp.write2stdout )
			delete writebuffer.buffer2[i];
	}
	delete [] writebuffer.buffer1;
	delete [] writebuffer.buffer2;

	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;
	delete [] kstat.dimer;
	delete [] kstat.pass;

	delete [] readA;
	delete [] readB;
	delete [] readA_data;
	delete [] readB_data;

	return 0;
}

int process_single_thread_PE_C( const ktrim_param &kp ) {
//	fprintf( stderr, "process_single_thread_PE_C\n" );
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	CPEREAD *read = new CPEREAD[ READS_PER_BATCH_ST ];
	register char *read_data = new char[ MEM_PE_READSET_ST ];
	
	for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
		read[i].id1   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq1  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual1 = read_data + j;
		j += MAX_READ_CYCLE;
		
		read[i].id2   = read_data + j;
		j += MAX_READ_ID;
		read[i].seq2  = read_data + j;
		j += MAX_READ_CYCLE;
		read[i].qual2 = read_data + j;
		j += MAX_READ_CYCLE;
	}

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ 1 ];
	kstat.real_adapter = new unsigned int [ 1 ];
	kstat.tail_adapter = new unsigned int [ 1 ];
	kstat.dimer        = new unsigned int [ 1 ];
	kstat.pass         = new unsigned int [ 1 ];
	kstat.dropped[0] = 0;
	kstat.real_adapter[0] = 0;
	kstat.tail_adapter[0] = 0;
	kstat.dimer[0] = 0;
	kstat.pass[0]  = 0;

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ 1 ];
	if( kp.write2stdout ) {
		writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST << 1 ];
		writebuffer.b1stored = new unsigned int	[ 1 ];
	} else {
		writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

		writebuffer.buffer2  = new char * [ 1 ];
		writebuffer.buffer2[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];
	}
	writebuffer.b1stored = new unsigned int	[ 1 ];
	writebuffer.b2stored = new unsigned int	[ 1 ];

// deal with multiple input files
	if( kp.R1s.size() != kp.R2s.size() ) {
		cerr << "\033[1;31mError: incorrect files!\033[0m\n";
		return 110;
	}
	unsigned int totalFiles = kp.R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = kp.R1s[fileCnt].size() - 3;
		register const char * p = kp.R1s[fileCnt].c_str();
		register const char * q = kp.R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		}

		while( true ) {
			// get fastq reads
			unsigned int loaded;
			if( file_is_gz ) {
				loaded = load_batch_data_PE_both_GZ( gfp1, gfp2, read, READS_PER_BATCH_ST );
			} else {
				loaded = load_batch_data_PE_both_C( fq1, fq2, read, READS_PER_BATCH_ST );
			}
			if( loaded == 0 ) break;

			write_thread = 0;
			workingThread_PE_C( 0, 0, loaded, read, &kstat, &writebuffer, kp );
			// write output and update fastq statistics
/*			if( ! kp.write2stdout ) {
				fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], kp.fout1 );
				fwrite( writebuffer.buffer2[0], sizeof(char), writebuffer.b2stored[0], kp.fout2 );
			}*/

			line += loaded;
			//cerr << '\r' << line << " reads loaded";

			if( file_is_gz ) {
				if( gzeof( gfp1 ) ) break;
			} else {
				if( feof( fq2 ) ) break;
			}
		}

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fclose( fq1 );
			fclose( fq2 );
		}
	}
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fprintf( kp.flog, "Total\t%u\nDropped\t%u\nAadaptor\t%u\nTailHit\t%u\nDimer\t%u\nPass\t%u\n",
				line, kstat.dropped[0], kstat.real_adapter[0], kstat.tail_adapter[0], kstat.dimer[0], kstat.pass[0] );

	delete writebuffer.buffer1[0];
	if( ! kp.write2stdout ) delete writebuffer.buffer2[0];
	delete [] read;
	delete [] read_data;

	return 0;
}

int process_two_thread_PE_C( const ktrim_param &kp ) {
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	// in this version, two data containers are used and auto-swapped for working and loading data
	CPEREAD *readA;
	register char *readA_data;
	ktrim_stat kstat;
	writeBuffer writebuffer;
//	vector<string> R1s, R2s;
	unsigned int totalFiles;

// now I use 2 theads for init
// prepare memory and file
	omp_set_num_threads( 2 );
	#pragma omp parallel
	{
		unsigned int tn = omp_get_thread_num();

		if( tn == 0 ) {
			readA = new CPEREAD[ READS_PER_BATCH ];
			readA_data = new char[ MEM_PE_READSET ];

			for( register int i=0, j=0; i!=READS_PER_BATCH; ++i ) {
				readA[i].id1   = readA_data + j;
				j += MAX_READ_ID;
				readA[i].seq1  = readA_data + j;
				j += MAX_READ_CYCLE;
				readA[i].qual1 = readA_data + j;
				j += MAX_READ_CYCLE;

				readA[i].id2   = readA_data + j;
				j += MAX_READ_ID;
				readA[i].seq2  = readA_data + j;
				j += MAX_READ_CYCLE;
				readA[i].qual2 = readA_data + j;
				j += MAX_READ_CYCLE;
			}
		} else {
			init_kstat_wrbuffer( kstat, writebuffer, kp.thread, kp.write2stdout );

			// deal with multiple input files
			totalFiles = kp.R1s.size();
			//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";
		}
	}

	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = kp.R1s[fileCnt].size() - 3;
		register const char * p = kp.R1s[fileCnt].c_str();
		register const char * q = kp.R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				cerr << "\033[1;31mError: open fastq file failed!\033[0m\n";
				return 104;
			}
		}
//		fprintf( stderr, "Loading files:\nRead1: %s\nRead2: %s\n", p, q );

		// initialization
		// get first batch of fastq reads

		// there is NO need to do rotation in 2-threads
		unsigned int loaded1, loaded2;
		bool metEOF;
		while( true ) {
			// load data
			omp_set_num_threads( 2 );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				if( tn == 0 ) {
					if( file_is_gz ) {
						loaded1 = load_batch_data_PE_GZ( gfp1, readA, READS_PER_BATCH, true );
						metEOF = gzeof( gfp1 );
					} else {
						loaded1 = load_batch_data_PE_C( fq1, readA, READS_PER_BATCH, true );
						metEOF = feof( fq1 );
					}
				} else {
					if( file_is_gz ) {
						loaded2 = load_batch_data_PE_GZ( gfp2, readA, READS_PER_BATCH, false );
					} else {
						loaded2 = load_batch_data_PE_C( fq2, readA, READS_PER_BATCH, false );
					}
				}
			}
			if( loaded1 != loaded2 ) {
				cerr << "ERROR: unequal read number for read1 (" << loaded1 << ") and read2 (" << loaded2 << ")!\n";
				return 1;
			}
			if( loaded1 == 0 ) break;

			// start analysis
			write_thread = 0;
			omp_set_num_threads( 2 );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				register int middle = loaded1 >> 1;
				if( tn == 0 ) {
					workingThread_PE_C( 0,      0,  middle, readA, &kstat, &writebuffer, kp );
				} else {
					workingThread_PE_C( 1, middle, loaded1, readA, &kstat, &writebuffer, kp );
				}
			} // parallel body
			// write output and update fastq statistics
			/*if( ! kp.write2stdout ) {
				fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], kp.fout1 );
				fwrite( writebuffer.buffer1[1], sizeof(char), writebuffer.b1stored[1], kp.fout1 );
				fwrite( writebuffer.buffer2[0], sizeof(char), writebuffer.b2stored[0], kp.fout2 );
				fwrite( writebuffer.buffer2[1], sizeof(char), writebuffer.b2stored[1], kp.fout2 );
			}*/
			line += loaded1;
			if( metEOF )break;
		}

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fclose( fq1 );
			fclose( fq2 );
		}
	} // all input files are loaded
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fprintf( kp.flog, "Total\t%u\nDropped\t%u\nAadaptor\t%u\nTailHit\t%u\nDimer\t%u\nPass\t%u\n",
				line, kstat.dropped[0]+kstat.dropped[1], kstat.real_adapter[0]+kstat.real_adapter[1],
				kstat.tail_adapter[0]+kstat.tail_adapter[1], kstat.dimer[0]+kstat.dimer[1], kstat.pass[0]+kstat.pass[1] );

	//free memory
	for(unsigned int i=0; i!=kp.thread; ++i) {
		delete writebuffer.buffer1[i];
		delete writebuffer.buffer2[i];
	}
	delete [] writebuffer.buffer1;
	delete [] writebuffer.buffer2;
	delete [] kstat.dropped;
	delete [] kstat.real_adapter;
	delete [] kstat.tail_adapter;
	delete [] kstat.dimer;
	delete [] kstat.pass;

	delete [] readA;
	delete [] readA_data;

	return 0;
}


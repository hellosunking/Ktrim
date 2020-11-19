/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Feb, 2020
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
/*
	if( a=='A' )
		return b=='T';
	if( a=='C' )
		return b=='G';
	if( a=='G' )
		return b=='C';
	if( a=='T' )
		return b=='A';

	//TODO: consider how to deal with N, call it positive or negative???
	return false;
*/
	switch( a ) {
		case 'A': return b=='T';
		case 'C': return b=='G';
		case 'G': return b=='C';
		case 'T': return b=='A';
		default : return false;
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
					ktrim_stat * kstat, writeBuffer * writebuffer, const ktrim_param & kp ) {

	writebuffer->b1stored[tn] = 0;
	writebuffer->b2stored[tn] = 0;

//	vector<unsigned int> seed;
//	vector<unsigned int> :: iterator it, end_of_seed;
	register int *seed = new int[ MAX_SEED_NUM ];
	register int hit_seed;
	register int *it, *end_of_seed;

	register CPEREAD *wkr = workingReads + start;
	for( unsigned int iii=end-start; iii; --iii, ++wkr ) {
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
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				CPEREAD_resize( wkr, *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
			i = wkr->size - 2;
			const char *p = wkr->seq1;
			const char *q = wkr->seq2;
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] &&
					q[i]==kp.adapter_r2[0] && q[i+1]==kp.adapter_r2[1] ) {	// a possible hit
				// if it is a real adapter, then Read1 and Read2 should be complimentary
				// in real data, the heading 5 bp are of poor quality, therefoe we test the 6th, 7th
				if( is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) ) {
					++ kstat->tail_adapter[tn];
					if( i < kp.min_length ) {
						++ kstat->dropped[tn];
						continue;
					}
					CPEREAD_resize( wkr, i );
				}
			} else {	// tail 2 is not good, check tail 1
				++ i;
				if( p[i]==kp.adapter_r1[0] && q[i]==kp.adapter_r2[0] ) {
					if( is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) && 
						is_revcomp(p[6], q[i-7]) && is_revcomp(q[6], p[i-7]) ) {
						++ kstat->tail_adapter[tn];
						if( i < kp.min_length ) {
							++ kstat->dropped[tn];
							continue;
						}
						CPEREAD_resize( wkr, i );
					}
				}
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
											"%s%s\n+\n%s\n", wkr->id1, wkr->seq1, wkr->qual1 );
		writebuffer->b2stored[tn] += sprintf( writebuffer->buffer2[tn]+writebuffer->b2stored[tn],
											"%s%s\n+\n%s\n", wkr->id2, wkr->seq2, wkr->qual2 );
	}

	delete [] seed;
}

int process_multi_thread_PE_C( const ktrim_param &kp ) {
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	// in this version, two data containers are used and auto-swapped for working and loading data
	CPEREAD *readA = new CPEREAD[ READS_PER_BATCH ];
	CPEREAD *readB = new CPEREAD[ READS_PER_BATCH ];
	register char *readA_data = new char[ MEM_PE_READSET ];
	register char *readB_data = new char[ MEM_PE_READSET ];
	
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

	CPEREAD *workingReads, *loadingReads, *swapReads;

	ktrim_stat kstat;
	kstat.dropped	   = new unsigned int [ kp.thread ];
	kstat.real_adapter = new unsigned int [ kp.thread ];
	kstat.tail_adapter = new unsigned int [ kp.thread ];

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ kp.thread ];
	writebuffer.buffer2  = new char * [ kp.thread ];
	writebuffer.b1stored = new unsigned int	[ kp.thread ];
	writebuffer.b2stored = new unsigned int [ kp.thread ];

	for(unsigned int i=0; i!=kp.thread; ++i) {
		writebuffer.buffer1[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];
		writebuffer.buffer2[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];

		kstat.dropped[i] = 0;
		kstat.real_adapter[i] = 0;
		kstat.tail_adapter[i] = 0;
	}

// deal with multiple input files
	vector<string> R1s, R2s;
	extractFileNames( kp.FASTQ1, R1s );
	extractFileNames( kp.FASTQ2, R2s );

	if( R1s.size() != R2s.size() ) {
		fprintf( stderr, "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n" );
		return 110;
	}
	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	string fileName = kp.outpre;
	fileName += ".read1.fq";
	FILE *fout1 = fopen( fileName.c_str(), "wt" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	FILE *fout2 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL || fout2==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	register unsigned int line = 0;
	unsigned int threadCNT = kp.thread - 1;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		register const char * q = R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		}

		// initialization
		// get first batch of fastq reads

		unsigned int loaded;
		bool metEOF;
		if( file_is_gz ) {
			loaded = load_batch_data_PE_GZ( gfp1, gfp2, readA, READS_PER_BATCH );
			metEOF = gzeof( gfp1 );
		} else {
			loaded = load_batch_data_PE_C( fq1, fq2, readA, READS_PER_BATCH );
			metEOF = feof( fq1 );
		}
		if( loaded == 0 ) break;

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
//				clock_t start, end;
//				start = clock();

				unsigned int tn = omp_get_thread_num();
				// if EOF is met, then all threads are used for analysis
				// otherwise 1 thread will do data loading
				if( metEOF ) {
					NumWkThreads = kp.thread;
					unsigned int start = loaded * tn / kp.thread;
					unsigned int end   = loaded * (tn+1) / kp.thread;
					workingThread_PE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					nextBatch = false;
				} else {	// use 1 thread to load files, others for trimming
					NumWkThreads = threadCNT;

					if( tn == threadCNT ) {
						if( file_is_gz ) {
							threadLoaded = load_batch_data_PE_GZ( gfp1, gfp2, loadingReads, READS_PER_BATCH );
							metEOF = gzeof( gfp1 );
						} else {
							threadLoaded = load_batch_data_PE_C( fq1, fq2, loadingReads, READS_PER_BATCH );
							metEOF = feof( fq1 );
						}

						nextBatch = (threadLoaded!=0);
				//cerr << "Loading thread: " << threadLoaded << ", " << metEOF << ", " << nextBatch << '\n';
					} else {
						unsigned int start = loaded * tn / threadCNT;
						unsigned int end   = loaded * (tn+1) / threadCNT;
						workingThread_PE_C( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					}
				}
//				end = clock();
//				float duration = (end - start)*1000.0 / CLOCKS_PER_SEC;
//				fprintf( stderr, "Thread %d, runtime %.1f\n", tn, duration );
			} // parallel body
			// swap workingReads and loadingReads for next loop
			swapReads	 = loadingReads;
			loadingReads = workingReads;
			workingReads = swapReads;
			// write output and update fastq statistics
			for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
				fwrite( writebuffer.buffer1[ii], sizeof(char), writebuffer.b1stored[ii], fout1 );
			}
			for( unsigned int ii=0; ii!=NumWkThreads; ++ii ) {
				fwrite( writebuffer.buffer2[ii], sizeof(char), writebuffer.b2stored[ii], fout2 );
			}
			line += loaded;
			loaded = threadLoaded;
			//cerr << '\r' << line << " reads loaded";
		}

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fclose( fq1 );
			fclose( fq2 );
		}
	}

	fclose( fout1 );
	fclose( fout2 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=kp.thread; ++i ) {
		dropped_all += kstat.dropped[i];
		real_all += kstat.real_adapter[i];
		tail_all += kstat.tail_adapter[i];
	}
	fout << "Total\t"    << line		<< '\n'
		 << "Dropped\t"  << dropped_all << '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t"  << tail_all	<< '\n';
	fout.close();

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
	register char *read_data = new char[ MEM_PE_READSET ];
	
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
	kstat.dropped[0] = 0;
	kstat.real_adapter[0] = 0;
	kstat.tail_adapter[0] = 0;

	// buffer for storing the modified reads per thread
	writeBuffer writebuffer;
	writebuffer.buffer1  = new char * [ 1 ];
	writebuffer.buffer2  = new char * [ 1 ];
	writebuffer.b1stored = new unsigned int	[ 1 ];
	writebuffer.b2stored = new unsigned int	[ 1 ];
	writebuffer.buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];
	writebuffer.buffer2[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

// deal with multiple input files
	vector<string> R1s, R2s;
	extractFileNames( kp.FASTQ1, R1s );
	extractFileNames( kp.FASTQ2, R2s );

	if( R1s.size() != R2s.size() ) {
		fprintf( stderr, "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n" );
		return 110;
	}
	unsigned int totalFiles = R1s.size();
	//cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	FILE *fout1, *fout2;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wt" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	fout2 = fopen( fileName.c_str(), "wt" );
	if( fout1==NULL || fout2==NULL ) {
		fprintf( stderr, "\033[1;31mError: write file failed!\033[0m\n" );
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		FILE *fq1, *fq2;
		gzFile gfp1, gfp2;
		register unsigned int i = R1s[fileCnt].size() - 3;
		register const char * p = R1s[fileCnt].c_str();
		register const char * q = R2s[fileCnt].c_str();
		if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
			file_is_gz = true;
			gfp1 = gzopen( p, "r" );
			gfp2 = gzopen( q, "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		} else {
			fq1 = fopen( p, "rt" );
			fq2 = fopen( q, "rt" );
			if( fq1==NULL || fq2==NULL ) {
				fprintf( stderr, "\033[1;31mError: open fastq file failed!\033[0m\n" );
				fclose( fout1 );
				fclose( fout2 );
				return 104;
			}
		}

		register unsigned int last_seed;
		vector<unsigned int> seed;
		vector<unsigned int> :: iterator it;

		while( true ) {
			// get fastq reads
			unsigned int loaded;

			if( file_is_gz ) {
				loaded = load_batch_data_PE_GZ( gfp1, gfp2, read, READS_PER_BATCH_ST );
			} else {
				loaded = load_batch_data_PE_C( fq1, fq2, read, READS_PER_BATCH_ST );
			}

			if( loaded == 0 )
				break;
			
			workingThread_PE_C( 0, 0, loaded, read, &kstat, &writebuffer, kp );

			// write output and update fastq statistics
			fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], fout1 );
			fwrite( writebuffer.buffer2[0], sizeof(char), writebuffer.b2stored[0], fout2 );

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
	fclose( fout1 );
	fclose( fout2 );
	//cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) {
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		return 105;
	}

	fout << "Total\t"	 << line					<< '\n'
		 << "Dropped\t"  << kstat.dropped[0]		<< '\n'
		 << "Aadaptor\t" << kstat.real_adapter[0]	<< '\n'
		 << "TailHit\t"  << kstat.tail_adapter[0]	<< '\n';
	fout.close();

	delete [] read;
	delete [] read_data;

	return 0;
}


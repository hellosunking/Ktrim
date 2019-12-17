/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Dec, 2019
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <thread>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include "common.h"
using namespace std;

void workingThread_PE( unsigned int tn, unsigned int start, unsigned int end, CppPERead *workingReads,
					ktrim_stat * kstat, writeBuffer * writebuffer, const ktrim_param & kp ) {
	writebuffer->b1stored[tn] = 0;
	writebuffer->b2stored[tn] = 0;

	register int i, j;
	register unsigned int last_seed;
	vector<unsigned int> seed;
	vector<unsigned int> :: iterator it;
	const char *p, *q;

	for( unsigned int ii=start; ii!=end; ++ii ) {
		// quality control
		p = workingReads[ii].qual1.c_str();
		q = workingReads[ii].qual2.c_str();
		j = workingReads[ii].seq1.length();
		for( i=j-1; i; --i ) {
			if( p[i]>=kp.quality && q[i]>=kp.quality ) break;
		}
		++ i;
		if( i < kp.min_length ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != j ) {  // quality-trim occurs
			workingReads[ii].seq1.resize(  i );
			workingReads[ii].seq2.resize(  i );
			workingReads[ii].qual1.resize( i );
			workingReads[ii].qual2.resize( i );
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		seed.clear();
		for( i=0; (i=workingReads[ii].seq1.find(kp.adapter_index1, i)) != string::npos; ++i )
			seed.push_back( i );
		for( i=0; (i=workingReads[ii].seq2.find(kp.adapter_index2, i)) != string::npos; ++i )
			seed.push_back( i );

		sort( seed.begin(), seed.end() );

		last_seed = impossible_seed;	// a position which cannot be in seed
		for( it=seed.begin(); it!=seed.end(); ++it ) {
			if( *it != last_seed ) {
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
				if( check_mismatch_dynamic_PE( workingReads[ii].seq1, workingReads[ii].seq2, *it, kp ) )
					break;
					last_seed = *it;
				}
		}
		if( it != seed.end() ) {	// adapter found
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				workingReads[ii].seq1.resize(  *it );
				workingReads[ii].seq2.resize(  *it );
				workingReads[ii].qual1.resize( *it );
				workingReads[ii].qual2.resize( *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
			i = workingReads[ii].seq1.length() - 2;
			p = workingReads[ii].seq1.c_str();
			q = workingReads[ii].seq2.c_str();
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] &&
					q[i]==kp.adapter_r2[0] && q[i+1]==kp.adapter_r2[1] ) {
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				workingReads[ii].seq1.resize(  i );
				workingReads[ii].seq2.resize(  i );
				workingReads[ii].qual1.resize( i );
				workingReads[ii].qual2.resize( i );
			} else {	// tail 2 is not good, check tail 1
				++ i;
				if( p[i] == kp.adapter_r1[0] && q[i] == kp.adapter_r2[0] ) {
					++ kstat->tail_adapter[tn];
					if( i < kp.min_length ) {
						++ kstat->dropped[tn];
						continue;
					}
					workingReads[ii].seq1.resize(  i );
					workingReads[ii].seq2.resize(  i );
					workingReads[ii].qual1.resize( i );
					workingReads[ii].qual2.resize( i );
				}
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
											"%s/1\n%s\n+\n%s\n",
											workingReads[ii].id.c_str(),
											workingReads[ii].seq1.c_str(),
											workingReads[ii].qual1.c_str() );
		writebuffer->b2stored[tn] += sprintf( writebuffer->buffer2[tn]+writebuffer->b2stored[tn],
											"%s/2\n%s\n+\n%s\n",
											workingReads[ii].id.c_str(),
											workingReads[ii].seq2.c_str(),
											workingReads[ii].qual2.c_str() );
	}
}

int process_multi_thread_PE( const ktrim_param &kp ) {
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	// in this version, two data containers are used and auto-swapped for working and loading data
	CppPERead *readA = new CppPERead[ READS_PER_BATCH ];
	CppPERead *readB = new CppPERead[ READS_PER_BATCH ];

	CppPERead *workingReads, *loadingReads, *swapReads;

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
		cerr << "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n";
		return 110;
	}
	unsigned int totalFiles = R1s.size();
	cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	FILE *fout1, *fout2;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wb" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	fout2 = fopen( fileName.c_str(), "wb" );
	if( fout1==NULL || fout2==NULL ) {
		cout << "\033[1;31mError: write file failed!\033[0m\n";
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	ifstream fq1, fq2;
	register unsigned int line = 0;
	unsigned int threadCNT = kp.thread - 1;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		fq1.open( R1s[fileCnt].c_str() );
		fq2.open( R2s[fileCnt].c_str() );
		if( fq1.fail() || fq2.fail() ) {
			cout << "\033[1;31mError: open fastq file failed!\033[0m\n";
			fq1.close();
			fq2.close();
			fclose( fout1 );
			fclose( fout2 );
			return 104;
		}
		// initialization
		// get first batch of fastq reads
		unsigned int loaded = load_batch_data_PE( fq1, fq2, readA, READS_PER_BATCH );
		if( loaded == 0 ) break;
		bool metEOF = fq1.eof();

		loadingReads = readB;
		workingReads = readA;
		bool nextBatch = true;
		unsigned int threadLoaded;
		while( nextBatch ) {
			// start parallalization
			omp_set_num_threads( kp.thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				// if EOF is met, then all threads are used for analysis
				// otherwise 1 thread will do data loading
				if( metEOF ) {
					unsigned int start = loaded * tn / kp.thread;
					unsigned int end   = loaded * (tn+1) / kp.thread;
					workingThread_PE( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					nextBatch = false;
				} else {
					if( tn == threadCNT ) {
						threadLoaded = load_batch_data_PE( fq1, fq2, loadingReads, READS_PER_BATCH );
						metEOF = fq1.eof();
						nextBatch = (threadLoaded!=0);
						//cerr << "Loading thread: " << threadLoaded << ", " << metEOF << ", " << nextBatch << '\n';
					} else {
						unsigned int start = loaded * tn / threadCNT;
						unsigned int end   = loaded * (tn+1) / threadCNT;
						workingThread_PE( tn, start, end, workingReads, &kstat, &writebuffer, kp );
					}
				}
			} // parallel body
			// swap workingReads and loadingReads for next loop
			swapReads	= loadingReads;
			loadingReads = workingReads;
			workingReads = swapReads;
			// write output and update fastq statistics
			for( unsigned int ii=0; ii!=kp.thread; ++ii ) {
				fwrite( writebuffer.buffer1[ii], sizeof(char), writebuffer.b1stored[ii], fout1 );
			}
			for( unsigned int ii=0; ii!=kp.thread; ++ii ) {
				fwrite( writebuffer.buffer2[ii], sizeof(char), writebuffer.b2stored[ii], fout2 );
			}
			line += loaded;
			loaded = threadLoaded;
			cerr << '\r' << line << " reads loaded";
		}
		fq1.close();
		fq2.close();
	}

	fclose( fout1 );
	fclose( fout2 );
	cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		cerr << "\033[1;34mError: cannot write log file!\033[0m\n";
		return 105;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=kp.thread; ++i ) {
		dropped_all += kstat.dropped[i];
		real_all += kstat.real_adapter[i];
		tail_all += kstat.tail_adapter[i];
	}
	fout << "Total: "	<< line		<< '\n'
		 << "Dropped : " << dropped_all << '\n'
		 << "Aadaptor: " << real_all	<< '\n'
		 << "Tail Hit: " << tail_all	<< '\n';
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

	return 0;
}

int process_single_thread_PE( const ktrim_param &kp ) {
	// IO speed-up
	ios::sync_with_stdio( false );
//	cin.tie( NULL );

	CppPERead *read = new CppPERead[ READS_PER_BATCH_ST ];

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
		cerr << "\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n";
		return 110;
	}
	unsigned int totalFiles = R1s.size();
	cout << "\033[1;34mINFO: " << totalFiles << " paired fastq files will be loaded.\033[0m\n";

	FILE *fout1, *fout2;
	string fileName = kp.outpre;
	fileName += ".read1.fq";
	fout1 = fopen( fileName.c_str(), "wb" );
	fileName[ fileName.size()-4 ] = '2';	// read1 -> read2
	fout2 = fopen( fileName.c_str(), "wb" );
	if( fout1==NULL || fout2==NULL ) {
		cout << "\033[1;31mError: write file failed!\033[0m\n";
		fclose( fout1 );
		fclose( fout2 );
		return 103;
	}

	ifstream fq1, fq2;
	register unsigned int line = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		fq1.open( R1s[fileCnt].c_str() );
		fq2.open( R2s[fileCnt].c_str() );
		if( fq1.fail() || fq2.fail() ) {
			cout << "\033[1;31mError: open fastq file failed!\033[0m\n";
			fq1.close();
			fq2.close();
			fclose( fout1 );
			fclose( fout2 );
			return 104;
		}
		register int i, j;
		register unsigned int last_seed;
		vector<unsigned int> seed;
		vector<unsigned int> :: iterator it;
		const char *p, *q;

		while( true ) {
			// get fastq reads
			unsigned int loaded = load_batch_data_PE( fq1, fq2, read, READS_PER_BATCH_ST );
			if( loaded == 0 )
				break;
			
			workingThread_PE( 0, 0, loaded, read, &kstat, &writebuffer, kp );

			// write output and update fastq statistics
			fwrite( writebuffer.buffer1[0], sizeof(char), writebuffer.b1stored[0], fout1 );
			fwrite( writebuffer.buffer2[0], sizeof(char), writebuffer.b2stored[0], fout2 );

			line += loaded;
			cerr << '\r' << line << " reads loaded";

			if( fq1.eof() ) break;
		}
		fq1.close();
		fq2.close();
	}

	fclose( fout1 );
	fclose( fout2 );
	cerr << "\rDone: " << line << " lines processed.\n";

	// write trim.log
	fileName = kp.outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) {
		cerr << "\033[1;34mError: cannot write log file!\033[0m\n";
		return 105;
	}

	fout << "Total: "	 << line					<< '\n'
		 << "Dropped : " << kstat.dropped[0]		<< '\n'
		 << "Aadaptor: " << kstat.real_adapter[0]	<< '\n'
		 << "Tail Hit: " << kstat.tail_adapter[0]	<< '\n';
	fout.close();

	return 0;
}

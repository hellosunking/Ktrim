/*
 * common.h
 *
 * This header file records the constants used in Ktrim
 *
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date:   Mar, 2020
 * This program is part of the Ktrim package
 *
 **/

/*
 *  Sequencing model:
 *
 *                *read1 -->
 *  5' adapter - sequence sequence sequence - 3' adapter
 *                                <-- read2*
 *
 *  so read1 may contains 3' adapter, name it adapter_r1,
 *  read2 may contains reversed and ACGT-paired 5' adapter, name it adapter_r2
 *
 **/

#ifndef _KTRIM_COMMON_
#define _KTRIM_COMMON_

#include <string>
#include <vector>
using namespace std;

const char * VERSION = "1.4.0 (Oct 2021)";

// 1.4.0, now we use 2 threads for reading files in PE mode, ~33% faster
// 1.3.1 fixed the bug in dimers when working on SE data processing using single-thread 
// 1.2.2 fixed the bug when "-o" is NOT present but the program does not quit
// 1.2.1 fixed the bug in multi-file handling

// structure of a READ
const unsigned int MAX_READ_ID    = 128;
const unsigned int MAX_READ_CYCLE = 512;

// maxmimum insert size to call a dimer
const unsigned int DIMER_INSERT   = 1;

typedef struct {
	char *id;
	char *seq;
	char *qual;
	unsigned int size;
} CSEREAD;

typedef struct {
	char *id1;
	char *seq1;
	char *qual1;
	unsigned int size;

	char *id2;
	char *seq2;
	char *qual2;
	unsigned int size2;
} CPEREAD;

typedef struct {
	unsigned int *dropped;
	unsigned int *real_adapter;
	unsigned int *tail_adapter;
	unsigned int *dimer;
} ktrim_stat;

typedef struct {
	char ** buffer1;
	char ** buffer2;
	unsigned int *b1stored;
	unsigned int *b2stored;
} writeBuffer;

// built-in adapters
const unsigned int MIN_ADAPTER_SIZE = 8;
const unsigned int MAX_ADAPTER_SIZE = 64;
const unsigned int ADAPTER_BUFFER_SIZE = 128;
const unsigned int ADAPTER_INDEX_SIZE = 3;
const unsigned int OFFSET_INDEX3 = 3;
// illumina TruSeq kits adapters
//const char * illumina_adapter_r1 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";
//const char * illumina_adapter_r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
//const unsigned int illumina_adapter_len = 33;
const char * illumina_adapter_r1 = "AGATCGGAAGAGC";
const char * illumina_adapter_r2 = illumina_adapter_r1;
const unsigned int illumina_adapter_len = 13;
const char * illumina_index1 = "AGA";			// use first 3 as index
const char * illumina_index2 = illumina_index1;	// use first 3 as index
const char * illumina_index3 = "TCG";			// for single-end data
// Nextera kits and AmpliSeq for Illumina panels
const char * nextera_adapter_r1 = "CTGTCTCTTATACACATCT";
const char * nextera_adapter_r2 = nextera_adapter_r1;
const unsigned int nextera_adapter_len = 19;
const char * nextera_index1 = "CTG";			// use first 3 as index
const char * nextera_index2 = nextera_index1;	// use first 3 as index
const char * nextera_index3 = "TCT";			// for single-end data
// Nextera transposase adapters
const char * transposase_adapter_r1 = "TCGTCGGCAGCGTC";
const char * transposase_adapter_r2 = "GTCTCGTGGGCTCG";
const unsigned int transposase_adapter_len = 14;
const char * transposase_index1 = "TCG";		// use first 3 as index
const char * transposase_index2 = "GTC";		// use first 3 as index
const char * transposase_index3 = "TCG";		// for single-end data
// BGI adapters
//const char * bgi_adapter_r1 = "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA";
//const char * bgi_adapter_r2 = "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGC";
//const unsigned int bgi_adapter_len = 32;
const char * bgi_adapter_r1 = "AAGTCGGAGGCCAAGCGGTC";
const char * bgi_adapter_r2 = "AAGTCGGATCGTAGCCATGT";
const unsigned int bgi_adapter_len = 20;
const char * bgi_index1 = "AAG";		// use first 3 as index
const char * bgi_index2 = "AAG";		// use first 3 as index
const char * bgi_index3 = "TCG";		// for single-end data

// seed and error configurations
const unsigned int impossible_seed = 1000000;
const unsigned int MAX_READ_LENGTH = 1024;
const unsigned int MAX_SEED_NUM    = 128;

//configurations for parallelization, which is highly related to memory usage
//but seems to have very minor effect on running time
const int READS_PER_BATCH            = 1 << 15;	// process 128 K reads per batch (for parallelization)
const int BUFFER_SIZE_PER_BATCH_READ = 1 << 25;	// 128 MB buffer for each thread to store FASTQ
const int MEM_SE_READSET = READS_PER_BATCH * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE);
const int MEM_PE_READSET = READS_PER_BATCH * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE) * 2;

// enlarge the buffer for single-thread run
//const int READS_PER_BATCH_ST = READS_PER_BATCH << 1;	// process 256 K reads per batch (for parallelization)
//const int BUFFER_SIZE_PER_BATCH_READ_ST = BUFFER_SIZE_PER_BATCH_READ << 1;	// 256 MB buffer for each thread to store FASTQ
const int READS_PER_BATCH_ST            = 1 << 15;	// No. of reads per-batch for single-thread
const int BUFFER_SIZE_PER_BATCH_READ_ST = 1 << 25;	// buffer for single-thread to store FASTQ
const int MEM_SE_READSET_ST = READS_PER_BATCH_ST * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE);
const int MEM_PE_READSET_ST = READS_PER_BATCH_ST * (MAX_READ_ID+MAX_READ_CYCLE+MAX_READ_CYCLE) * 2;

const char FILE_SEPARATOR = ',';

// paramaters related
typedef struct {
	char *FASTQ1, *FASTQ2, *FASTQU;
	char *outpre;

	unsigned int thread;
	unsigned int min_length;
	unsigned int phred;
	unsigned int minqual;
	unsigned int quality;
	unsigned int window;

	const char *seqKit, *seqA, *seqB;
	const char *adapter_r1, *adapter_r2;
	unsigned int adapter_len;
	const char *adapter_index1, *adapter_index2, *adapter_index3;

	bool use_default_mismatch;
	float mismatch_rate;
} ktrim_param;

const char * param_list = "1:2:U:o:t:k:s:p:q:w:a:b:m:hv";

// definition of functions
void usage();
void init_param( ktrim_param &kp );
int  process_cmd_param( int argc, char * argv[], ktrim_param &kp );
void print_param( const ktrim_param &kp );
void extractFileNames( const char *str, vector<string> & Rs );

// C-style
unsigned int load_batch_data_PE_C( FILE * fq1, FILE * fq2, CPEREAD *loadingReads, unsigned int num );
bool check_mismatch_dynamic_PE_C( const CPEREAD *read, unsigned int pos, const ktrim_param &kp );
int process_single_thread_PE_C( const ktrim_param &kp );
int process_multi_thread_PE_C(  const ktrim_param &kp );

unsigned int load_batch_data_SE_C( FILE * fp, CSEREAD *loadingReads, unsigned int num );
bool check_mismatch_dynamic_SE_C( const char * p, unsigned int pos, const ktrim_param &kp );
int process_single_thread_SE_C( const ktrim_param &kp );
int process_multi_thread_SE_C(  const ktrim_param &kp );

#endif


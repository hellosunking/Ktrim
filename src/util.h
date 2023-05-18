/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Feb, 2020
 * This program is part of the Ktrim package
**/

#ifndef _KTRIM_UTIL_
#define _KTRIM_UTIL_

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <math.h>
#include <zlib.h>
#include "common.h"

using namespace std;

// extract file names
void extractFileNames( const char *str, vector<string> & Rs ) {
	string fileName = "";
	for(unsigned int i=0; str[i]!='\0'; ++i) {
		if( str[i] == FILE_SEPARATOR ) {
			Rs.push_back( fileName );
			fileName.clear();
		} else {
			fileName += str[i];
		}
	}
	if( fileName.size() )   // in case there is a FILE_SEPARATOR at the end
		Rs.push_back( fileName );
}

void loadFQFileNames( ktrim_param & kp ) {
	ifstream fin;
	fin.open( kp.filelist );
	if( fin.fail() ) {
		cerr << "Error: load file " << kp.filelist << " failed!\n";
		exit(10);
	}

	string line;
	stringstream ss;
	string fq1, fq2;

	// read the first valid line, determine SE/PE data
	int lineCnt = 0;
	bool pe_data;
	while( true ) {
		getline(fin, line);
		if( fin.eof() )break;

		if( line[0] != '#' ) {	// lines starts with "#" are ignored as comments
			ss.str( line );
			ss.clear();
			ss >> fq1;

			if( ss.rdbuf()->in_avail() ) {
				ss >> fq2;
				pe_data = true;

				kp.R1s.push_back( fq1 );
				kp.R2s.push_back( fq2 );
			} else {
				pe_data = false;
				kp.R1s.push_back( fq1 );
			}

			++ lineCnt;

			break;
		}
	}

	bool inconsistent = false;
	while( true ) {
		getline(fin, line);
		if( fin.eof() )break;

		if( line[0] == '#' ) {
			continue;
		}

		ss.str( line );
		ss.clear();
		ss >> fq1;
		
		if( pe_data ) {
			if( ss.rdbuf()->in_avail() ) {
				ss >> fq2;
				kp.R1s.push_back( fq1 );
				kp.R2s.push_back( fq2 );
			} else {
				inconsistent = true;
				break;
			}
		} else {
			if( ss.rdbuf()->in_avail() ) {
				inconsistent = true;
				break;
			}
			kp.R1s.push_back( fq1 );
		}

		++ lineCnt;
	}

	fin.close();
	if( inconsistent ) {
		cerr << "\033[1;31mERROR: inconsistent PE/SE in '" << kp.filelist << "'!\033[0m\n";
		exit(3);
	}
	if( pe_data ) {
		if( kp.R1s.size() != kp.R2s.size() ) {
			cerr << "\033[1;31mError: incorrect pairs in file list!\033[0m\n";
			exit(110);
		}
	}
//	cerr << "Done: " << lineCnt << " lines of " << (pe_data) ? 'P' : 'S' << "E data loaded.\n";
	kp.paired_end_data = pe_data;
}

//load 1 batch of data, using purely C-style
unsigned int load_batch_data_SE_C( FILE *fp, CSEREAD *loadingReads, unsigned int num ) {
	register unsigned int loaded = 0;
	register CSEREAD *p = loadingReads;
	register CSEREAD *q = p + num;
	while( p != q ) {
		if( fgets( p->id, MAX_READ_ID,  fp ) == NULL ) break;
		fgets( p->seq,  MAX_READ_CYCLE, fp );
		fgets( p->qual, MAX_READ_CYCLE, fp );	// this line is useless
		fgets( p->qual, MAX_READ_CYCLE, fp );

		// remove the tail '\n'
		register unsigned int i = strlen( p->seq ) - 1;
		p->size = i;
		p->seq[i]  = 0;
		p->qual[i] = 0;

		++ p;
	}
	return p-loadingReads;
}

// update in v1.2: support Gzipped input file
unsigned int load_batch_data_SE_GZ( gzFile gfp, CSEREAD *loadingReads, unsigned int num ) {
	register unsigned int loaded = 0;
//	string unk;
	register CSEREAD *p = loadingReads;
	register CSEREAD *q = p + num;
	register char c;
	register unsigned int i;

	while( p != q ) {
		if( gzgets( gfp, p->id, MAX_READ_ID ) == NULL ) break;
		gzgets( gfp, p->seq,  MAX_READ_CYCLE );
		gzgets( gfp, p->qual, MAX_READ_CYCLE );	// this line is useless
		gzgets( gfp, p->qual, MAX_READ_CYCLE );

		// remove the tailing '\n'
		i = strlen( p->seq ) - 1;
		p->size = i;
		p->seq[i]  = 0;
		p->qual[i] = 0;

		++ p;
	}

	return p-loadingReads;
}

unsigned int load_batch_data_PE_C( FILE *fq, CPEREAD *loadingReads, const unsigned int num, const bool isRead1 ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	// load read1
	if( isRead1 ) {
		while( p != q ) {
			if( fgets( p->id1, MAX_READ_ID,  fq ) == NULL ) break;
			fgets( p->seq1,  MAX_READ_CYCLE, fq );
			fgets( p->qual1, MAX_READ_CYCLE, fq );	// this line is useless
			fgets( p->qual1, MAX_READ_CYCLE, fq );
			p->size = strlen( p->seq1 );

			++ p;
		}
	} else {
		while( p != q ) {
			if( fgets( p->id2, MAX_READ_ID,  fq ) == NULL ) break;
			fgets( p->seq2,  MAX_READ_CYCLE, fq );
			fgets( p->qual2, MAX_READ_CYCLE, fq );	// this line is useless
			fgets( p->qual2, MAX_READ_CYCLE, fq );
			p->size2 = strlen( p->seq2 );

			++ p;
		}
	}
	return p - loadingReads;
}

// update in v1.2: support Gzipped input file
unsigned int load_batch_data_PE_GZ( gzFile gfp, CPEREAD *loadingReads, const unsigned int num, const bool isRead1 ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	if( isRead1 ) {
		while( p != q ) {
			if( gzgets( gfp, p->id1, MAX_READ_ID  ) == NULL ) break;
			gzgets( gfp, p->seq1,  MAX_READ_CYCLE );
			gzgets( gfp, p->qual1, MAX_READ_CYCLE );	// this line is useless
			gzgets( gfp, p->qual1, MAX_READ_CYCLE );
			p->size = strlen( p->seq1 );

			++ p;
		}
	} else {
		while( p != q ) {
			if( gzgets( gfp, p->id2, MAX_READ_ID  ) == NULL ) break;
			gzgets( gfp, p->seq2,  MAX_READ_CYCLE );
			gzgets( gfp, p->qual2, MAX_READ_CYCLE );	// this line is useless
			gzgets( gfp, p->qual2, MAX_READ_CYCLE );
			p->size2 = strlen( p->seq2 );

			++ p;
		}
	}
	return p - loadingReads;
}

unsigned int load_batch_data_PE_both_C( FILE *fq1, FILE *fq2, CPEREAD *loadingReads, unsigned int num ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	register CPEREAD *s = p;
	// load read1
	while( p != q ) {
		if( fgets( p->id1, MAX_READ_ID,  fq1 ) == NULL ) break;
		fgets( p->seq1,  MAX_READ_CYCLE, fq1 );
		fgets( p->qual1, MAX_READ_CYCLE, fq1 );	// this line is useless
		fgets( p->qual1, MAX_READ_CYCLE, fq1 );
		p->size = strlen( p->seq1 );

		++ p;
	}

	// load read2
	while( s != p ) {
		fgets( s->id2,   MAX_READ_ID,    fq2 );
		fgets( s->seq2,  MAX_READ_CYCLE, fq2 );
		fgets( s->qual2, MAX_READ_CYCLE, fq2 );	// this line is useless
		fgets( s->qual2, MAX_READ_CYCLE, fq2 );
		s->size2 = strlen( s->seq2 );

		++ s;
	}

	return p - loadingReads;
}

// update in v1.2: support Gzipped input file
unsigned int load_batch_data_PE_both_GZ( gzFile gfp1, gzFile gfp2, CPEREAD *loadingReads, unsigned int num ) {
	//register unsigned int loaded = 0;
	register CPEREAD *p = loadingReads;
	register CPEREAD *q = p + num;
	register CPEREAD *s = p;
	while( p != q ) {
		if( gzgets( gfp1, p->id1, MAX_READ_ID  ) == NULL ) break;
		gzgets( gfp1, p->seq1,  MAX_READ_CYCLE );
		gzgets( gfp1, p->qual1, MAX_READ_CYCLE );	// this line is useless
		gzgets( gfp1, p->qual1, MAX_READ_CYCLE );
		p->size = strlen( p->seq1 );

		++ p;
	}
	while( s != p ) {
		gzgets( gfp2, s->id2,   MAX_READ_ID    );
		gzgets( gfp2, s->seq2,  MAX_READ_CYCLE );
		gzgets( gfp2, s->qual2, MAX_READ_CYCLE );	// this line is useless
		gzgets( gfp2, s->qual2, MAX_READ_CYCLE );
		s->size2 = strlen( s->seq2 );

		++ s;
	}

	return p - loadingReads;
}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/8
*/
bool check_mismatch_dynamic_SE_C( CSEREAD *read, unsigned int pos, const ktrim_param & kp ) {
	register unsigned int mis=0;
	register unsigned int i, len;
	len = read->size - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;

	register unsigned int max_mismatch_dynamic;
	// update in v1.1.0: allows the users to set the proportion of mismatches
	if( kp.use_default_mismatch ) {
		max_mismatch_dynamic = len >> 3;
		if( (max_mismatch_dynamic<<3) != len )
			++ max_mismatch_dynamic;
	} else {
		max_mismatch_dynamic = ceil( len * kp.mismatch_rate );
	}

	register const char *p = read->seq;
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != kp.adapter_r1[i] ) {
			if( mis == max_mismatch_dynamic )
			  return false;

			++ mis;
		}
	}

	return true;
}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/4 for read1 + read2
*/
bool check_mismatch_dynamic_PE_C( const CPEREAD *read, unsigned int pos, const ktrim_param &kp ) {
	register unsigned int mis1=0, mis2=0;
	register unsigned int i, len;
	len = read->size - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;

	register unsigned int max_mismatch_dynamic;
	// update in v1.1.0: allows the users to set the proportion of mismatches
	// BUT it is highly discouraged
	if( kp.use_default_mismatch ) {
		// each read allows 1/8 mismatches of the total comparable length
		max_mismatch_dynamic = len >> 3;
		if( (max_mismatch_dynamic<<3) != len )
			++ max_mismatch_dynamic;
	} else {
		max_mismatch_dynamic = ceil( len * kp.mismatch_rate );
	}

	// check mismatch for each read
	register const char * p = read->seq1 + pos;
	register const char * q = kp.adapter_r1;
	for( i=0; i!=len; ++i, ++p, ++q ) {
		if( *p != *q ) {
			if( mis1 == max_mismatch_dynamic )
				return false;

			++ mis1;
		}
	}
	p = read->seq2 + pos;
	q = kp.adapter_r2;
	for( i=0; i!=len; ++i, ++p, ++q ) {
		if( *p != *q ) {
			if( mis2 == max_mismatch_dynamic )
				return false;

			++ mis2;
		}
	}

	return true;
}

// update in v1.2: support window check
int get_quality_trim_cycle_se( const char *p, const int total_size, const ktrim_param &kp ) {
	register int i, j, k;
	register int stop = kp.min_length-1;
	for( i=total_size-1; i>=stop; ) {
		if( p[i] >= kp.quality ) {
			k = i - kp.window;
			for( j=i-1; j!=k; --j ) {
				if( j<0 || p[j]<kp.quality ) {
					break;
				}
			}
			if( j == k ) { // find the quality trimming position
				break;
			} else {	// there is a low-quality base in the middle
				i = j - 1;
			}
		} else {
			-- i;
		}
	}

	if( i >= stop )
		return i + 1;
	else
		return 0;
}

int get_quality_trim_cycle_pe( const CPEREAD *read, const ktrim_param &kp ) {
	register int i, j, k;
	register int stop = kp.min_length - 1;
	const char *p = read->qual1;
	const char *q = read->qual2;
	for( i=read->size-1; i>=stop; ) {
		if( p[i]>=kp.quality && q[i]>=kp.quality ) {
			k = i - kp.window;
			for( j=i-1; j!=k; --j ) {
				if( j<0 || p[j]<kp.quality || q[j]<kp.quality ) {
					break;
				}
			}
			if( j == k ) { // find the quality trimming position
				break;
			} else {
				i = j - 1;
			}
		} else {
			-- i;
		}
	}
	
	if( i >= stop )
		return i + 1;
	else
		return 0;
}

#endif


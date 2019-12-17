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

// extract file names
void extractFileNames( const char *str, vector<string> & Rs ) {
	string fileName="";
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

//load 1 batch of data
unsigned int load_batch_data_SE( ifstream & fq1, CppSERead *loadingReads, unsigned int num ) {
	register unsigned int loaded = 0;
	string unk;
	register CppSERead *p = loadingReads;
	register CppSERead *q = p + num;
	while( p != q ) {
		getline( fq1, p->id );
		if( fq1.eof() )break;
		getline( fq1, p->seq );
		getline( fq1, unk );
		getline( fq1, p->qual );

		++ p;
	}
	return p-loadingReads;
}

unsigned int load_batch_data_PE( ifstream &fq1, ifstream &fq2, CppPERead *loadingReads, unsigned int num ) {
	//register unsigned int loaded = 0;
	string unk;
	register CppPERead *p = loadingReads;
	register CppPERead *q = p + num;
	register CppPERead *s = p;
	while( p != q ) {
		getline( fq1, p->id );
		if( fq1.eof() )break;
		getline( fq1, p->seq1 );
		getline( fq1, unk );
		getline( fq1, p->qual1 );

		++ p;
	}
	while( s != p ) {
		getline( fq2, unk );
		getline( fq2, s->seq2 );
		getline( fq2, unk );
		getline( fq2, s->qual2 );

		++ s;
	}
	return p-loadingReads;
}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/8
*/
bool check_mismatch_dynamic_SE( const string & s, unsigned int pos, const ktrim_param & kp ) {
	register unsigned int mis=0;
	register unsigned int i, len;
	len = s.length() - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;

	register unsigned int max_mismatch_dynamic = len >> 3;
	if( (max_mismatch_dynamic<<3) != len )
		++ max_mismatch_dynamic;

	const char * p = s.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != kp.adapter_r1[i] ) {
			++ mis;
			if( mis > max_mismatch_dynamic )
			  return false;
		}
	}

	return true;
}

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * here the maximum mismatch allowed is LEN/4 for read1 + read2
*/
bool check_mismatch_dynamic_PE( const string & s1, string & s2, unsigned int pos, const ktrim_param &kp ) {
	register unsigned int mis1=0, mis2=0;
	register unsigned int i, len;
	len = s1.length() - pos;
	if( len > kp.adapter_len )
		len = kp.adapter_len;

	register unsigned int max_mismatch_dynamic = len >> 2;
	// here use 1/8 of the total length, use roof(len/4)
	if( (max_mismatch_dynamic<<2) != len )
		++ max_mismatch_dynamic;

	// check mismatch for each read
	register const char * p = s1.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != kp.adapter_r1[i] ) {
			++ mis1;
			if( mis1 > max_mismatch_dynamic )
				return false;
		}
	}
	p = s2.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != kp.adapter_r2[i] ) {
			++ mis2;
			if( mis2 > max_mismatch_dynamic )
				return false;
		}
	}
	// check total mismatch
	//max_mismatch_dynamic = (len+1) >> 1;
	if( mis1 + mis2 > max_mismatch_dynamic )
		return false;

	return true;
}


/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Feb, 2020
 * Main program of Ktrim
**/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>
#include "common.h"
#include "util.h"
#include "param_handler.h"
#include "pe_handler.h"
#include "se_handler.h"

using namespace std;

int main( int argc, char *argv[] ) {
	// process the command line parameters
	static ktrim_param kp;
	init_param( kp );
	int retValue = process_cmd_param( argc, argv, kp );
	if( retValue == 100 )	// help or version
		return 0;
	else if( retValue != 0 )
		return retValue;

	if( kp.FASTQ2 == NULL ) {  // single-end data
		if( kp.thread == 1 )
			retValue = process_single_thread_SE_C( kp );
		else
			retValue = process_multi_thread_SE_C( kp );
	} else {
		if( kp.thread == 1 )
			retValue = process_single_thread_PE_C( kp );
		else if( kp.thread == 2 )
			retValue = process_two_thread_PE_C( kp );
		else
			retValue = process_multi_thread_PE_C( kp );
	}

	return retValue;
}


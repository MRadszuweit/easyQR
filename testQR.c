#include "easyQR.h"

/**
 * \brief Test program to test the QR algorithm
 * \param arg #1: the number of nodes to use (default value is n=200)
 * \param arg #2: the number of threads to start
 * 
 * This little test program calls the test routine "testQR" from EasyQR.c.
 * It computes the eigensystem of some matrix and plots it (see testQR)
 */
 
int main(int argc, char* argv[]){
	const int default_dim = 200;
	const int default_proc = 1;
	
	if (argc>2) testQR(atoi(argv[1]),atoi(argv[2])); 
	else if (argc>1) testQR(atoi(argv[1]),default_proc); 
	else testQR(default_dim,default_proc);	
	return 0;
}


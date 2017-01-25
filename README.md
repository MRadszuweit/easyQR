# easyQR

This is a small library for finding all eigenvalues and eigenvectors of a symmetric (nxn)-matrix. 
For this, the QR algorithm is chosen. It is suited for small (n<1000) but dense matrices. 
There are no dependencies. Just copy the *.h and .c file in your project folder. Additionaly there is a small 
test program called "testQR" and a makefile. The library can optionally make use of openmp. Compile options
via the parameter OPT can be passed:

OPT=debug   for debugging
OPT=serial    no openmp
OPT=O3    with O3 and openmp

example: make -B OPT=O3

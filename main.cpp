#include <iostream>                                                                            
#include <stdlib.h>                                                                            
#include <stdio.h>                                                                             
#include "icoord.h"                                                                            
                                                                                               
using namespace std;                                                                           
                                                                                               
int main(int argc, char* argv[]){                                                              
  string xyzfile, gradfile;                                                                             
  if (argc < 3){                                                                               
    cout << "Must include xyzfile" <<endl;                                                     
    return -1;                                                                                 
  }                                                                                            
  xyzfile=argv[1];                                                                            
 	gradfile=argv[2];                                                                                              
                                                                                                                   
  ICoord test;                                                                                 
	test.isOpt=1;
	test.farBond=1.;
  test.init(xyzfile,gradfile);                                                                         
/*
 	test.bmat_alloc();                                                                                              
	test.bmatp_create();
	test.bmatp_to_U();
	test.bmat_create();
*/	
  return 0;                                                                                    
}                             

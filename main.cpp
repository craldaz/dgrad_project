#include <iostream>                                                                            
#include <stdlib.h>                                                                            
#include <stdio.h>                                                                             
#include "icoord.h"                                                                            
                                                                                               
using namespace std;                                                                           
                                                                                               
int main(int argc, char* argv[]){                                                              
  string xyzfile1;                                                                             
  if (argc < 2){                                                                               
    cout << "Must include xyzfile" <<endl;                                                     
    return -1;                                                                                 
  }                                                                                            
  xyzfile1=argv[1];                                                                            
                                                                                               
                                                                                                                   
  ICoord test;                                                                                 
	test.isOpt=1;
	test.farBond=1.;
  test.init(xyzfile1);                                                                         
                                                                                               
  return 0;                                                                                    
}                             

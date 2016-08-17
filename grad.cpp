#include "grad.h"

using namespace std;

void Gradient::grad_read(string gradfile,double* grad)
{ 
   
  cout <<" Reading gradient" << endl;
  
  ifstream infile;
  infile.open(gradfile.c_str());
  if (!infile){
    cout << "!!!!Error opening grad file!!!!" << endl;
    exit(-1);
  } 
  
  cout <<"  -reading file..." << endl;
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  cout <<"  natoms: " << natoms << endl;
  
  success=getline(infile, line);
  
  anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];
    
  cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    anumbers[i]=PTable::atom_number(anames[i]);
    amasses[i]=PTable::atom_mass(anumbers[i]);
  }
  
  infile.close();
  
  cout <<"  -Reading coordinates...";
  cout << "Opening the grad file" << endl;
  infile.open(gradfile.c_str());
  fflush(stdout);
  cout << "  gradfile opened" << endl;
  fflush(stdout);
  
  
    success=getline(infile, line);
    success=getline(infile, line);
    for (int j=0;j<natoms;j++){
      success=getline(infile, line);
      int length=StringTools::cleanstring(line);
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      grad[3*j+0]=atof(tok_line[1].c_str());
      grad[3*j+1]=atof(tok_line[2].c_str());
      grad[3*j+2]=atof(tok_line[3].c_str());
    }
  
   
  infile.close();
  cout << " Finished reading information from structure file" << endl;
}



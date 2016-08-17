#ifndef ICOORD_H
#define ICOORD_H

#include "stringtools.h"
#include "pTable.h"
#include "grad.h"

class ICoord {

  private:
    double* amasses;              //array of atomic masses
    double* amasses3;

    int nfrags;
    int* frags;
    int** imptor;
    int max_bonds;
    int max_angles;
    int max_torsions;
    int max_imptor;

    int max_nonbond;
    int n_nonbond;
    int** nonbond;
    double* nonbondd;
  	void structure_read(string xyzfile);
  	void alloc_mem();
  	void make_bonds();
  	void coord_num();
		void make_imptor();
  	void make_imptor_nobonds(); 

  	int make_nonbond();
  	void make_angles();
  	void make_torsions();

  	void make_frags();
  	void bond_frags();
  	void hbond_frags();
  	void linear_ties();

	public:
  	int natoms;
  	double* coords;
  	double* coords0;
  	string* anames;               //array of atomic symbols 
  	int* anumbers;                //array of atomic indices 
  	int* coordn;                  //coordination number
  	int nimptor;

  	Gradient grad1;
  	double farBond;
		
    double* bondd;
    double* anglev;
    double* torv;
    double* torv0;
    double* torfix;
    double* imptorv;

  	int** bonds;
  	int nbonds;
  	int** angles;
  	int nangles;
  	int** torsions;
  	int ntor;
		int isOpt;
		void print_xyz();
  	void print_ic();
		int init(string xyzfile,string gradfile);
		int ic_create();	
  	void update_ic();
  	void update_bonds();
  	void update_angles();
  	void update_torsion();
  	void update_imptor();
  	void update_nonbond();
		double getR(int i);	
  	int nicd;
  	int ixflag;
    int nicd0; //before constraint applied
  	double distance(int i, int j);
  	double angle_val(int i, int j, int k);
  	double torsion_val(int i, int j, int k, int l); // for imptor and torsion
  	int bond_exists(int b1, int b2);
  	int bond_num(int b1, int b2);
  	int angle_num(int b1, int b2, int b3);
  	int tor_num(int b1, int b2, int b3, int b4);
  	int isTM(int anum);
  	// Gradient terms
  	int bmat_alloc();
		int bmatp_create();
  	void bmatp_dqbdx(int a1, int a2, double* dqbdx);
  	void bmatp_dqadx(int a1, int a2, int a3, double* dqadx);
  	void bmatp_dqtdx(int a1, int a2, int a3, int a4, double* dqtdx);
  	int bmatp_to_U();
		int bmat_create();
		int grad_to_q();
		void print_gradq();
		void project_grad();
		int ic_to_xyz();
  	double* Ut;
  	double* Ut0;
  	double* bmat;
  	double* bmatti;
  	double* bmatp; // in primitives
  	double* grad;
		double* dgrad;
		 double* pgrad;
  	double* gradq;
		double* dgradq;
  	double* pgradq;
  	double* gradqprim;
  	double* pgradqprim;
		double gradrms;
		double pgradrms;
  	double* dq0;
  	double* dqprim;
  	double* q;
  	double MAXAD;
  	double DMAX;
  	double DMIN0;
};


#endif

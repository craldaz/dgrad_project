#include "icoord.h"


void ICoord::alloc_mem()
{
 grad = new double[3*natoms];
 nbonds = 0;
 max_bonds=natoms*natoms+8;
 bonds = new int*[max_bonds];
 for (int i=0;i<max_bonds;i++)
   bonds[i]=new int[2];
 bondd = new double[max_bonds];

 nangles = 0;
 max_angles=natoms*24;
 angles = new int*[max_angles];
 for (int i=0;i<max_angles;i++)
   angles[i]=new int[3];
 anglev = new double[max_angles];

 ntor = 0;
 max_torsions=natoms*60;
 torsions = new int*[max_torsions];
 for (int i=0;i<max_torsions;i++)
   torsions[i] = new int[4];
 torv = new double[max_torsions];

 nimptor = 0;
 max_imptor=natoms*10;
 imptor = new int*[max_imptor];
 for (int i=0;i<max_imptor;i++)
   imptor[i] = new int[4];
 imptorv = new double[max_imptor];


 n_nonbond = 0;
 max_nonbond = natoms*natoms;
 nonbond = new int*[max_nonbond];
 for (int i=0;i<max_nonbond;i++)
   nonbond[i]=new int[2];
 nonbondd = new double[max_nonbond];

 coordn = new int[natoms];
 for (int i=0;i<natoms;i++)
    coordn[i]=0;


 return;
}

#include "icoord.h"
#include "utils.h"
#include <mkl.h>

using namespace std; 

#define THRESH 1E-3

int ICoord::bmat_alloc()
{
	int size_ic	= nbonds+nangles+ntor +150; //buffer of 150 for new primitives
	int size_xyz = 3*natoms; 
  printf(" in bmat_alloc, size_ic: %i size_xyz: %i \n",size_ic-150,size_xyz);
  printf(" max_bonds: %i max_angles: %i max_torsions: %i \n",max_bonds,max_angles,max_torsions);
  bmat = new double[size_ic*size_xyz+100];
  bmatp = new double[size_ic*size_xyz+1000];
  bmatti = new double[size_ic*size_xyz+100];
  torv0 = new double[max_torsions+100];
  for (int i=0;i<max_torsions;i++) torv0[i] = 0;
  torfix = new double[max_torsions+100];
  Ut = new double[size_ic*size_ic+100];
  Ut0 = new double[size_ic*size_ic+100];
  q = new double[size_ic+100];
  dq0 = new double[size_ic+100];
  dqprim = new double[size_ic+100];
  pgradq = new double[size_ic+100];
  gradq = new double[size_ic+100];
	dgradq = new double[size_ic+100]; 
  for (int i=0;i<size_ic+100;i++) pgradq[i]=0.;
  for (int i=0;i<size_ic+100;i++) gradq[i]=0.;
	for (int i=0;i<size_ic+100;i++) dgradq[i]=0.;
  pgradqprim = new double[size_ic+100];
  gradqprim = new double[size_ic+100];
  for (int i=0;i<size_ic+100;i++) pgradqprim[i]=0.;
  for (int i=0;i<size_ic+100;i++) gradqprim[i]=0.;
 
  MAXAD = 0.075; //max along one coordinate (was using 0.1)
  DMAX = 0.1; //max of step magnitude (was using 0.125)
  DMIN0 = DMAX/5.; //was 5.
#if USE_PRIMA
  DMAX = 0.025;
#endif


  return 0;
}



int ICoord::bmatp_create() 
{

  printf(" in bmatp_create \n");

  int len = nbonds+nangles+ntor;
  int N3 = 3*natoms;
  int max_size_ic = len;
  int size_xyz = N3;

  for (int i=0;i<max_size_ic*N3;i++)
    bmatp[i] = 0.;

  double* dqbdx = new double[6];
  for (int i=0;i<nbonds;i++)
  {
    for (int j=0;j<6;j++) dqbdx[j] = 0.;
    int a1=bonds[i][0];
    int a2=bonds[i][1];
    bmatp_dqbdx(a1,a2,dqbdx);
    double FACTOR = 1;
    bmatp[i*N3+3*a1+0] = dqbdx[0]*FACTOR;
    bmatp[i*N3+3*a1+1] = dqbdx[1]*FACTOR;
    bmatp[i*N3+3*a1+2] = dqbdx[2]*FACTOR;
    bmatp[i*N3+3*a2+0] = dqbdx[3]*FACTOR;
    bmatp[i*N3+3*a2+1] = dqbdx[4]*FACTOR;
    bmatp[i*N3+3*a2+2] = dqbdx[5]*FACTOR;
  }

  double* dqadx = new double[9];
  for (int i=nbonds;i<nbonds+nangles;i++)
  {
    for (int j=0;j<9;j++) dqadx[j] = 0.;
    int a1=angles[i-nbonds][0];
    int a2=angles[i-nbonds][1];
    int a3=angles[i-nbonds][2];
    bmatp_dqadx(a1,a2,a3,dqadx);
    bmatp[i*N3+3*a1+0] = dqadx[0];
    bmatp[i*N3+3*a1+1] = dqadx[1];
    bmatp[i*N3+3*a1+2] = dqadx[2];
    bmatp[i*N3+3*a2+0] = dqadx[3];
    bmatp[i*N3+3*a2+1] = dqadx[4];
    bmatp[i*N3+3*a2+2] = dqadx[5];
    bmatp[i*N3+3*a3+0] = dqadx[6];
    bmatp[i*N3+3*a3+1] = dqadx[7];
    bmatp[i*N3+3*a3+2] = dqadx[8];
  }
  double* dqtdx = new double[12];
  for (int i=nbonds+nangles;i<nbonds+nangles+ntor;i++)
  {
    for (int j=0;j<12;j++) dqtdx[j] = 0.;
    int a1=torsions[i-nbonds-nangles][0];
    int a2=torsions[i-nbonds-nangles][1];
    int a3=torsions[i-nbonds-nangles][2];
    int a4=torsions[i-nbonds-nangles][3];
    bmatp_dqtdx(a1,a2,a3,a4,dqtdx);
    bmatp[i*N3+3*a1+0] = dqtdx[0]; //*1.8897
    bmatp[i*N3+3*a1+1] = dqtdx[1];
    bmatp[i*N3+3*a1+2] = dqtdx[2];
    bmatp[i*N3+3*a2+0] = dqtdx[3];
    bmatp[i*N3+3*a2+1] = dqtdx[4];
    bmatp[i*N3+3*a2+2] = dqtdx[5];
    bmatp[i*N3+3*a3+0] = dqtdx[6];
    bmatp[i*N3+3*a3+1] = dqtdx[7];
    bmatp[i*N3+3*a3+2] = dqtdx[8];
    bmatp[i*N3+3*a4+0] = dqtdx[9];
    bmatp[i*N3+3*a4+1] = dqtdx[10];
    bmatp[i*N3+3*a4+2] = dqtdx[11];
  }

  printf(" \n after creating bmatp \n");
#if 1
  printf(" printing bond contributions \n");
  for (int i=0;i<nbonds;i++)
  {
    for (int j=0;j<natoms;j++)
    {
      for (int k=0;k<3;k++)
        printf(" %1.3f",bmatp[i*size_xyz+3*j+k]);
    }
    printf(" \n");
  }
#endif
#if 1
  printf(" printing angle contributions \n");
  for (int i=nbonds;i<nbonds+nangles;i++)
  {
    for (int j=0;j<natoms;j++)
    {
      for (int k=0;k<3;k++)
        printf(" %1.3f",bmatp[i*size_xyz+3*j+k]);
    }
    printf(" \n");
  }
#endif

#if 1
  int nztor = 0;
  double* x = new double[3];
  printf(" printing torsion contributions \n");
  for (int i=nbonds+nangles;i<nbonds+nangles+ntor;i++)
  {
    x[0] = x[1] = x[2] = 0.;
    for (int j=0;j<natoms;j++)
    {
//      for (int k=0;k<3;k++)
//        printf(" %7.3f",bmatp[i*size_xyz+3*j+k]);
      for (int k=0;k<3;k++) 
        x[k] += bmatp[i*N3+3*j+k]*bmatp[i*N3+3*j+k];
    }
    printf("   mag: %8.4f %8.4f %8.4f ",x[0],x[1],x[2]);
    printf(" \n");
    if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]>0.001)
      nztor++;
  }
  printf("  non-zero torsions: %2i \n",nztor);
  delete [] x;
#endif


  delete [] dqbdx;
  delete [] dqadx;
  delete [] dqtdx;

  return 0;
}
void ICoord::bmatp_dqbdx(int i, int j, double* dqbdx) 
{

  double* u = new double[3];  
  u[0] = coords[3*i+0]-coords[3*j+0];
  u[1] = coords[3*i+1]-coords[3*j+1];
  u[2] = coords[3*i+2]-coords[3*j+2];

  double norm = distance(i,j);
  //double norm = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

  u[0] = u[0]/norm;
  u[1] = u[1]/norm;
  u[2] = u[2]/norm;

  dqbdx[0] = u[0];
  dqbdx[1] = u[1];
  dqbdx[2] = u[2];
  dqbdx[3] = -u[0];
  dqbdx[4] = -u[1];
  dqbdx[5] = -u[2];

  delete [] u;

  return;
}
void ICoord::bmatp_dqadx(int i, int j, int k, double* dqadx) 
{

  double angle = angle_val(i,j,k) *3.14159/180; // in radians
  //printf(" angle_val: %1.2f \n",angle);
#if 0
  if (angle>3.0)
  {
   // printf(" near-linear angle, using finite difference \n");
    double fstep=0.0001;
    int a1=i;
    int a2=j;
    int a3=k;

    double b0,b1,b2,b3;

    b0 = angle_val(a1,a2,a3)*3.14159/180;
    for (int j=0;j<3;j++)
    {
      coords[3*a1+j]+=fstep;
      b1 = angle_val(a1,a2,a3)*3.14159/180;
      coords[3*a1+j]-=fstep;
      coords[3*a2+j]+=fstep;
      b2 = angle_val(a1,a2,a3)*3.14159/180;
      coords[3*a2+j]-=fstep;
      coords[3*a3+j]+=fstep;
      b3 = angle_val(a1,a2,a3)*3.14159/180;
      coords[3*a3+j]-=fstep;

      dqadx[3*0+j] = (b1-b0)/fstep;
      dqadx[3*1+j] = (b2-b0)/fstep;
      dqadx[3*2+j] = (b3-b0)/fstep;
    }
    return;
  }
#endif

  double* u = new double[3];
  double* v = new double[3];
  u[0] = coords[3*i+0]-coords[3*j+0];
  u[1] = coords[3*i+1]-coords[3*j+1];
  u[2] = coords[3*i+2]-coords[3*j+2];
  v[0] = coords[3*k+0]-coords[3*j+0];
  v[1] = coords[3*k+1]-coords[3*j+1];
  v[2] = coords[3*k+2]-coords[3*j+2];

  double n1 = distance(i,j);
  double n2 = distance(j,k);

  u[0] = u[0]/n1;  
  u[1] = u[1]/n1;  
  u[2] = u[2]/n1;  
  v[0] = v[0]/n2;  
  v[1] = v[1]/n2;  
  v[2] = v[2]/n2;  

  double* w = new double[3];
  w[0]=w[1]=w[2] = 0.;

  cross(w,u,v);
  double mag = (w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  if (mag<THRESH)
  { 
//    printf(" Linear angle detected, w: %1.6f %1.6f %1.6f \n",w[0],w[1],w[2]);
    double* vn = new double[3];
    vn[0]=0.; vn[1]=0.; vn[2]=1.;
    cross(w,u,vn);
    mag = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if (mag<THRESH)
    {
//      printf(" Linear angle(b) detected, w: %1.6f %1.6f %1.6f \n",w[0],w[1],w[2]);
      vn[0]=0.; vn[1]=1.; vn[2]=0.;
      cross(w,u,vn);
    }
    delete [] vn;
  }
//  if (angle>3.0)
//    printf(" w: %1.3f %1.3f %1.3f \n",w[0],w[1],w[2]);

  double n3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  //printf(" n3: %1.3f \n",n3);
  w[0] = w[0]/n3;
  w[1] = w[1]/n3;
  w[2] = w[2]/n3;
 
  double* uw = new double[3]; 
  double* wv = new double[3]; 
  cross(uw,u,w);
  cross(wv,w,v);

  dqadx[0] = uw[0]/n1;
  dqadx[1] = uw[1]/n1;
  dqadx[2] = uw[2]/n1;
  dqadx[3] = -uw[0]/n1 + -wv[0]/n2;
  dqadx[4] = -uw[1]/n1 + -wv[1]/n2;
  dqadx[5] = -uw[2]/n1 + -wv[2]/n2;
  dqadx[6] = wv[0]/n2;
  dqadx[7] = wv[1]/n2;
  dqadx[8] = wv[2]/n2;

#if 0
  if (angle>3.0)
  {
    printf(" uw: %1.3f %1.3f %1.3f, n1: %1.3f \n",uw[0],uw[1],uw[2],n1);
    printf(" wv: %1.3f %1.3f %1.3f, n2: %1.3f \n",wv[0],wv[1],wv[2],n2);
    printf(" dqadx: %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f \n",dqadx[0],dqadx[1],dqadx[2],dqadx[3],dqadx[4],dqadx[5],dqadx[6],dqadx[7],dqadx[8]);
  }
#endif

  delete [] u;
  delete [] v;
  delete [] w;
  delete [] uw;
  delete [] wv;

  return;
}
void ICoord::bmatp_dqtdx(int i, int j, int k, int l, double* dqtdx) 
{

  //printf(" \n beginning torsion bmat: %i %i %i %i \n",i,j,k,l);

  double angle1 = angle_val(i,j,k) *3.14159/180; // in radians
  double angle2 = angle_val(j,k,l) *3.14159/180; // in radians
  if (angle1>3.0 || angle2>3.0)
  {
    //printf(" near-linear angle, skipping bmat element \n");
    return;
  }
  //printf(" angle1,2: %1.2f %1.2f (%1.1f %1.1f) \n",angle1,angle2,angle1*180/3.14159,angle2*180/3.14159);

//u is between first and second atoms
//w is between third and second atoms
//v is between fourth and third atoms

  double* u = new double[3];
  double* w = new double[3];
  double* v = new double[3];
  u[0] = coords[3*i+0]-coords[3*j+0];
  u[1] = coords[3*i+1]-coords[3*j+1];
  u[2] = coords[3*i+2]-coords[3*j+2];
  w[0] = coords[3*k+0]-coords[3*j+0];
  w[1] = coords[3*k+1]-coords[3*j+1];
  w[2] = coords[3*k+2]-coords[3*j+2];
  v[0] = coords[3*l+0]-coords[3*k+0];
  v[1] = coords[3*l+1]-coords[3*k+1];
  v[2] = coords[3*l+2]-coords[3*k+2];

  double n1 = distance(i,j);
  double n2 = distance(j,k);
  double n3 = distance(k,l);
 
  //printf(" n1,n2,n3: %1.3f %1.3f %1.3f \n",n1,n2,n3);

  u[0] = u[0]/n1;  
  u[1] = u[1]/n1;  
  u[2] = u[2]/n1;  
  w[0] = w[0]/n2;  
  w[1] = w[1]/n2;  
  w[2] = w[2]/n2;  
  v[0] = v[0]/n3;  
  v[1] = v[1]/n3;  
  v[2] = v[2]/n3;  

  double* uw = new double[3]; 
  double* vw = new double[3]; 
  cross(uw,u,w);
  cross(vw,v,w);

  //double n4 = sqrt(uw[0]*uw[0]+uw[1]*uw[1]+uw[2]*uw[2]);
  //double n5 = sqrt(vw[0]*vw[0]+vw[1]*vw[1]+vw[2]*vw[2]);
  //do not normalize uw and vw

  double cosphiu = u[0]*w[0] + u[1]*w[1] + u[2]*w[2];
  double cosphiv = -v[0]*w[0] - v[1]*w[1] - v[2]*w[2];
//  double sin2phiu = sqrt(1-cosphiu*cosphiu);
//  double sin2phiv = sqrt(1-cosphiv*cosphiv);
  double sin2phiu = 1-cosphiu*cosphiu;
  double sin2phiv = 1-cosphiv*cosphiv;

  //printf(" cos's: %1.4f %1.4f vs %1.4f %1.4f \n",cosphiu,cosphiv,cos(angle1),cos(angle2));
  //printf(" sin2's: %1.4f %1.4f vs %1.4f %1.4f \n",sin2phiu,sin2phiv,sin(angle1)*sin(angle1),sin(angle2)*sin(angle2));

  if (sin2phiu<THRESH || sin2phiv<THRESH)
  { 
  //  printf(" sin2phi too small, not creating element \n");
    return;
  }

  //printf(" angle1,2: %1.2f %1.2f \n",angle1,angle2);
  //printf(" n1,n2,n3: %1.2f %1.2f %1.2f sin2phiu,v: %1.2f %1.2f cosphiu,v: %1.2f %1.2f \n",n1,n2,n3,sin2phiu,sin2phiv,cosphiu,cosphiv);


//CPMZ possible error in uw calc
  dqtdx[0]  = uw[0]/(n1*sin2phiu);
  dqtdx[1]  = uw[1]/(n1*sin2phiu);
  dqtdx[2]  = uw[2]/(n1*sin2phiu);
#if 0
//according to Helgaker, but doesn't work
  dqtdx[3]   = -uw[0]/(n1*sin2phiu) + ( uw[0]*cosphiu/(n2*sin2phiu) - vw[0]*cosphiv/(n2*sin2phiv) );
  dqtdx[4]   = -uw[1]/(n1*sin2phiu) + ( uw[1]*cosphiu/(n2*sin2phiu) - vw[1]*cosphiv/(n2*sin2phiv) );
  dqtdx[5]   = -uw[2]/(n1*sin2phiu) + ( uw[2]*cosphiu/(n2*sin2phiu) - vw[2]*cosphiv/(n2*sin2phiv) );
  dqtdx[6]   =  vw[0]/(n3*sin2phiv) - ( uw[0]*cosphiu/(n2*sin2phiu) - vw[0]*cosphiv/(n2*sin2phiv) );
  dqtdx[7]   =  vw[1]/(n3*sin2phiv) - ( uw[1]*cosphiu/(n2*sin2phiu) - vw[1]*cosphiv/(n2*sin2phiv) );
  dqtdx[8]   =  vw[2]/(n3*sin2phiv) - ( uw[2]*cosphiu/(n2*sin2phiu) - vw[2]*cosphiv/(n2*sin2phiv) );
#endif
#if 1
  dqtdx[3]   = -uw[0]/(n1*sin2phiu) + ( uw[0]*cosphiu/(n2*sin2phiu) + vw[0]*cosphiv/(n2*sin2phiv) );
  dqtdx[4]   = -uw[1]/(n1*sin2phiu) + ( uw[1]*cosphiu/(n2*sin2phiu) + vw[1]*cosphiv/(n2*sin2phiv) );
  dqtdx[5]   = -uw[2]/(n1*sin2phiu) + ( uw[2]*cosphiu/(n2*sin2phiu) + vw[2]*cosphiv/(n2*sin2phiv) );
  dqtdx[6]   =  vw[0]/(n3*sin2phiv) - ( uw[0]*cosphiu/(n2*sin2phiu) + vw[0]*cosphiv/(n2*sin2phiv) );
  dqtdx[7]   =  vw[1]/(n3*sin2phiv) - ( uw[1]*cosphiu/(n2*sin2phiu) + vw[1]*cosphiv/(n2*sin2phiv) );
  dqtdx[8]   =  vw[2]/(n3*sin2phiv) - ( uw[2]*cosphiu/(n2*sin2phiu) + vw[2]*cosphiv/(n2*sin2phiv) );
#endif
  dqtdx[9]   = -vw[0]/(n3*sin2phiv);
  dqtdx[10]  = -vw[1]/(n3*sin2phiv);
  dqtdx[11]  = -vw[2]/(n3*sin2phiv);

//  for (int i=0;i<12;i++)
//    dqtdx[i] = dqtdx[i]/10;
//  for (int i=0;i<12;i++)
//    printf(" dqtdx[%i]: %1.4f \n",i,dqtdx[i]);

  delete [] u;
  delete [] v;
  delete [] w;
  delete [] uw;
  delete [] vw;

  return;
}


int ICoord::bmatp_to_U()
{
  printf(" in bmatp_to_U \n");
  fflush(stdout);
  int len = nbonds+nangles+ntor;
  int N3 = 3*natoms;
  int max_size_ic = len;
  int size_xyz = N3;
  printf(" bmatp_to_U. nbonds: %2i nangles: %2i ntorsions: %2i  total: %3i \n",nbonds,nangles,ntor,len);

  int len_d;
  double* e = new double[len];
  double* U = new double[len*len]; //why was this turned off?
  double* tmp = new double[len*N3];
  for (int i=0;i<len*N3;i++)
    tmp[i] = bmatp[i];

  double* G = new double[len*len];
#if 1
  mat_times_mat_bt(G,bmatp,bmatp,len,len,N3);
#else
  for (int i=0;i<len*len;i++) G[i] = 0.;
  for (int i=0;i<len;i++) 
  for (int j=0;j<len;j++)
  for (int k=0;k<N3;k++)
    G[i*len+j] += bmatp[i*N3+k]*bmatp[j*N3+k];
#endif

#if 1
  printf(" G: \n");
  for (int i=0;i<len;i++)
  {
    for (int j=0;j<len;j++)
      printf(" %1.2f",G[i*len+j]);
    printf("\n");
  }
#endif

  printf("\n using diagonalize(G) \n");
//  printf(" before diagonalize: mkl_threads: %i \n",mkl_get_max_threads());
//  fflush(stdout);
  Diagonalize(G,e,len);
  len_d = N3-6;
  printf(" after diagonalize \n");
//  fflush(stdout);

#if 1
  printf(" eigenvalues:");
  for (int i=0;i<len;i++)
    printf(" %10.8f",e[i]);
  printf("\n");
#endif
#if 1
  int lowev = 0;
  for (int i=0;i<len_d;i++)
  if (e[len-1-i]<0.001)
  {
#if 1
    printf(" small ev: %10.8f \n",e[len-1-i]);
    int i1 = len-1-i;
  //  for (int j=0;j<len;j++)
  //    printf(" %8.5f",G[i1*len+j]);
  //  printf("\n");
#endif
    lowev++;
  }
  if (lowev>0)
    printf(" lowev: %i",lowev);
  len_d -= lowev;
  if (lowev>3)
  {
    printf("\n\n ERROR: optimization space less than 3N-6 DOF \n");
    printf("  probably need more IC's \n");
    printf("  check fragmentation or linear angles \n");
    exit(-1);
  }
#endif

  int redset = len - len_d;
  for (int i=0;i<len_d;i++)
  for (int j=0;j<len;j++)
    Ut[i*len+j] = G[(i+redset)*len+j];
  for (int i=0;i<redset;i++)
  for (int j=0;j<len;j++)
    Ut[(len_d+i)*len+j] = G[i*len+j];

#if 1
  printf(" Ut eigenvalues:");
  for (int i=0;i<len_d;i++)
    printf(" %1.4f",e[len-1-i]);
  printf(" \n");
#endif
#if 1
  if (lowev)
  {
    for (int i=0;i<nangles;i++)
    if (anglev[i]>175.)
      printf(" angle %i: %i %i %i: %1.1f \n",i+1,angles[i][0],angles[i][1],angles[i][2],anglev[i]);
    printf(" \n");
  }
#endif

#if 1
  //printf(" checking orthonormality of U vectors \n");
  trans(U,Ut,len,len);
  double dot;
  for (int i=0;i<len;i++)
  for (int j=0;j<len;j++)
  {
    dot = 0.;
    for (int k=0;k<len;k++)
      dot += U[k*len+i]*U[k*len+j];
    if (i!=j && abs(dot)>0.0001)
      printf(" WARNING: dot of %i %i: %1.3f \n",i,j,dot);
  }
#endif

#if 1
  printf(" Printing %i nonredundant (column) vectors of U \n",len_d);
  for (int i=0;i<len;i++)
  {
    for (int j=0;j<len_d;j++)
      printf(" %2.3f",U[i*len+j]);
    printf("\n");
  }
#endif

#if 0
  double* weights = new double[len];
  for (int i=0;i<len;i++) weights[i] = 0.;
  for (int i=0;i<len;i++)
  for (int j=0;j<nbonds;j++)
    weights[i] += Ut[i*len+j]*Ut[i*len+j];
  for (int i=0;i<len;i++)
    printf(" coord %i has bond weight %1.1f \n",i,weights[i]*100);
  double tweight = 0.;
  for (int i=0;i<len_d;i++)
    tweight += weights[i];
  printf(" total in nonred set: %1.3f \n",tweight*100);
 
  delete [] weights;

#endif

#if 0
//CPMZ previous
  printf("\n using SVD \n");
  //double* tmp2 = new double[len*N3];
  //trans(tmp2,tmp,N3,len);
  //SVD(tmp2,U,e,N3,len);
  SVD(tmp,U,e,len,N3);
  trans(Ut,U,len,len);
#endif

//  if (len_d>N3-6 && len_d>1) len_d = N3-6;
//  else if (len_d<N3-6 && len_d>1) len_d = N3-6;
  nicd = len_d;

  for (int i=0;i<ntor;i++)
    torv0[i] = torv[i];

  for (int i=0;i<len*len;i++)
    Ut0[i] = Ut[i];
  nicd0 = nicd;

  delete [] tmp;
  delete [] G;
  delete [] e;
  delete [] U;

  return 0;
}


int ICoord::bmat_create() 
{
  printf(" in bmat_create() \n");
 // fflush(stdout);

  int len = nbonds+nangles+ntor;
  int N3 = 3*natoms;
  int max_size_ic = len;
  int size_xyz = N3;
  
  int len_d = nicd0;

  printf(" determining q in delocalized internals \n");
  printf(" nicd: %i \n",nicd);
  update_ic();
  for (int i=0;i<len_d;i++)
    q[i] = 0.;

#if 0
  for (int i=0;i<len_d;i++)
    for (int j=0;j<nbonds;j++)
      q[i] += Ut[len*i+j]*bondd[j];
  for (int i=0;i<len_d;i++)
    for (int j=0;j<nangles;j++)
      q[i] += Ut[len*i+nbonds+j]*anglev[j]*3.14159/180;
  for (int i=0;i<len_d;i++)
    for (int j=0;j<ntor;j++)
      q[i] += Ut[len*i+nbonds+nangles+j]*torv[j]*3.14159/180;
#endif

#if 1
  for (int i=0;i<len_d;i++)
    for (int j=0;j<nbonds;j++)
      q[i] += Ut[len*i+j]*distance(bonds[j][0],bonds[j][1]);
  for (int i=0;i<len_d;i++)
    for (int j=0;j<nangles;j++)
      q[i] += Ut[len*i+nbonds+j]*angle_val(angles[j][0],angles[j][1],angles[j][2])*3.14159/180;
  for (int j=0;j<ntor;j++) torfix[j] = 0.;
  for (int j=0;j<ntor;j++)
  {
    double tordiff = torv0[j] - torsion_val(torsions[j][0],torsions[j][1],torsions[j][2],torsions[j][3]);
    if (tordiff>180.)
      torfix[j] = 360.;
    else if (tordiff<-180)
      torfix[j] = -360.;
    else torfix[j] = 0;
//    if (abs(tordiff)>180)
//      printf(" tordiff: %1.1f, effective tor_val: %1.1f ",tordiff,torsion_val(torsions[j][0],torsions[j][1],torsions[j][2],torsions[j][3])+torfix[j]);
  }
  //printf(" torfix: "); for (int j=0;j<ntor;j++) printf(" %1.1f",torfix[j]); printf("\n");
  for (int i=0;i<len_d;i++)
    for (int j=0;j<ntor;j++)
      q[i] += Ut[len*i+nbonds+nangles+j]*(torfix[j]+torsion_val(torsions[j][0],torsions[j][1],torsions[j][2],torsions[j][3]))*3.14159/180;
#endif

#if 1
  printf(" printing q: \n");
  for (int i=0;i<len_d;i++)
    printf(" %1.4f",q[i]);
  printf(" \n");
#endif
 
#if 0
  printf(" verifying q \n");
  double* q0 = new double[len];
  for (int i=0;i<len;i++)
    q0[i] = 0.;
  for (int i=0;i<len;i++)
    for (int j=0;j<len_d;j++)
      q0[i] += Ut[j*len+i]*q[j];
#if 0
  printf(" printing q0: \n");
  for (int i=0;i<len;i++)
    printf(" %1.2f",q0[i]);
  printf(" \n");
#endif
  printf(" printing q0 vs. bonds: \n");
  for (int i=0;i<nbonds;i++)
    printf(" %1.2f",q0[i]-bondd[i]);
  printf(" \n");

  delete [] q0;
#endif

  printf(" now making bmat in delocalized internals (len: %i len_d: %i) \n",len,len_d);
#if 1
  mat_times_mat(bmat,Ut,bmatp,len_d,N3,len);
#else
  for (int i=0;i<len_d*N3;i++) bmat[i] = 0.;
  for (int i=0;i<len_d;i++)
  for (int j=0;j<N3;j++)
  for (int k=0;k<len;k++)
    bmat[i*N3+j] += Ut[i*len+k]*bmatp[k*N3+j];
#endif

#if 1
  printf(" printing bmat in coordinates U \n");
  for (int i=0;i<len_d;i++)
  {
    for (int j=0;j<natoms;j++)
    {
      for (int k=0;k<3;k++)
        printf(" %1.3f",bmat[i*N3+3*j+k]);
    }
    printf(" \n");
  }
#endif

  double* bbt = new double[len_d*len_d];
  double* bbti = new double[len_d*len_d];

#if 1
  mat_times_mat_bt(bbt,bmat,bmat,len_d,len_d,N3);
#else
  for (int i=0;i<len_d*len_d;i++) bbt[i] = 0.;
  for (int i=0;i<len_d;i++)
  for (int j=0;j<len_d;j++)
  for (int k=0;k<N3;k++)
    bbt[i*len_d+j] += bmat[i*N3+k]*bmat[j*N3+k];
#endif

  for (int i=0;i<len_d*len_d;i++)
    bbti[i] = bbt[i];

  //need to invert bbt, then bbt-1 * bmat = bt-1
 // printf(" before invert bbti \n");
 // fflush(stdout);
  Invert(bbti,len_d);

#if 0
  //Checked inverse, it is okay only when 3N-6 vectors are present

  double* tmp2 = new double[len_d*len_d];
  for (int i=0;i<len_d*len_d;i++)
    tmp2[i] = 0.;
  for (int i=0;i<len_d;i++)
  for (int j=0;j<len_d;j++)
  for (int k=0;k<len_d;k++)
    tmp2[i*len_d+j] += bbti[i*len_d+k]*bbt[k*len_d+j];

  printf(" debug: bbti*bbt diagonals \n");
  for (int i=0;i<len_d;i++)
    printf(" %1.2f",tmp2[i*len_d+i]);
#if 0
  printf(" debug: bbti*bbt \n");
  for (int i=0;i<len_d;i++)
  {
    for (int j=0;j<len_d;j++)
      printf(" %1.3f",tmp2[i*len_d+j]);
    printf("\n");
  }
#endif
  printf("\n");
  delete [] tmp2;
#endif

 //printf(" bmatti formation \n");
#if 1
  mat_times_mat(bmatti,bbti,bmat,len_d,N3,len_d);
#else
  for (int i=0;i<len_d*N3;i++)
    bmatti[i] = 0.;
  for (int i=0;i<len_d;i++)
  for (int j=0;j<N3;j++)
  for (int k=0;k<len_d;k++)
    bmatti[i*N3+j] += bbti[i*len_d+k]*bmat[k*N3+j];
#endif

// printf(" dealloc \n");
  delete [] bbt;
  delete [] bbti;

  return 0;
}


int ICoord::grad_to_q() 
{


  printf(" in grad_to_q \n");

  int N3 = 3*natoms;
  int len_d = nicd0;
  int len0 = nbonds+nangles+ntor;
  for (int i=0;i<len_d;i++)
    pgradq[i] = gradq[i];
  for (int i=0;i<len_d;i++)
    gradq[i] = 0.0;


  for (int i=0;i<len_d;i++)
  for (int j=0;j<N3;j++)
    gradq[i] += bmatti[i*N3+j] * grad[j];

  gradrms = 0.;
  for (int i=0;i<nicd;i++)
    gradrms+=gradq[i]*gradq[i];
  gradrms = sqrt(gradrms/nicd);

#if 0
// for Hessian update
  for (int i=0;i<len0;i++) pgradqprim[i] = gradqprim[i];
  for (int i=0;i<len0;i++) gradqprim[i] = 0.;
  for (int i=0;i<len0;i++)
  for (int j=0;j<len_d;j++)
    gradqprim[i] += Ut[j*len0+i]*gradq[j];
#endif

  print_gradq();


  return 0;
}
void ICoord::print_gradq()
{

  printf(" Gradient in delocalized IC:\n");
  int len = nicd0;
  for (int i=0;i<len;i++)
    printf(" %12.10f",gradq[i]);
  printf("\n");

  return;
}


void ICoord::project_grad()
{
	
  int len = nbonds+nangles+ntor;
	printf(" In project_grad\n"); 

  nicd = nicd0;
	double* gradq0 = new double[nicd];
	for (int i=0;i<nicd;i++)
		gradq0[i] = gradq[i];

  //take gradq vector, project it out of all Ut
  double norm = 0.;
  for (int i=0;i<len;i++)
    norm += gradq0[i]*gradq0[i];
  norm = sqrt(norm);
  for (int j=0;j<len;j++)
    gradq0[j] = gradq0[j]/norm;

  double* dots = new double[len];
  for (int i=0;i<len;i++) dots[i] =0.;

  double* Cn = new double[len];
  for (int i=0;i<len;i++) Cn[i] =0.;

  for (int i=0;i<len;i++)
  for (int j=0;j<len;j++)
    dots[i] += gradq0[j]*Ut[i*len+j];
 
  for (int i=0;i<nicd0;i++) //CPMZ fix? right now subspace projection
  for (int j=0;j<len;j++)
    Cn[j] += dots[i]*Ut[i*len+j];

  norm = 0.;
  for (int i=0;i<len;i++)
    norm += Cn[i]*Cn[i];
  norm = sqrt(norm);
  //printf(" Cn norm: %1.2f \n",norm);
  for (int j=0;j<len;j++)
    Cn[j] = Cn[j]/norm;

#if 0
	//Save this to take the overlap with (before Schmidt'ing)
	for (int j=0;j<len;j++)
		C[j]=Cn[j];
#endif
#if 1
  printf(" projected gradient constraint: ");
  for (int i=0;i<len;i++)
    printf(" %1.3f",gradq0[i]);
  printf("\n");
#endif
  for (int i=0;i<len;i++) dots[i] =0.;
  for (int i=0;i<len;i++) 
  for (int j=0;j<len;j++)
    dots[i] += Cn[j]*Ut[i*len+j];

//  for (int i=0;i<nicd0;i++)
//    printf(" dots[%i]: %1.2f \n",i,dots[i]);

  for (int i=0;i<nicd0;i++)
  {
    if (i!=nicd0-1)
    for (int j=0;j<len;j++)
      Ut[i*len+j] -= dots[i] * Cn[j];

    for (int k=0;k<i;k++)
    {
      double dot2 = 0.;
      for (int j=0;j<len;j++)
        dot2 += Ut[i*len+j] * Ut[k*len+j];

      for (int j=0;j<len;j++)
        Ut[i*len+j] -= dot2 * Ut[k*len+j];
    } // loop k over previously formed vectors
 
    double norm = 0.;
    for (int j=0;j<len;j++)
      norm += Ut[i*len+j] * Ut[i*len+j];
    norm = sqrt(norm);
    if (abs(norm)<0.00001) norm = 1;
    if (abs(norm)<0.00001) printf(" WARNING: small norm: %1.7f \n",norm);
    for (int j=0;j<len;j++)
      Ut[i*len+j] = Ut[i*len+j]/norm;
  }

  for (int j=0;j<len;j++)
    Ut[nicd*len+j] = Cn[j];
#if 0
  printf(" printing Cn vs. Ut[nicd*len]\n");
  for (int j=0;j<len;j++)
    printf(" %1.2f/%1.2f\n",Cn[j],Ut[nicd*len+j]);
#endif
#if 1
  printf(" printing orthonormalized vectors \n");
  for (int i=0;i<nicd0;i++)
  {
    for (int j=0;j<len;j++)
      printf(" %1.3f",Ut[i*len+j]);
    printf("\n");
  }
#endif

  delete [] dots;
  delete [] Cn;
	delete [] gradq0; 

  return;
}

int ICoord::ic_to_xyz() 
{
	cout << "In ic_to_xyz" << endl;
  int MAX_STEPS = 10; //was 6

  int success = 1;

  int N3 = 3*natoms;
  //int len0 = min(N3,nbonds+nangles+ntor);
  int len = nicd0;
  double** xyzall = new double*[MAX_STEPS+2];
  for (int i=0;i<MAX_STEPS+2;i++)
    xyzall[i] = new double[N3];
  double* magall = new double[MAX_STEPS+2];
  for (int i=0;i<MAX_STEPS+2;i++)
    magall[i] = 100.;
  double* xyz1 = new double[N3];
  double* xyzd = new double[N3];
  //double* xyzd0 = new double[N3];
  double* btit = new double[N3*len];
  double* dq = new double[len];
  double* qn = new double[len]; //target IC values

  for (int i=0;i<N3;i++)
    xyzall[0][i] = coords[i];

  update_ic();
  bmatp_create();
  bmat_create();

  for (int i=0;i<len;i++)
    dq[i] = dq0[i];
  for (int i=0;i<len;i++)
    qn[i] = q[i] + dq[i];

#if 1
    printf(" qn:");
    for (int i=0;i<len;i++)
      printf(" %1.4f",qn[i]);
    printf("\n");
#endif
#if 0
  printf(" dq at start: \n");
  for (int i=0;i<len;i++)
    printf(" %1.4f",dq[i]);
  printf("\n");
#endif

  double MAX_MOVE = 0.75;
  double mag = 0.;
  double magp = 100.;
  double SCALEBT = 1.5;
  for (int n=0;n<MAX_STEPS;n++) 
  {
    trans(btit,bmatti,N3,len);

    for (int i=0;i<N3;i++) xyzd[i] = 0.;
    for (int i=0;i<N3;i++)
    for (int j=0;j<len;j++) 
      xyzd[i] += btit[len*i+j] * dq[j];

    mag = 0.;
    for (int i=0;i<N3;i++)
      mag += xyzd[i]*xyzd[i];

#if 0
   //was on, why does it exist?
    if (mag>1.*natoms)
    { 
      printf(" rsc");
      double rsc=1.*natoms/mag;
      for (int i=0;i<N3;i++)
        xyzd[i] = xyzd[i]*rsc;
    }
#endif
#if 0
    for (int i=0;i<N3;i++)
      if (abs(xyzd[i])>MAX_MOVE)
        xyzd[i] = sign(xyzd[i])*MAX_MOVE;
#endif

    for (int i=0;i<N3;i++)
      xyz1[i] = coords[i] + xyzd[i]/SCALEBT;

    //if(n==0)
    //  printf(" diff in xyz mag (start) is: %1.6f \n",sqrt(mag));
    for (int i=0;i<N3;i++)
      xyzall[n+1][i] = xyz1[i];
    magall[n] = mag;

    if (mag>magp)
      SCALEBT *= 1.5;
    magp = mag;

    for (int i=0;i<N3;i++)
      coords[i] = xyz1[i];
    update_ic();
    bmatp_create();
    bmat_create();

    for (int i=0;i<len;i++)
      dq[i] = qn[i] - q[i];

#if 0
    printf(" dq: \n");
    for (int i=0;i<len;i++)
      printf(" %1.4f",dq[i]);
    printf("\n");
#endif

    if (mag<0.00005) break;
  }

#if 0
  for (int i=0;i<natoms;i++)
    printf(" dX: %1.3f %1.3f %1.3f \n",xyzd[3*i+0],xyzd[3*i+1],xyzd[3*i+2]);
  printf("\n");
#endif

  //printf(" diff in xyz mag (end) is: %1.4f \n",sqrt(mag)); 
  double MAXMAG = 0.025*natoms;
  if (sqrt(mag)>MAXMAG)
  {
    //printf(" WARNING: diff in xyz mag (end) is: %1.4f, using first step, mag: %1.4f \n",sqrt(mag),sqrt(mag0));
    ixflag++;
    double maglow = 100.;
    int nlow = -1;
    for (int n=0;n<MAX_STEPS+2;n++)
    if (magall[n]<maglow)
    {
      nlow = n+1;
      maglow = magall[n];
    }
    if (maglow<MAXMAG)
    {
      for (int i=0;i<N3;i++)
        coords[i] = xyzall[nlow][i];
      printf("W(%6.5f/%i)",maglow,nlow);
    }
    else
    {
      for (int i=0;i<N3;i++)
        coords[i] = xyzall[0][i];
      printf("Wf(%6.5f/%i)",maglow,nlow);
      success = 0;
    }
  }
  else if (ixflag>0)
    ixflag = 0;

  delete [] btit;
  delete [] xyz1;
  delete [] xyzd;
  for (int i=0;i<MAX_STEPS+2;i++)
    delete [] xyzall[i];
  delete [] magall;
  delete [] xyzall;
  delete [] dq;
  delete [] qn;

  return success;
}



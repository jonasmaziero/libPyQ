#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

Erho(a, b, D, t)
{
  double a, b, D, t

  return E
}

double concurrence(int *d, double _Complex rho[][*d]){  // Returns the entanglement concurrence, for two-qubit states
  void ArrayDisplayC();
  int j, ds = 2;
  double _Complex R[*d][*d], rhot[*d][*d], rhoc[*d][*d], A[*d][*d];
  //double _Complex s2[2][2] = {{0, -I},{I,0}}, kp[*d][*d];  void kronecker();  kronecker(&ds, &ds, &ds, &ds, s2, s2, kp);
  double _Complex kp[4][4] = {{0,0,0,-1},{0,0,1,0},{0,1,0,0},{-1,0,0,0}};
  void matconj();  matconj(d, d, rho, rhoc);
  void matmulC();  matmulC(d, d, d, kp, rhoc, A);  matmulC(d, d, d, A, kp, rhot);  matmulC(d, d, d, rho, rhot, R);
  double egval[*d];  double _Complex egvalR[*d];
  char jobz = 'N';  void lapacke_zgeev();  lapacke_zgeev(&jobz, d, R, egvalR);
  for(j = 0; j < (*d); j++){
    egval[j] = creal(egvalR[j]);
    //printf("%f \t", egval[j]);
  }
  double maxArray1D(); double egvalMax;  egvalMax = maxArray1D(d, egval); // printf("%f \n", egvalMax);
  double cc = 2.0*sqrt(egvalMax) - sqrt(egval[1]) - sqrt(egval[2]) - sqrt(egval[3]) - sqrt(egval[4]);
  double conc = 0.0;  if(cc > 0.0){conc = cc;}
  //printf("%f \n",conc);
  return conc;
}

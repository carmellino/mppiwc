#include <stdio.h>
#include <math.h>
#include <string.h>

// rozd. siatki
#define LX 128
#define LY 128

// stale modelu
#define F 5
#define cs_sq (1./3.)
#define wC (1./3.)
#define wNC (1./6.)

enum LatVel
{ _C, _N, _S, _W, _E };

// x od C (central), E, N, W, S
double cx[F] = {0, 0, 0, -1, 1};
// y od C, E, N, W, S
double cy[F] = {0, 1, -1, 0, 0};

// wagi
double w[F] = {wC, wNC, wNC, wNC, wNC};

// funkcje dystrybucji
double cells[LX*LY*F];
double cells_tmp[LX*LY*F];
// x w nawiasie bo makra
#define ARR(x,y,f) ( f + F*(x) + F*(y)*LX)

double tau;
double uLB;
double vLB;

double fEq(int i, double phi, double u, double v)
{
  return w[i] * phi * (1 + (cx[i]*u + cy[i]*v) / cs_sq);
}

void collideAndStream(double* c, double* ct)
{
  int x,y,i,xp,xm,yp,ym;
  double Fstart[F];
  double phi;

  for(x=0; x<LX; x++)
  {
    for(y=0; y<LY; y++)
    {
      // 1. wielkosci makroskopowe
      phi = 0;
      for(i=0; i<F; i++)
      {
        // (tablica jednowymiarowa zamiast trojwymiarowej)
        phi += c[ARR(x,y,i)];
      }

      // 2. kolizja
      for(i=0; i<F; i++)
      {
        Fstart[i] = (1 - 1./tau) * c[ARR(x,y,i)] + 1./tau * fEq(i, phi, uLB, vLB);
      }

      // periodic BC
      xp = ((x == LX-1) ? 0 : x+1);
      xm = ((x == 0) ? LX-1 : x-1);
      yp = ((y == LY-1) ? 0 : y+1);
      ym = ((y == 0) ? LY-1 : y-1);

      ct[ARR(x,y,_C)] = Fstart[_C];
      ct[ARR(x,yp,_N)] = Fstart[_N];
      ct[ARR(x,ym,_S)] = Fstart[_S];
      ct[ARR(xm,y,_W)] = Fstart[_W];
      ct[ARR(xp,y,_E)] = Fstart[_E];
    }
  }
}

// warunki poczatkowe
void setIC(double* c, double phii, double uLBi, double vLBi)
{
  int x,y,i;
  for(x=0; x<LX; x++)
  {
    for(y=0; y<LY; y++)
    {
      for(i=0; i<F; i++)
      {
        c[ARR(x,y,i)] = fEq(i, phii, uLBi, vLBi);
      }
    }
  }

  for(i=0; i<F; i++)
  {
    c[ARR(LX/2,LY/2,i)] = fEq(i, 1., uLBi, vLBi);
  }
}


// ParaView
// VisIt

void dumpStateVTK(char *fname) {
 FILE *fp;
 int x,y;

 fp = fopen(fname, "w");
 fprintf(fp,"# vtk DataFile Version 2.0\n");
 fprintf(fp,"2D-ADE data file \n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET RECTILINEAR_GRID\n");
 fprintf(fp,"DIMENSIONS %d %d %d\n", LX, LY, 1);
 fprintf(fp,"X_COORDINATES %d int\n", LX);
 for (x = 0; x < LX; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Y_COORDINATES %d int\n", LY);
 for (x = 0; x < LY; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Z_COORDINATES 1 int\n");
 fprintf(fp, "0\n");
 fprintf(fp,"POINT_DATA %d \n", LX*LY);
 fprintf(fp,"SCALARS temperature double 1\n");
 fprintf(fp,"LOOKUP_TABLE default\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++) {
    double rho = 0;
    int f;
    for (f = 0; f < F; f++) rho += cells[ARR(x,y,f)];
    fprintf(fp, "%e \n", rho);
  }
 fprintf(fp,"\n");
 fclose(fp);
}




int main(int argc, char** argv)
{
  // liczba iteracji
  long int iter = 0;
  int ITERMAX = 1e3;

  // warunki poczatkowe
  double phiI = 0;
  double uLBI = 0;
  double vLBI = 0;
  double DlB = 0.01;

  tau = DlB / cs_sq + .5;

  setIC(cells, phiI, uLBI, vLBI);
  dumpStateVTK("state0.vtk");
  

  do
  {
    if(iter%2)
    {
      collideAndStream(cells_tmp, cells);
    }
    else 
    {
      collideAndStream(cells, cells_tmp);
    }
    iter ++;
  }
  while (iter < ITERMAX);

  // dumpStateVTK("state.vtk");

  return 0;
}
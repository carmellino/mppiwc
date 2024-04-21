#include <stdio.h>
#include <math.h>
#include <string.h>

// rozd. siatki
#define LX 256
#define LY 256

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
double cells[LX * LY * F];
double cells_tmp[LX * LY * F];
// x w nawiasie bo makra
#define ARR(x,y,f) ( f + F * (x) + F * (y) * LX)

double uLBm[LX * LY];
double vLBm[LX * LY];
double phim[LX * LY];
#define ARRM(x,y) (x + (y) * LX)

double tau;

double fEq(int i, double phi, double u, double v)
{
  // printf("fEq(%d, %lf, %lf, %lf)\n", i, phi, u, v);
  return w[i] * phi * (1 + (cx[i] * u + cy[i] * v) / cs_sq);
}

void collideAndStream(double* c, double* ct)
{
  // printf("collideAndStream()\n");
  int x,y,i,xp,xm,yp,ym;
  double Fstart[F];
  double phi, uLB, vLB;

  for(x=0; x<LX; x++)
    for(y=0; y<LY; y++)
    {
      // 1. wielkosci makroskopowe
      phi = 0;
      for(i=0; i<F; i++)
        // (tablica jednowymiarowa zamiast trojwymiarowej)
        phi += c[ARR(x,y,i)];

      phim[ARRM(x,y)] = phi;
      uLB = uLBm[ARRM(x,y)];
      vLB = vLBm[ARRM(x,y)];

      // 2. kolizja
      for(i=0; i<F; i++)
        Fstart[i] = (1 - 1. / tau) * c[ARR(x,y,i)] + 1. / tau * fEq(i, phi, uLB, vLB);

      // periodic BC
      xp = ((x == LX - 1) ? 0 : x + 1);
      xm = ((x == 0) ? LX - 1 : x - 1);
      yp = ((y == LY - 1) ? 0 : y + 1);
      ym = ((y == 0) ? LY - 1 : y - 1);

      ct[ARR(x,y,_C)] = Fstart[_C];
      ct[ARR(x,yp,_N)] = Fstart[_N];
      ct[ARR(x,ym,_S)] = Fstart[_S];
      ct[ARR(xm,y,_W)] = Fstart[_W];
      ct[ARR(xp,y,_E)] = Fstart[_E];
    }
}

// warunki poczatkowe
void setIC(double* c, double phii, double uLBi, double vLBi)
{
  printf("setIC()\n");
  int x,y,i;
  for(x=0; x<LX; x++)
    for(y=0; y<LY; y++)
    {
      for(i=0; i<F; i++)
      {
        double phi = 0;
        if((pow((x / ((double)LX - 1) - .5), 2.) + pow((y / ((double)LY - 1) - .75), 2.)) <= (.15 * .15) &&
           ((fabs(x/((double)LX-1)-.5) >= .025) || (y / ((double)LY - 1) >= .85)))
          phi++;
        
        double u = uLBi *(.5 - y / ((double)LY - 1) - .5);
        double v = vLBi *(x / ((double)LX - 1) - .5);
        c[ARR(x,y,i)] = fEq(i, 0, uLBi, vLBi);
        uLBm[ARRM(x,y)] = u;
        vLBm[ARRM(x,y)] = v;
        phim[ARRM(x,y)] = phi;
      }
    }

  // to tu zeruje wszystko
  for(x=0; x<LX; x++)
    for(y=0; y<LY; y++)
    {
      for(i=0; i<F; i++)
      {
        c[ARR(x,y,i)] = fEq(i, 1., uLBi, vLBi);
      }
    }

  for(x=32-4; x<32+4; x++)
    for(y=64-4; y<64+4; y++)
    {
      for(i=0; i<F; i++)
      {
        c[ARR(x,y,i)] = fEq(i, 1., uLBi, vLBi);
      }
    }

  for(y = 0; y < LY; y++)
    for(i=0; i<F; i++)
    {
      c[ARR(32,LY/2,i)] = fEq(i, 1., uLBi, vLBi);
    }

  printf("setIC(): end\n");
}


// ParaView
// VisIt
void dumpStateVTK(char *fname) {
  printf("dumpStateVTK()\n");
  FILE *fp = fopen(fname, "w");
  int x,y;

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
    for (x = 0; x < LX; x++)
    {
      double rho = 0;
      int f;
      for (f = 0; f < F; f++) rho += cells[ARR(x,y,f)];
      fprintf(fp, "%e \n", rho);
    }

  fprintf(fp,"\n");
  fprintf(fp,"Vectors velocity double\n");

  for(y=0; y<LY; y++)
    for(x=0; x<LX; x++)
      fprintf(fp,"%e %e 0.\n", uLBm[ARRM(x,y)], vLBm[ARRM(x,y)]);
    
  fclose(fp);

  printf("dumpStateVTK(): end\n");
}


int main(int argc, char** argv)
{
  printf("main()\n");

  // liczba iteracji
  long int iter = 0;
  int ITERMAX = 5e3;

  // Ma & Pe def.
  double Ma = 0.01; // zawsze znacząco mniejsza od 1
  double Pe = 100;

  double uLB = Ma * sqrt(cs_sq);
  double vLB = uLB;

  double Dlb = uLB * LX / Pe;

  // warunki poczatkowe
  double phiI = 1;
  double uLBI = uLB;
  double vLBI = vLB;

  tau = Dlb / cs_sq + .5;
  double dt = 1. / ((double)LX) * uLB / 1.;
  
  ITERMAX = 2 * M_PI / dt;
  setIC(cells, phiI, uLBI, vLBI);
  dumpStateVTK("state0.vtk");

  int percentage = -1;
  do
  {
    int newPercentage = iter * 100 / ITERMAX;
    if(newPercentage > percentage)
    {
      percentage = newPercentage;
      printf("main(): loop: %d%%\n", percentage);
    }
      
    
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

  dumpStateVTK("state.vtk");

  return 0;
}
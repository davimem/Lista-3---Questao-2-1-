#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** LeMatriz(char *nome,int *m, int *n)
{
  FILE *arq;
  int i,j;
  double **a;
  
  arq=fopen(nome,"r");
  
  fscanf(arq,"%d",m);
  fscanf(arq,"%d",n);

  a=(double **) malloc(sizeof(double *) * *m);
  
  for(i=0; i<*m; i++) a[i]=(double *)malloc(*n * sizeof(double));
  
  for(i=0; i<*m; i++)
  {
    for(j=0;j<*n; j++)
    {
      fscanf(arq,"%lf",&a[i][j]);
    }
  }
  return a;
}

double* LeVetor(char *nome,int *m)
{
  FILE *fp;
  int i;
  double *v;
  
  fp=fopen(nome,"r");
  
  fscanf(fp,"%d",m);

  v=malloc(sizeof(double) * (*m));
  
  for(i=0;i<*m; i++)
  {
    fscanf(fp,"%lf",&v[i]);
  }
  
  return v;
}

void ImprimeVetor(double *V, int m)
{
  int i;
  for(i=0; i<m; i++) printf("%g\t",V[i]);
  puts("");
}

double NormaVetor(double *v,int m,double p)
{
  double norma=0;
  int i;
  
  if(p==0) 
  {
    norma=fabs(v[0]);
    for(i=1;i<m;i++)
    {
      if(fabs(v[i])>norma) norma=fabs(v[i]);
    }
  }
  else
  {
    for(i=0; i<m; i++)
    {
      norma+=pow(fabs(v[i]), p);
    }
    norma=pow(norma, 1/p);
  }
  
  return norma;
}

double *MultiplicaMatVet(double **M, int m1, int m2, double *N, int n1)
{
  double *S,t;
  int i, k;
  
  if(m2!=n1)
  {
    printf("A mutiplicacao nao pode ser feita!\n"); 
    return NULL;
  }
  
  S=malloc(m1*sizeof(double));

  for(i=0; i<m1; i++)
  {
    t=0;
    for(k=0; k<m2; k++)
    {
      t+=M[i][k]*N[k];
    }
    S[i]=t;
  }
  
  return S;
}

double ProdInt(double *V, int m, double *W)
{
  int i;
  double t=0;

  for(i=0; i<m; i++)
  {
    t+=V[i]*W[i];    
  }

  return t;
}

double Jacobi(double **M, int m, int n, double *v, int p)
{
  double *aux, *dif, *b, dx=0;
  int i, j;

  b=calloc(m, sizeof(double));
  for(i=0; i<m; i++) b[i]=M[i][n-1];
  
  dif=calloc(m, sizeof(double));
  aux=calloc(m, sizeof(double));
  
  for(i=0; i<m; i++) 
  {
    aux[i]=b[i];
    
    for(j=0; j<n-1; j++)
    {
      if(i==j)continue;
      aux[i]-= (M[i][j]* v[j]);
    }
    aux[i]/=M[i][i];
  }
  
  for(i=0; i<m; i++)
  { 
    dif[i]=fabs(aux[i]-v[i]); 
    v[i]=aux[i];
  }

  dx=NormaVetor(dif, m, p);
  return dx;
}

double Gauss(double **M, int m, int n, double *v, int p)
{
  double *aux, *dif, *b, dx=0;
  int i, j;

  b=calloc(m, sizeof(double));
  for(i=0; i<m; i++) b[i]=M[i][n-1];
  
  dif=calloc(m, sizeof(double));
  aux=calloc(m, sizeof(double));
  
  for(i=0; i<m; i++) 
  {
    aux[i]=v[i];
    v[i]=b[i];
    
    for(j=0; j<n-1; j++)
    {
      if(i==j)continue;
      v[i]-= (M[i][j]* v[j]);
    }
    v[i]/=M[i][i];
  }
  
  for(i=0; i<m; i++)
  { 
    dif[i]=fabs(aux[i]-v[i]); 
  }

  dx=NormaVetor(dif, m, p);
  return dx;
}

void Gradiente(double **M, int m, int n, double *x0, double tol)
{
  double *r, *b, l=0;
  int i, it=1;

  b=calloc(m, sizeof(double));
  for(i=0; i<m; i++) b[i]=M[i][n-1];

  r=calloc(m, sizeof(double));
  
  n-=1;
  
  do{
    r=MultiplicaMatVet(M, m, n, x0, m);
  
    for(i=0; i<m; i++)
    {
      r[i]=b[i]-r[i];    
    }
  
    l=ProdInt(r,m,r)/ProdInt(r,m,(MultiplicaMatVet(M, m, n, r, m)));

    for(i=0; i<m; i++)
    {
      x0[i]+=l*r[i];    
    }

    r=MultiplicaMatVet(M, m, n, x0, m);
  
    for(i=0; i<m; i++)
    {
      r[i]=b[i]-r[i];    
    }
    
    printf("%d ", it);
    for(i=0; i<m; i++) printf(" %g", x0[i]);
    puts("");
    
    it++;
  }while(NormaVetor(r,m,2)>tol);
}

void GradienteConjugado(double **M, int m, int n, double *x0, double tol)
{
  double *r, *d, *b, *Ax, *Ad, rr=0,a=0,bt=0;
  int i, it=1;

  b=calloc(m, sizeof(double));
  for(i=0; i<m; i++) b[i]=M[i][n-1];

  r=calloc(m, sizeof(double));
  d=calloc(m, sizeof(double));
  Ax=calloc(m, sizeof(double));
  Ad=calloc(m, sizeof(double));
  
  n-=1;
  
  Ax=MultiplicaMatVet(M, m, n, x0, m);
  
  for(i=0; i<m; i++)
  {
    r[i]=d[i]=b[i]-Ax[i];    
  }
  
  do{
    Ad=MultiplicaMatVet(M, m, n, d, m);
    rr=ProdInt(r,m,r);
    a=rr/ProdInt(d,m,Ad);

    for(i=0; i<m; i++)
    {
      x0[i]+=a*d[i];
    }

    for(i=0; i<m; i++)
    {
      r[i]-=a*Ad[i];
    }
  
    bt=ProdInt(r,m,r)/rr;
  
    for(i=0; i<m; i++)
    {
      d[i]=r[i]+bt*d[i];    
    }
    
    printf("%d ", it);
    for(i=0; i<m; i++) printf(" %g", x0[i]);
    puts("");
    
    it++;
  }while(NormaVetor(r,m,2)>tol);
}

int main(int argc, char **argv) 
{
  double **M, *v, *b, dx, tol=1e-8;
  int i, m, n, l, it=0, p=2;
  
  M=LeMatriz("m.dat", &m, &n);
  v=LeVetor("v.dat",&l);
  
  GradienteConjugado(M, m, n, v, tol);
  puts("");

  //Gradiente(M, m, n, v, tol);
  //puts("");
  
  /*
  //Gauss Seidel
  do{
    it++;
    dx=Gauss(M, m, n, v, p);
    printf("%d %g", it, dx);
    for(i=0; i<m; i++) printf(" %g",v[i]);
    puts("");
    
  }while (dx>tol);
  */  
  return 0;
}
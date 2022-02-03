#include <iostream>
#include<cassert>
#include<cmath>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<sstream>
using namespace std;

class SIR_KPP
{
 public :
  matrice laplacien0(int,int);
  matrice iden(int);
  matrice laplacien_2d(int,int);
  matrice A(int,int);
  vecteur G_1(int,int,double,vecteur&);
  vecteur G_2(int,int,double,double,vecteur&);
  vecteur G_3(int,int,double,vecteur&);
  vecteur G_KPP(int, int, double,double,vecteur&);
  vecteur RK4_KPP(int,int,int,double,double,double, double,vecteur&,const char*);
  vecteur PHI(int,int,double,double,double,double,vecteur&);
  matrice DG(int,int,double,double,vecteur&);
  matrice DPHI(int,int,double,double,double,double,vecteur&); 
  vecteur EI_KPP(int,int,int,double,double,double,double,vecteur&,double,const char*);   
};

  
// Matrice du laplacien en une dimension :

matrice SIR_KPP::laplacien0(int n,int m)
{
  matrice H(n,n);
  for(int i=0;i<n-1;i++){
    H(i,i)= 4;
    H(i,i+1)= -1;
    H(i+1,i)= -1;
  }
  H(n-1,m-1)= 4;
  H(n-1,m-2)= -1;
  return H;
}

// Matrice identité :

matrice SIR_KPP::iden(int n)
{
  matrice K(n,n);
  for(int i=0;i<n;i++){
    K(i,i)= 1;
  }
  return K;
}

// Matrce du laplacien en 2D :

matrice SIR_KPP::laplacien_2d(int Nx,int Ny)
{
  double dx=(float)1/(Nx+1); double dy=(float)1/(Ny+1);
  matrice P(Nx*Ny,Nx*Ny);
  int k=1;
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++)
      P(i,j)= laplacien0(Nx,Ny)(i,j)/pow(dx,2);
  }
  while(k < Nx){
    for(int i=0;i<Nx;i++){
      for(int j=0;j<Ny;j++)
	P(i+k*Nx,j+k*Ny)= laplacien0(Nx,Ny)(i,j)/pow(dx,2);
    }
    k=k+1;
  }
  for (int i=0;i <Nx*Ny-1;i++){
     if (i > Nx-1)
       P(i,i-Nx)=-(float)1/pow(dy,2);
     if (i <( Nx-1)*(Ny))
       P(i,i+Ny)=-(float)1/pow(dy,2);
   }
 
  P(Nx*Ny-1,Nx*Ny-1)= 4 /pow(dx,2);
  P((Nx*Ny)-1,((Nx-1)*Ny)-1)= -(float)1/pow(dy,2);
  return P;
}

  
// Construction de la matrice A par bloc de l'ennoncé :

matrice SIR_KPP::A(int n,int m)
{
  matrice T(3*n*m,3*n*m);
  for(int i=0;i<n*m;i++){
    for(int j=0;j<n*m;j++){
      T(i,j)=laplacien_2d(n,m)(i,j);
      T(i+n*m,j+n*m)=laplacien_2d(n,m)(i,j);
      T(i+2*n*m,j+2*n*m)=laplacien_2d(n,m)(i,j);
    }
  }

  return T;
}

//construction de la fonction G :

vecteur SIR_KPP::G_1(int Nx,int Ny,double k_s,vecteur& U0)
{
  vecteur s(Nx*Ny);
  for(int i=0;i<Nx*Ny;i++){
    s(i)=-k_s*U0(i)*U0(i+Nx*Ny);
  }
  return s;
}

vecteur SIR_KPP::G_2(int Nx,int Ny,double k_s,double k_i,vecteur& U0)
{
  vecteur s2(Nx*Ny);
  for(int i=0;i<Nx*Ny;i++){
    s2(i)= k_s*U0(i)*U0(i+Nx*Ny) - k_i*U0(i+Nx*Ny);
  }
  return s2;
}

vecteur SIR_KPP::G_3(int Nx,int Ny,double k_i,vecteur& U0)
{
  vecteur s3(Nx*Ny);
  for(int i=0;i<Nx*Ny;i++){
    s3(i)= k_i*U0(i+Nx*Ny);
  }
  return s3;
}
  
vecteur SIR_KPP::G_KPP(int Nx,int Ny,double k_s,double k_i,vecteur& U)
{
  vecteur v(3*Nx*Ny);
  for(int i=0;i<Nx*Ny;i++){
    v(i)=G_1(Nx,Ny,k_s,U)(i);
    v(i+Nx*Ny)=G_2(Nx,Ny,k_s,k_i,U)(i);
    v(i+2*Nx*Ny)=G_3(Nx,Ny,k_i,U)(i);
  }
  return v;
}

// solution schéma RK4 :

vecteur SIR_KPP::RK4_KPP(int Nx,int Ny,int T,double D,double h,double k_s,double k_i,vecteur& U,const char* Nomfichier)
{
  int Nt= T/h;
  vecteur x(3*Nx*Ny);
  vecteur x_new(3*Nx*Ny);
  ofstream file(Nomfichier);
  x=U;int k=0;
  for(int n=0;n < Nt;n++){
    x_new= x - D*h*A(Nx,Ny)*x + h*G_KPP(Nx,Ny,k_s,k_i,x);
    if((n<20) && (n%3==0)){
      for(int i=0;i<Nx*Ny;i++){
	file << x_new(i) << " " << x_new(i+Nx*Ny) << " " << x_new(i+2*Nx*Ny) << endl;
	file << endl << endl;
      }
    }
    x=x_new;
  }
  file.close();
  return x_new;
}

// cas du schema Euler implicite:

vecteur SIR_KPP::PHI(int Nx,int Ny,double h,double D,double k_s,double k_i,vecteur& U)
{
  vecteur x(3*Nx*Ny);
  x= (iden(3*Nx*Ny)+D*h*A(Nx,Ny))*U - h*G_KPP(Nx,Ny,k_s,k_i,U)-U;
  return x;
}
  
matrice SIR_KPP::DG(int Nx,int Ny,double k_s,double k_i,vecteur& U)
{
  matrice m(3*Nx*Ny,3*Nx*Ny);
  for(int i=0;i<Nx*Ny;i++){
    m(i,i)= -k_s*U(i+Nx*Ny);
    m(i,i+Nx*Ny)= -k_s*U(i);
    m(i,i+2*Nx*Ny)= 0;

    m(i+Nx*Ny,i)= k_s*U(i+Nx*Ny);
    m(i+Nx*Ny,i+Nx*Ny)= k_s*U(i)-k_i;
    m(i+Nx*Ny,i+2*Nx*Ny)=0;

    m(i+2*Nx*Ny,i)=0;
    m(i+2*Nx*Ny,i+Nx*Ny)= k_i;
    m(i+2*Nx*Ny,i+2*Nx*Ny)=0;  
  }
  return m;
}

matrice SIR_KPP::DPHI(int Nx,int Ny,double D,double h,double k_s,double k_i,vecteur& U)
{
  matrice m(3*Nx*Ny,3*Nx*Ny);
  m= iden(3*Nx*Ny) + D*h*A(Nx,Ny) - h*DG(Nx,Ny,k_s,k_i,U);
  return m;
}

vecteur SIR_KPP::EI_KPP(int Nx,int Ny,int T,double k_s,double k_i,double h,double D,vecteur& U,double esp,const char* fichier)
{
  int Nt= T/h;
  vecteur L_ei(3*Nx*Ny);
  vecteur x_ei(3*Nx*Ny);
  matrice dph_ei(3*Nx*Ny,3*Nx*Ny);
  matrice LU_ei(3*Nx*Ny,3*Nx*Ny);
  vecteur phi_ei(3*Nx*Ny);             
  vecteur piv_ei(3*Nx*Ny);
  int k=0; int iter_max=15; double cond;
  ofstream file(fichier);
  L_ei=U;
  for (int n=0; n < Nt; n++){
    //cout <<"l= " << L << endl;
    do{
      phi_ei=PHI(Nx,Ny,h,D,k_s,k_i,L_ei);
      dph_ei=DPHI(Nx,Ny,D,h,k_s,k_i,L_ei);
      LU_ei=dph_ei.lu(piv_ei);
      x_ei= L_ei - LU_ei.solvelu(phi_ei,piv_ei);
      if((n<20) && (n%3==0)){
	for(int i=0;i<Nx*Ny;i++){
	  file << x_ei(i) << " " << x_ei(i+Nx*Ny) << " " << x_ei(i+2*Nx*Ny) << endl;
	  file << endl << endl;
	}
      }
      cond=(x_ei - L_ei).norme_2();
      L_ei=x_ei;
      k=k+1;
    }while(k < iter_max && cond > esp && G_KPP(Nx,Ny,k_s,k_i,L_ei).norme_2() > esp && DG(Nx,Ny,k_s,k_i,L_ei).norme_inf()> esp);
  }
  file.close();
  return x_ei;
}


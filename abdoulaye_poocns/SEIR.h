#include <iostream>
#include<cassert>
#include<cmath>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<sstream>
#include"matrice.h"
using namespace std;

class SEIR
{ 
 public :
  vecteur F(double, double, double,double, double, double, vecteur&);
  matrice DF(double, double, double,double, double, vecteur&, double);
  vecteur phi(double, double, double,double, double, double, vecteur&, double);
  matrice dphi(double, double, double,double, double, double, vecteur&, double);
  vecteur L_demi(double, double, double, double, double,double, double, vecteur&);
  vecteur k_2(double, double, double, double, double,double, double, vecteur&);
  vecteur k2_rk(double b,double s,double g,double mu,double nu,double N,double pas,vecteur& L);
  vecteur k3_rk(double,double,double,double,double,double,double,vecteur&);
  vecteur k4_rk(double,double,double,double,double,double,double,vecteur&);
  vecteur RK4(double,double,double,double,double,double,int, double,vecteur&, const char*);
  vecteur Heun(double,double,double,double,double,double,int,double,vecteur&,const char*);
  vecteur Newton_cranck_nicolson(double,double,double,double, double, int, double,double, vecteur&,double,  const char*);
};

vecteur SEIR::F(double b, double s, double g,double mu,double nu,double N, vecteur& L)
{
  vecteur v(4);
  v(0) = - b *L(0)*L(2)/N - nu*L(0) ;
  v(1) = (b*L(0)*L(2)/N) - s*L(1) ;
  v(2) = s*L(1) -(g+mu)*L(2) ;
  v(3) = g*L(2) + nu*L(0);
  return v;
}

matrice SEIR::DF(double b, double s,double g,double mu,double nu, vecteur& L,double N)
{
  matrice p(4,4);
  p(0, 0) = -b*L(2)/N - nu;
  p(0, 1) = 0;
  p(0, 2) = -b*L(0)/N;
  p(0, 3) = 0;
  p(1, 0) = b*L(2)/N;
  p(1, 1) = -s;
  p(1, 2) = b*L(0)/N;
  p(1, 3) = 0;
  p(2, 0) = 0;
  p(2, 1) = s;
  p(2, 2) = -(g+mu);
  p(2, 3) = 0;
  p(3, 0) = nu;
  p(3, 1) = 0;
  p(3, 2) = g;
  p(3, 3) = 0;
  return p;
}
  
vecteur SEIR::phi(double b,double si,double g,double mu, double nu,double N,vecteur& L, double pas){
  vecteur s(4);
  s= L - pas*F(b,si,g,mu,nu,N, L)*0.5 - L - 0.5*pas*F(b,si,g,mu,nu,N, L);
  return s;
}

matrice SEIR::dphi(double b,double s,double g,double mu, double nu, double N,vecteur& L, double pas)
{
  matrice m(4,4);
  m(0, 0) = 1 -pas*DF(b,s,g,mu,nu,L,N)(0,0)*0.5;
  m(0, 1) = -pas*DF(b,s,g,mu,nu,L,N)(0,1)*0.5;
  m(0, 2) = -pas*DF(b,s,g,mu,nu,L,N)(0,2)*0.5;
  m(0, 3) = -pas*DF(b,s,g,mu,nu,L,N)(0,3)*0.5;
  
  m(1, 0) = -pas*DF(b,s,g,mu,nu,L,N)(1,0)*0.5;
  m(1, 1) = 1 -pas*DF(b,s,g,mu,nu,L,N)(1,1)*0.5;
  m(1, 2) = -pas*DF(b,s,g,mu,nu,L,N)(1,2)*0.5;
  m(1, 3) = -pas*DF(b,s,g,mu,nu,L,N)(1,3)*0.5;

  m(2, 0) = -pas*DF(b,s,g,mu,nu,L,N)(2,0)*0.5;
  m(2, 1) = -pas*DF(b,s,g,mu,nu,L,N)(2,1)*0.5;
  m(2, 2) = 1 -pas*DF(b,s,g,mu,nu,L,N)(2,2)*0.5;
  m(2, 3) = -pas*DF(b,s,g,mu,nu,L,N)(2,3)*0.5;

  m(3, 0) = -pas*DF(b,s,g,mu,nu,L,N)(3,0)*0.5;
  m(3, 1) = -pas*DF(b,s,g,mu,nu,L,N)(3,1)*0.5;
  m(3, 2) = -pas*DF(b,s,g,mu,nu,L,N)(3,2)*0.5;
  m(3, 3) = 1 -pas*DF(b,s,g,mu,nu,L,N)(3,3)*0.5;

  return m;
}

// schéma cranck nicolson pour le modèle SEIR :

vecteur SEIR::Newton_cranck_nicolson(double b, double s, double g,double mu,double nu, int Nt, double pas,double N, vecteur& L0,double esp, const char* Nomdufichier)
{
  vecteur L(4);
  vecteur x(4);
  matrice dph(4,4);
  matrice LU(4,4);
  vecteur phi_u(4);             
  vecteur piv(4);
  int k=0; int iter_max=25; double cond;
  ofstream file(Nomdufichier);
  file << L0(0) << " " << L0(1) << " " << L0(2) << " " << L0(3) << endl;
  L=L0;
  for (int n=0; n < Nt; n++){
    //cout <<"l= " << L << endl;
    do{
      phi_u=phi(b,s,g,mu,nu,N,L,pas);
      dph=dphi(b,s,g,mu,nu,N,L,pas);
      LU=dph.lu(piv);
      x= L - LU.solvelu(phi_u,piv);
      file << x(0) << " " << x(1) << " " << x(2) << " " << x(3)  << endl;
      cond=(x - L).norme_2();
      L=x;
      k=k+1;
    }while(k < iter_max && cond > esp && F(b,s,g,mu,nu,N,L).norme_2() > esp && DF(b,s,g,mu,nu,L,N).norme_inf()> esp);
    file << endl << endl;
  }
  file.close();
  return x;
}

// schéma Heun pour le modèle SEIR :

// construction de la fonction L_demi :
vecteur SEIR::L_demi(double b,double s,double g,double mu, double nu,double N,double pas,vecteur& L)
{
  vecteur d(4);
  d=L+pas*F(b,s,g,mu,nu,N,L);
  return d;
}

vecteur SEIR::k_2(double b,double s,double g,double mu,double nu,double N,double pas,vecteur& L)
{
  vecteur k(4);
  vecteur k_demi(4);
  k_demi=L_demi(b,s,g,mu,nu,N,pas,L);
  k=F(b,s,g,mu,nu,N,k_demi);
  return k;
}
vecteur SEIR::Heun(double b,double s,double g,double mu, double nu,double N,int NT, double pas,vecteur& L, const char* Nomfichier)
{
  vecteur y(4);
  vecteur y_new(4);
  ofstream file(Nomfichier);
  file << L(0) << " " << L(1) << " " << L(2) << " " << L(3) << endl;
  y=L;
  for (int n=0;n < NT;n++){
    y_new= y + pas*0.5*(F(b,s,g,mu ,nu,N,y) +k_2(b,s,g,mu,nu,N,pas,y));
    file << y_new(0) << " " << y_new(1) << " " << y_new(2) << " " << y_new(3)  << endl;
    y=y_new;
    file << endl << endl;
  }
  file.close();
  return y_new;
}

// schema RK4 pour le modèle SEIR
vecteur SEIR::k2_rk(double b,double s,double g,double mu,double nu,double N,double pas,vecteur& L)
{
  vecteur u(4);
  vecteur v(4);
  v= L + pas*0.5*F(b,s,g,mu,nu,N,L);
  u= F(b,s,g,mu,nu,N,v);
  return u;
}

vecteur SEIR::k3_rk(double b,double s,double g,double mu,double nu,double N,double pas,vecteur& L)
{
  vecteur r(4);
  vecteur t(4);
  r= L + pas*0.5*k2_rk(b,s,g,mu,nu,N,pas,L);
  t= F(b,s,g,mu,nu,N,r);
  return t;
}

vecteur SEIR::k4_rk(double b,double s,double g,double mu,double nu,double N,double pas,vecteur& L)
{
  vecteur m(4);
  vecteur n(4);
  m= L + pas*k3_rk(b,s,g,mu,nu,N,pas,L);
  n= F(b,s,g,mu,nu,N,m);
  return n;
}

vecteur SEIR::RK4(double b,double s,double g,double mu, double nu,double N,int NT, double pas,vecteur& L, const char* Nom_fichier)
{
  vecteur g1(4);
  vecteur g_new(4);
  ofstream file(Nom_fichier);
  file << L(0) << " " << L(1) << " " << L(2) << " " << L(3) << endl;
  g1=L;
  for (int n=0; n < NT; n++){
    //cout << "g1= " << g1 << endl;
    g_new= g1 + pas*(F(b,s,g,mu,nu,N,g1) + 2*k2_rk(b,s,g,mu,nu,N,pas,g1) + 2*k3_rk(b,s,g,mu,nu,N,pas,g1) + k4_rk(b,s,g,mu,nu,N,pas,g1)) /6;
    file << g_new(0) << " " << g_new(1) << " " << g_new(2) << " " << g_new(3) << endl;
    g1=g_new;
    file << endl << endl;
  }
  file.close();
  return g_new;
}
  
  

  
						    

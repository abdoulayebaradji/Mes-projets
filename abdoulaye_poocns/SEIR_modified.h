#include <iostream>
#include<cassert>
#include<cmath>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<sstream>
using namespace std;

class SEIR_modified
{
 public :
  
  vecteur G(double, double, double, double, double, double, double, double, double, double, double, double, double, vecteur&);

  matrice DG(double, double, double, double, double, double, double, double, double, double, double, double, double, vecteur&);

  vecteur phi_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  matrice dphi_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur L_demi_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur k_2_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur k2_rk_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur k3_rk_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur k4_rk_modified(double, double, double,double,double,double, double, double, double, double, double, double, double, vecteur&, double);

  vecteur RK4_modified(double, double, double,double,double,double, double, double, double, double, double, double, double,int, double,vecteur&, const char*);

  vecteur Heun_modified(double, double, double,double,double,double, double, double, double, double, double, double, double,int, double,vecteur&, const char*);

  vecteur Newton_cranck_modified(double, double, double,double,double,double, double, double, double, double, double, double, double,int, double,vecteur&,double, const char*);
};

vecteur SEIR_modified::G(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L)
{
  vecteur v(8);
  v(0) = -(b*c+c*q*(1-b))*L(0)*(L(2)+teta*L(3))+lambda*L(4) ;
  v(1) = b*c*(1-q)*L(0)*(L(2)+teta*L(3))-s*L(1) ;
  v(2) = s*rho*L(1)-(delta_I+alpha+gama_I)*L(2);
  v(3) = s*(1-rho)*L(1)-gama_A*L(3);
  v(4) = (1-b)*c*q*L(0)*(L(2)+teta*L(3))-lambda*L(4);
  v(5) = b*c*q*L(0)*(L(2)+teta*L(3))-delta_q*L(5);
  v(6) = delta_I*L(2)+delta_q*L(5)-(alpha+gama_H)*L(6);
  v(7) = gama_I*L(2)+gama_A*L(3)+gama_H*L(6);
  return v;
}

matrice SEIR_modified::DG(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L)
{
  matrice p(8,8);
  p(0, 0) = -(b*c+c*q*(1-b))*(L(2)+teta*L(3));
  p(0, 1) = 0;
  p(0, 2) = -(b*c+c*q*(1-b))*L(0);
  p(0, 3) = -teta*(b*c+c*q*(1-b))*L(0);
  p(0, 4) = lambda;
  p(0, 5) = 0;
  p(0, 6) = 0;
  p(0, 7) = 0;
  p(1, 0) = b*c*(1-q)*(L(2)+teta*L(3));
  p(1, 1) = -s;
  p(1, 2) = b*c*(1-q)*L(0);
  p(1, 3) = teta*b*c*(1-q)*L(0);
  p(1, 4) = 0;
  p(1, 5) = 0;
  p(1, 6) = 0;
  p(1, 7) = 0;
  p(2, 0) = 0;
  p(2, 1) = s*rho;
  p(2, 2) = -(delta_I+alpha+gama_I);
  p(2, 3) = 0;
  p(2, 4) = 0;
  p(2, 5) = 0;
  p(2, 6) = 0;
  p(2, 7) = 0;
  p(3, 0) = 0;
  p(3, 1) = s*(1-rho);
  p(3, 2) = 0;
  p(3, 3) = gama_A;
  p(3, 4) = 0;
  p(3, 5) = 0;
  p(3, 6) = 0;
  p(3, 7) = 0;
  p(4, 0) = (1-b)*c*q*(L(2)+teta*L(3));
  p(4, 1) = 0;
  p(4, 2) = (1-b)*c*q*L(0);
  p(4, 3) = teta*(1-b)*c*q*L(0);
  p(4, 4) = -lambda;
  p(4, 5) = 0;
  p(4, 6) = 0;
  p(4, 7) = 0;
  p(5, 0) = b*c*q*(L(2)+teta*L(3));
  p(5, 1) = 0;
  p(5, 2) = b*c*q*L(0);
  p(5, 3) = teta*b*c*q*L(0);
  p(4, 4) = 0;
  p(5, 5) = -delta_q;
  p(5, 6) = 0;
  p(5, 7) = 0;
  p(6, 0) = 0;
  p(6, 1) = 0;
  p(6, 2) = delta_I;
  p(6, 3) = 0;
  p(6, 4) = 0;
  p(6, 5) = delta_q;
  p(6, 6) = -(alpha+gama_H);
  p(6, 7) = 0;
  p(7, 0) = 0;
  p(7, 1) = 0;
  p(7, 2) = gama_I;
  p(7, 3) = gama_A;
  p(7, 4) = 0;
  p(7, 5) = 0;
  p(7, 6) = gama_H;
  p(7, 7) = 0;
  return p;
}
  
vecteur SEIR_modified::phi_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L, double pas)
{
  vecteur o(8);
  o= L - pas*G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)*0.5 - L - 0.5*pas*G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L);
  return o;
}

matrice SEIR_modified::dphi_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  matrice m(8,8);
  m(0, 0) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,0)*0.5;
  m(0, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,1)*0.5;
  m(0, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,2)*0.5;
  m(0, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,3)*0.5;
  m(0, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,4)*0.5;
  m(0, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,5)*0.5;
  m(0, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,6)*0.5;
  m(0, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(0,7)*0.5;

  m(1, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,0)*0.5;
  m(1, 1) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,1)*0.5;
  m(1, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,2)*0.5;
  m(1, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,3)*0.5;
  m(1, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,4)*0.5;
  m(1, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,5)*0.5;
  m(1, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,6)*0.5;
  m(1, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(1,7)*0.5;
  
  m(2, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,0)*0.5;
  m(2, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,1)*0.5;
  m(2, 2) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,2)*0.5;
  m(2, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,3)*0.5;
  m(2, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,4)*0.5;
  m(2, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,5)*0.5;
  m(2, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,6)*0.5;
  m(2, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(2,7)*0.5;

  
  m(3, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,0)*0.5;
  m(3, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,1)*0.5;
  m(3, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,2)*0.5;
  m(3, 3) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,3)*0.5;
  m(3, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,4)*0.5;
  m(3, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,5)*0.5;
  m(3, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,6)*0.5;
  m(3, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(3,7)*0.5;

  
  m(4, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,0)*0.5;
  m(4, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,1)*0.5;
  m(4, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,2)*0.5;
  m(4, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,3)*0.5;
  m(4, 4) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,4)*0.5;
  m(4, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,5)*0.5;
  m(4, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,6)*0.5;
  m(4, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(4,7)*0.5;

  m(5, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,0)*0.5;
  m(5, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,1)*0.5;
  m(5, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,2)*0.5;
  m(5, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,3)*0.5;
  m(5, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,4)*0.5;
  m(5, 5) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,5)*0.5;
  m(5, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,6)*0.5;
  m(5, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(5,7)*0.5;

  m(6, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,0)*0.5;
  m(6, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,1)*0.5;
  m(6, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,2)*0.5;
  m(6, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,3)*0.5;
  m(6, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,4)*0.5;
  m(6, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,5)*0.5;
  m(6, 6) = 1-pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,6)*0.5;
  m(6, 7) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(6,7)*0.5;

  m(7, 0) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,0)*0.5;
  m(7, 1) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,1)*0.5;
  m(7, 2) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,2)*0.5;
  m(7, 3) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,3)*0.5;
  m(7, 4) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,4)*0.5;
  m(7, 5) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,5)*0.5;
  m(7, 6) = -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,6)*0.5;
  m(7, 7) = 1 -pas*DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L)(7,7)*0.5;
  
  return m;
}

// schéma cranck nicolson pour le modèle SEIR Modified :

vecteur SEIR_modified::Newton_cranck_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta,int Nt, double pas, vecteur& L0,double esp, const char* Nomdufichier)
{
  vecteur L(8);
  vecteur x(8);
  matrice dph(8,8);
  matrice LU(8,8);
  vecteur phi_u(8);
  vecteur piv(8);
  int k=0; int iter_max=25; double cond;
  ofstream file(Nomdufichier);
  file << L0(0) << " " << L0(1) << " " << L0(2) << " " << L0(3) << " " << L0(4) << " " <<  L0(5) << " " <<  L0(6) <<  " " << L0(7) << endl;
  L=L0;
  for (int n=0; n < Nt; n++){
    //cout <<"l= " << L << endl;
    do{
      phi_u=phi_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L,pas);
      dph=dphi_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L,pas);
      LU=dph.lu(piv);
      x= L - LU.solvelu(phi_u,piv);
      file << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << x(4) << " " <<  x(5) << " " <<  x(6) << " " <<  x(7) << endl;
      cond=(x - L).norme_2();
      L=x;
      k=k+1;
    }while(k < iter_max && cond > esp && G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L).norme_2() > esp && DG(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L).norme_inf()> esp);
    file << endl << endl;
  }
  file.close();
  return x;
}

// schéma Heun pour le modèle SEIR Modified :

// construction de la fonction L_demi :
vecteur SEIR_modified::L_demi_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  vecteur d(8);
  d=L+pas*G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L);
  return d;
}

vecteur SEIR_modified::k_2_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  vecteur k(8);
  vecteur k_demi(8);
  k_demi=L_demi_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L,pas);
  k=G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,k_demi);
  return k;
}
vecteur SEIR_modified::Heun_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta,int NT, double pas,vecteur& L, const char* Nomfichier)
{
  vecteur y(8);
  vecteur y_new(8);
  ofstream file(Nomfichier);
  file << L(0) << " " << L(1) << " " << L(2) << " " << L(3) << " "  << L(4) << " "  << L(5) << " "  << L(6) << " " << L(7) <<  endl;
  y=L;
  for (int n=0;n < NT;n++){
    y_new= y + pas*0.5*(G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,y) +k_2_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,y,pas));
    file << y_new(0) << " " << y_new(1) << " " << y_new(2) << " " << y_new(3) << " "  << y_new(4) << " "  << y_new(5) << " " << y_new(6) << " " << y_new(7)  << endl;
    y=y_new;
    file << endl << endl;
  }
  file.close();
  return y_new;
}

// schema RK4 pour le modèle SEIR Modified
vecteur SEIR_modified::k2_rk_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  vecteur u(8);
  vecteur v(8);
  v= L + pas*0.5*G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L);
  u= G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,v);
  return u;
}

vecteur SEIR_modified::k3_rk_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  vecteur r(8);
  vecteur t(8);
  r= L + pas*0.5*k2_rk_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L,pas);
  t= G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,r);
  return t;
}

vecteur SEIR_modified::k4_rk_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta, vecteur& L,double pas)
{
  vecteur m(8);
  vecteur n(8);
  m= L + pas*k3_rk_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,L,pas);
  n= G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,m);
  return n;
}

vecteur SEIR_modified::RK4_modified(double c, double b, double q,double s,double lambda,double rho, double delta_I, double delta_q, double gama_I, double gama_A, double gama_H, double alpha, double teta,int NT, double pas,vecteur& L, const char* Nom_fichier)
{
  vecteur g1(8);
  vecteur g_new(8);
  ofstream file(Nom_fichier);
  file << L(0) << " " << L(1) << " " << L(2) << " " << L(3) << " "  << L(4) << " "  << L(5) << " " <<  L(6) << " " <<  L(7) << endl;
  g1=L;
  for (int n=0; n < NT; n++){
    //cout << "g1= " << g1 << endl;
    g_new= g1 + pas*(G(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,g1) + 2*k2_rk_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,g1,pas) + 2*k3_rk_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,g1,pas) + k4_rk_modified(c,b,q,s,lambda,rho,delta_I,delta_q,gama_I,gama_A,gama_H,alpha,teta,g1,pas))/6;
    file << g_new(0) << " " << g_new(1) << " " << g_new(2) << " " << g_new(3) << " " <<  g_new(4) << " "  << g_new(5) << " " << g_new(6) << " " << g_new(7) <<  endl;
    g1=g_new;
    file << endl << endl;
  }
  file.close();
  return g_new;
}
  

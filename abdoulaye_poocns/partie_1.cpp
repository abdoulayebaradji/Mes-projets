#include <iostream>
#include<cassert>
#include<cmath>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<sstream>
#include<complex>
#include"SEIR.h"
#include"SEIR_modified.h"
#include"SIR_KPP.h"
using namespace std;

int main()
{
  // choix du schema numerique
  int c;
  cout << "quel schema numera utilisé ? " << endl;
  cout << " 1 pour cranck nicolson" << endl;
  cout << " 2 pour heun" << endl;
  cout << " 3 pour RK4" << endl;
  cin >> c;
  // nombre total de population initial
  double  N= 1.0; double N1=60800000.00; double N_F=67e+6;double N_I=60e+6;double N_S=47e+6;
  
  // nombre total de jour de simulation T
  int T=30; int T1=90;

  // nombre de subdivison (points)
  int Nt=1000; int Nt1=900;

  // pas de l'intervalle
  double h=(float)T/(Nt-1); double h0=(float)T1/(Nt1-1);
  cout << "h= " << h << endl;
  cout << "h0= " << h0 << endl;

  // solution initiale
  
  vecteur L0(4); vecteur L11(4);vecteur L0_F(4);vecteur L0_I(4);vecteur L0_S(4);

  L0(0)=0.9; L0(1)=0.1; L0(2)=0; L0(3)=0;
  cout << "L0= " << L0 << endl;

  L11(0)=N1-0.0001*N1; L11(1)=0; L11(2)=0.0001*N1; L11(3)=0;
  cout << "solution intiale wuhan= " << L11 << endl;

  L0_F(0)=N_F-100;L0_F(1)=0;L0_F(2)=100;L0_F(3)=0;
  cout << "solution initiale france : " << L0_F << endl;

  L0_I(0)=N_I-132;L0_I(1)=0;L0_I(2)=132;L0_I(3)=0;
  cout << "solution initale Italy :" << L0_I << endl;

  L0_S(0)=N_S-136;L0_S(1)=0;L0_S(2)=136;L0_S(3)=0;
  cout << "solution initiale Spain: " << L0_S << endl;
  
  // definition des paramètres
  double beta=0.9; double gama=0.2; double sigma=0.5; double mu=0.0; double nu=0.0;

  // défintitions des paramètres pour la population de wuhan
  double beta1=1.2;double gama1=0.11;double sigma1=0.143; double mu1=0.03;double nu1=0; double beta1_news=0.25;

  // conditions d'arrêt d'itérations pour newton.
  double esp= 10e-10;

 
  SEIR f1,df1,fi1,dfi1,ncranck1,f,df,fi,dfi,ncranck,k_demi,k2,heun_seir,k_demi1,k2w,heun_seir_wuhan,RK4_seir,RK4_seir_wuhan,rk4_ref,RK4_seir_france,RK4_seir_italy,RK4_seir_spain,RK4_spain_isoleted,RK4_italy_isoleted,RK4_france_isoleted;
  
  //création du fichier
  ofstream file1("cranck_nicolson.txt",ios::out);
 
  // affichage des résultats lorsque la population totale vaut 1 :
  /*
   cout << "f= " <<  endl;
   cout << f.F(beta,sigma,gama,mu,nu,N,L0) << endl;
   cout << "df= " << endl;
   cout << df.DF(beta,sigma,gama,mu,nu,L0,N);
   cout << "fi= " << endl;
   cout << fi.phi(beta,sigma,gama,mu,nu,N,L0,h);
   cout << "dfi= " << endl;
   cout << dfi.dphi(beta,sigma,gama,mu,nu,N,L0, h);
  */
  
   cout << "solution cranck nicolson  " << endl;
   cout << ncranck.Newton_cranck_nicolson(beta,sigma,gama,mu,nu,Nt,h,N,L0,esp, "cranck_nicolson.txt");

  
  //affichage des résultats pour cranck nicolson population wuhan

  ofstream file2("cranck_nicolson_wuhan.txt", ios::out);
  /*
  cout << "f1= " <<  endl;
  cout << f1.F(beta1,sigma1,gama1,mu1,nu1,N1,L11) << endl;
  cout << "df1= " << endl;
  cout << df1.DF(beta1,sigma1,gama1,mu1,nu1,L11,N1);
  cout << "fi1= " << endl;
  cout << fi1.phi(beta1,sigma1,gama1,mu1,nu1,N1,L11,h0);
  cout << "dfi1= " << endl;
  cout << dfi1.dphi(beta1,sigma1,gama1,mu1,nu1,N1,L11, h0);
  */
  cout << "solution pour le schema cranck nicolson population de wuhan " << endl;
  cout << ncranck1.Newton_cranck_nicolson(beta1,sigma1,gama1,mu1,nu1,Nt1,h0,N1,L11,esp, "cranck_nicolson_wuhan.txt");
  
  // affichage des solutions du modèle SEIR pour le schéma de Heun

  // affichage pour nombre de population totale qui vaut 1 : 
  ofstream file3("heun.txt",ios::out);
  /*
  cout << "k_demi= " << endl;
  cout << k_demi.L_demi(beta,sigma,gama,mu,nu,N,h,L0);
  cout << "k2= " << endl;
  cout << k2.k_2(beta,sigma,gama,mu,nu,N,h,L0);
  */
  cout << "solution pour le schema heun : " << endl;
  cout << heun_seir.Heun(beta,sigma,gama,mu,nu,N,Nt,h,L0,"heun.txt");

  // affichage des solutions pour la population de wuhan :
  ofstream file4("heun_wuhan.txt",ios::out);
  /*
  cout << "k_demi1 :  " << endl;
  cout << k_demi1.L_demi(beta1,sigma1,gama1,mu1,nu1,N1,h0,L11);
  cout << "k2w : " << endl;
  cout << k2w.k_2(beta1,sigma1,gama1,mu1,nu1,N1,h0,L11);
  */
  cout << "solution pour le schema heun population de wuhan : " << endl;
  cout << heun_seir_wuhan.Heun(beta1,sigma1,gama1,mu1,nu1,N1,Nt1,h0,L11,"heun_wuhan.txt");
  
  // affichage des solutions du modèle SEIR pour le schema RK4 :

  // affichage pour nombre de population totale qui vaut 1 :
  ofstream file5("RK4.txt",ios::out);
  cout << "solution pour le schéma RK4 : " << endl;
  cout << RK4_seir.RK4(beta,sigma,gama,mu,nu,N,Nt,h,L0,"RK4.txt");

  // affichage des solutions pour la population de wuhan :
  
  ofstream file6("RK4_wuhan.txt",ios::out);
  cout << "solution pour le schéma RK4 population de wuhan : " << endl;
  cout << RK4_seir_wuhan.RK4(beta1,sigma1,gama1,mu1,nu1,N1,Nt1,h0,L11,"RK4_wuhan.txt");

  // affichage des solutions pour la population de france :
  ofstream file7("RK4_france.txt",ios::out);
  cout << "solution pour le schema RK4 population de france :" << endl;
  cout << RK4_seir_france.RK4(beta1,sigma1,gama1,mu1,nu1,N_F,Nt1,h0,L0_F,"RK4_france.txt");
  
  //cas avec l'isolement :
  ofstream file8("RK4_france_isolement.txt",ios::out);
  cout << "solution schema RK4 population de france isolement :" << endl;
  cout << RK4_france_isoleted.RK4(beta1_news,sigma1,gama1,mu1,nu1,N_F,Nt1,h0,L0_F,"RK4_france_isolement.txt");

  // affichage des solutions pour la population d'Italie :
  ofstream file9("RK4_italy.txt",ios::out);
  cout << "solution pour le schema RK4 population d'italy :" << endl;
  cout << RK4_seir_italy.RK4(beta1,sigma1,gama1,mu1,nu1,N_I,Nt1,h0,L0_I,"RK4_italy.txt");
  // cas avec l'isolement :
  ofstream file10("RK4_italy_isolement.txt",ios::out);
  cout << "solution schema RK4 population d'italy isolement :" << endl;
  cout << RK4_italy_isoleted.RK4(beta1_news,sigma1,gama1,mu1,nu1,N_I,Nt1,h0,L0_I,"RK4_italy_isolement.txt");

  // affichage des solutions pour la population de spain :
  ofstream file11("RK4_spain.txt",ios::out);
  cout << "solution pour le schema RK4 population de spain :" << endl;
  cout << RK4_seir_spain.RK4(beta1,sigma1,gama1,mu1,nu1,N_S,Nt1,h0,L0_S,"RK4_spain.txt");

  // cas avec l'isolement :
  ofstream file12("RK4_spain_isolement.txt",ios::out);
  cout << "solution pour le schema RK4 population spain isolement :" << endl;
  cout << RK4_spain_isoleted.RK4(beta1_news,sigma1,gama1,mu1,nu1,N_S,Nt1,h0,L0_S,"RK4_spain_isolement.txt");
  
  // solution de référence RK4 pour un pas de l'ordre de 0.001 ou 0.0001 :

  ofstream file13("RK4_ref.txt",ios::out);
  cout << "solution de reference pour le schéma RK4 : " << endl;
  vecteur S0(4);
  S0=  rk4_ref.RK4(beta,sigma,gama,mu,nu,N,30000,0.001,L0,"RK4_ref.txt");
  cout << "S0= " << S0;

  // solution de référence RK4 wuhan pour un pas de l'ordre de 0.001 ou 0.0001 :

  ofstream file14("RK4_wuhan_ref.txt",ios::out);
  cout << "solution de reference pour le schéma RK4 population wuhan : " << endl;
  vecteur S0_wuhan(4);
  S0_wuhan=  rk4_ref.RK4(beta1,sigma1,gama1,mu1,nu1,N1,90000,0.001,L11,"RK4_wuhan_ref.txt");
  cout << "S0_wuhan= " << S0_wuhan;
  
  //calcul de l'erreur et  solutions obtenues avec des h plus grands puis on le divise par 2  :
  
  vecteur err(4);
  double e_inf;
  vecteur S(4);
  int nt_1;
  double h_new;
  ofstream file15 ("solution_erreur.txt",ios::out);
  for (int i=0; i<10;i++){
    h_new=(float) 1/pow(2,i);
    nt_1= T1/h_new;
    switch (c){
      case 1:{
	S=ncranck1.Newton_cranck_nicolson(beta1,sigma1,gama1,mu1,nu1,nt_1,h_new,N1,L11,esp,"cranck_nicolson_wuhan.txt");
	
	break;
       }
	case 2:{
	  S=heun_seir_wuhan.Heun(beta1,sigma1,gama1,mu1,nu1,N1,nt_1,h_new,L11,"heun_wuhan.txt");
	  break;
	}
	  case 3:{
	    S=RK4_seir_wuhan.RK4(beta1,sigma1,gama1,mu1,nu1,N1,nt_1,h_new,L11,"RK4_wuhan.txt");
	    break;
	  }
	    default:
	      cout << " choix non prise en compte" << endl;
    }
    err= S - S0_wuhan;
    e_inf=err.norme_inf();
    cout << "e_inf= "  << e_inf <<  endl;
    file15 << e_inf << endl;
    file15 << endl << endl;
  }
  file15.close();
  
  /////////////////// Schema SEIR modifed ///////////////////////////////////

  // choix du schema numérique
  
  int ch;
  cout << "quel schema numera utilisé ? " << endl;
  cout << " 4 pour cranck modified" << endl;
  cout << " 5 pour heun modified" << endl;
  cout << " 6 pour RK4 modified" << endl;
  cin >> ch;

  // nombre total de population initial
  double  N00= 11081000;

  // nombre total de jour de simulation T
  int T00= 60;

  // nombre de subdivison (points)
  int Nt00=1000;

  // pas de l'intervalle
  double h00=(float)T00/Nt00;
  
  cout << "h00= " << h << endl;
  // solution initiale
  
  vecteur L00(8);
  L00(0)=11081000; L00(1)=105; L00(2)=28; L00(3)=54; L00(4)=739; L00(5)=1; L00(6)=1; L00(7)=2;
  cout << "L00= " << L00 << endl;

  // definition des paramètres
  double c0=14.781; double beta0=2.1011e-8; double q0=1.8887e-7; double sigma0=(float)1/7; double lambda0=(float)1/14; double rho0=0.86834; double delta_i0=0.13266; double delta_q0=0.1259; double gama_i0=0.33029; double gama_a0=0.13978; double gama_h0=0.11624; double alpha0=1.7826e-5; double theta0=2.5;

  SEIR_modified g0,dg0,fi0,dfi0,l_demi0,k_20,k2_rk0,k3_rk0,k4_rk0,rk40,heun0,newton_cranck0,rk40_ref;

  // test de nos fonction g,dg, phi, dg

  cout << "g0= " << endl;
  cout << g0.G(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,L00);
  cout << "dg0= " << endl;
  cout << dg0.DG(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,L00);
  cout << "fi0= " << endl;
  cout << fi0.phi_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,L00,h00);
  cout << "dfi0= " << endl;
  cout << dfi0.dphi_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,L00,h00);

  // affichage des solutions pour cranck modified
  
  ofstream file16("cranck_modified.txt",ios::out);
  cout << "solution pour le schema cranck nicolson modified " << endl;
  cout << newton_cranck0.Newton_cranck_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,Nt00,h00,L00,esp,"cranck_modified.txt");

  // affichage des solutions pour le schema heun modified
  
  ofstream file17("heun_modified.txt",ios::out);
  cout << "solution pour le schema heun modified " << endl;
  cout << heun0.Heun_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,Nt00,h00,L00,"heun_modified.txt");

  // affichage des solutions pour le schema RK4 modified

  ofstream file18("RK4_modified.txt",ios::out);

  cout << "solution pour le schema RK4 modified" << endl;
  cout << rk40.RK4_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,Nt00,h00,L00,"RK4_modified.txt");


  /////// calcul de l'erreur pour le modèle SEIR modified ////////////////////

  // solution de reference pour le modèle SEIR modified

  ofstream file19("Rk4_ref_modified.txt",ios::out);
  vecteur S0_modified(8);
  S0_modified=rk40_ref.RK4_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,60000,0.001,L00,"Rk4_ref_modified.txt");
  cout << "solution de reférence pour le schema SEIR modified" << endl;
  cout << S0_modified;

  // calcul de l'erreur et solution pour des pas plus grand :
  
  vecteur err_modified(8);
  double e_inf_modified;
  vecteur S_modified(8);
  int nt_00;
  double h_new00;
  ofstream file20 ("solution_erreur_modified.txt",ios::out);
  for (int i=0; i<10;i++){
    h_new00=(float) 1/pow(2,i);
    nt_00= T00/h_new00;
    switch (ch){
      case 4:{
	S_modified= newton_cranck0.Newton_cranck_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,nt_00,h_new00,L00,esp,"cranck_modified.txt");
	
	break;
       }
	case 5:{
	  S_modified=heun0.Heun_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,nt_00,h_new00,L00,"heun_modified.txt");
	  break;
	}
	  case 6:{
	    S_modified= rk40.RK4_modified(c0,beta0,q0,sigma0,lambda0,rho0,delta_i0,delta_q0,gama_i0,gama_a0,gama_h0,alpha0,theta0,nt_00,h_new00,L00,"rk4_modified.txt");
	    break;
	  }
	    default:
	      cout << " choix non prise en compte" << endl;
    }
    err_modified= S_modified - S0_modified;
    e_inf_modified=err_modified.norme_inf();
    cout << "e_inf_modified= "  << e_inf_modified <<  endl;
    file20 << e_inf_modified << endl;
    file20 << endl << endl;
  }
  file20.close();
  
  
   ////////////////// Partie SIR_FISHER_KPP //////////////////////////////

  // définitions des paramètres :

  int nx=10; int ny=10;double dx=(float)1/(nx+1);double dy=(float)1/(ny+1);int L=1;double r=(float) L/2;double T_kpp=30; double D=0.001; double h_kpp=0.005; double k_si=0.25; double k_ir=0.1;

  // test pour la matrice du laplacien en 2d
  
  SIR_KPP LA,TA,a,b,rk4_kpp,g_kpp,g1,g2,g3,phi_kpp,dphi_kpp,ei_kpp,dg_kpp,test0;
  /*
  cout << a.laplacien0(4,4);
  cout << b.iden(4);
  */
  cout << "laplacien en 2d : " << endl;
  
  cout << LA.laplacien_2d(nx,ny);
  
  cout << "A : " << endl;
  cout << TA.A(nx,ny);
  
  
  // création du vecteur initiale :
 
  vecteur U0(3*nx*ny);
  for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
      if (-r<=  (i+1)*dx &&  r>= (i+1)*dx && -r <= (j+1)*dy &&  r >= (j+1)*dy){
	U0(i*(nx)+j)=0.9;
        U0(pow(nx,2)+i*(nx)+j)=0.1;
      }else{
	U0(i*(nx)+j)=1;
        U0(pow(nx,2)+i*(nx)+j)=0;
      }
    }
  }
  cout << "U0 : " << endl;
  cout << U0;
  
  ofstream file21("solution_initiale_kpp.txt",ios::out);
  for(int i=0;i<3*nx*ny;i++){
    file21 << U0(i);
    file21 << endl << endl;
  }
  file21.close();
  
  // test des fonctions g1,g2,g3 et g :
  
  cout << "g1= " << endl;
  cout << g1.G_1(nx,ny,k_si,U0);
  
  cout << "g2= " << endl;
  cout << g2.G_2(nx,ny,k_si,k_ir,U0);
  
  cout << "g3= " << endl;
  cout << g3.G_3(nx,ny,k_ir,U0);
  
  cout << "g= " << endl;
  cout << g_kpp.G_KPP(nx,ny,k_si,k_ir,U0);
  
  // affichage de la solution de l'equation 4 :
  
  ofstream file22("rk4_fisher_kpp.txt",ios::out);
  cout << "solution du schema RK4_kPP pour la resolution l'equation 4" << endl;
  cout << rk4_kpp.RK4_KPP(nx,ny,T_kpp,D,h_kpp,k_si,k_ir,U0,"rk4_fisher_kpp.txt"); 
  
  /////////////// EULER IMPLICITES : //////////////////////////////
  
  cout << "phi : " << endl;
  cout << phi_kpp.PHI(nx,ny,h_kpp,D,k_si,k_ir,U0);

  cout << "dg : " << endl;
  cout << dg_kpp.DG(nx,ny,k_si,k_ir,U0);

  cout << "dphi : " << endl;
  cout << dphi_kpp.DPHI(nx,ny,D,h_kpp,k_si,k_ir,U0);
  ofstream file23("ei_fisher_kpp.txt",ios::out);
  cout << "schema euler implite : " << endl;
  cout << ei_kpp.EI_KPP(nx,ny,T_kpp,k_si,k_ir,h_kpp,D,U0,esp,"ei_fisher_kpp.txt");
  
}

  

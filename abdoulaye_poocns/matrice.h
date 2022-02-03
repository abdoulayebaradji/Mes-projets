/**********************************************
 *  TP 2 : Creation de la classe matrice
 *  Caterina Calgaro, POOCN 2019-2020
 *  Compilateur g++
 ***********************************************/

#if !defined (__IOSTREAM_H)
#include <iostream>
#endif
#if !defined (__FSTREAM_H)
#include <fstream>
#endif
#if !defined (__CASSERT_H)
#include <cassert>
#endif
#if !defined(__CMATH_H)
#include <cmath>
#endif
#include "vecteur.h"

using namespace std;

class matrice : public vecteur
// class matrice                // c'est la meme chose
{
private:
    int      size1;  // dimension 1 de la matrice (nb de lignes)
    int      size2;  // dimension 2 de la matrice (nb de colonnes)
    vecteur* matrix; // chaque ligne est un vecteur
  
public:
    matrice() {this->size1 = this->size2 = 0;};  // constructeur par defaut
    matrice(int,int);                  // constructeur en donnant les 2 dimensions
    matrice(const matrice&);           // constructeur par copie
    ~matrice();                        // destructeur

    int dim1() const {return this->size1;}    // retourne la 1ere dimension de la matrice
    int dim2() const {return this->size2;}    // retourne la 2eme dimension de la matrice
    
    vecteur& operator [] (int) const;         // retourne une ligne
    double&  operator () (int,int) const;     // retourne un coefficient
    matrice& operator = (const matrice&);     // surcharge de = (affectation)
	
    matrice transpose(const matrice&);    // transposee d'une matrice
    vecteur colonne(const int);           // retourne une colonne de la matrice

    matrice operator + (const matrice&); // somme de 2 matrices
    matrice operator - (const matrice&); // difference de 2 matrices
    matrice operator * (const matrice&); // produit de 2 matrices
  
    matrice operator * (const double&); // matrice x cte
    friend matrice operator * (const double&, const matrice&); // cte x matrice
 
    vecteur operator * (const vecteur&); // produit A*v
    friend vecteur operator * (vecteur&, const matrice&); // produit v^t * A
	
    friend ostream& operator << (ostream&, const matrice &);
    friend istream& operator >> (istream&, const matrice &);
    
    friend void save_matr(const char*, const matrice &);   // ecriture dans un fichier
    friend matrice read_matr(const char *);                // lecture dans un fichier

    matrice lu(vecteur&);                                // factorisation lu
    vecteur solvelu(vecteur&, vecteur&);                 // resolution lu
    double norme_inf(); // norme infinie d'une matrice
};

//------------------------------------------------------
//------------------------------------------------------
// Constructeur v1
matrice::matrice(int n, int m)
{
  assert((n>0)&&(m>0));
  this->size1 = n;
  this->size2 = m;
  this->matrix = new vecteur[n];
  for (int i=0; i<n ; i++)
    this->matrix[i]=vecteur(m);
}
//-----------------------------------------------------
// Constructeur par copie
matrice::matrice(const matrice & mat)
{
  assert((mat.size1)&&(mat.size2));
  this->size1=mat.size1;
  this->size2=mat.size2;
  this->matrix=new vecteur[this->size1];
  for (int i=0; i<mat.size1 ; i++)
    this->matrix[i]=mat.matrix[i];
}
//-----------------------------------------------------
// Destructeur
matrice::~matrice()
{
  if (this->size1||this->size2) delete[] this->matrix;
  this->size1=0;
  this->size2=0;
}
//-------------------------------------------------------
// Retourne la ligne i
vecteur& matrice::operator [] (int i) const
{
  assert((i>=0)&&(i<this->size1));
  return this->matrix[i];
}
//------------------------------------------------------
// Retourne le coefficient (i,j)
double& matrice::operator () (int i, int j) const
{
  assert((i>=0)&&(i<this->size1)&&(j>=0)&&(j<this->size2));
  return this->matrix[i](j);
}
//-------------------------------------------------------
// Surcharge = (affectation)
matrice& matrice::operator = (const matrice& mat)
{
  assert((mat.size1>0)&&(mat.size2>0));
  if ((!this->size1)||(!this->size2))
    {this->size1=mat.size1;
     this->size2=mat.size2;
     this->matrix=new vecteur[this->size1];}
  assert((this->size1==mat.size1)&&(this->size2==mat.size2));
  for (int i=0; i<this->size1 ; i++)
    this->matrix[i]=mat.matrix[i];
  return (*this);
  }
//-------------------------------------------------------
// Transposee d'une matrice
matrice matrice::transpose (const matrice& A)
{
	assert((A.size1>0)&&(A.size2>0));
	matrice C(A.size2,A.size1);
	for (int i=0; i<A.size2; i++)
		for (int j=0; j<A.size1; j++)
			C(i,j) = A(j,i);
	return C;
}
//-------------------------------------------------------
// Retourner la colonne j d'une matrice
vecteur matrice::colonne (const int j)
{
    assert((j>=0)&&(j<this->size2));
    vecteur vcol(this->size1);
    for (int i=0; i<this->size1; i++)
        vcol(i) = this->matrix[i](j);
    // vcol(i) = this->matrix[i][j]; // ici probleme de compilation
    return vcol;
}
//-------------------------------------------------------
// Surcharge operateur + : somme de 2 matrices
matrice matrice::operator + (const matrice& A)
{
	assert((this->size1>0)&&(this->size1==A.size1));
	assert((this->size2>0)&&(this->size2==A.size2));
	matrice C(this->size1,this->size2);
	for (int i=0; i<this->size1; i++)
		C[i] = this->matrix[i]+A[i];
	return C;
}
//-------------------------------------------------------
// Surcharge operateur - : difference de 2 matrices
matrice matrice::operator - (const matrice& A)
{
	assert((this->size1>0)&&(this->size1==A.size1));
	assert((this->size2>0)&&(this->size2==size2));
	matrice C(this->size1,this->size2);
	for (int i=0; i<this->size1; i++)
		C[i] = this->matrix[i]-A[i];
	return C;
}
//-------------------------------------------------------
// Surcharge operateur * : produit de 2 matrices
matrice matrice::operator * (const matrice& A)
{
    assert((this->size1>0)&&(this->size2>0)&&(A.size1>0)&&(A.size2>0));
	assert(this->size2==A.size1);
	matrice TA(A.size2,A.size1);
	matrice C(this->size1,A.size2);
	for (int i=0; i<this->size1; i++)
	     for (int j=0; j<A.size2; j++)
		  C(i,j) = this->matrix[i]*TA.transpose(A)[j];
	return C;
}
//-------------------------------------------------------
// Surcharge operateur * : matrice x cte
matrice matrice::operator * (const double& cte)
{
    assert((this->size1>0)&&(this->size2>0));
    matrice C(this->size1,this->size2);
    for (int i=0; i<this->size1; i++)
         C[i] = this->matrix[i]*cte;
    return C;
}
//-------------------------------------------------------
// Surcharge operateur * : cte x matrice
matrice operator * (const double& cte, const matrice& A)
{
	assert((A.size1>0)&&(A.size2>0));
	matrice C(A.size1,A.size2);
	for (int i=0; i<A.size1; i++)
	     C[i] = cte*A[i];
	return C;
}
//-------------------------------------------------------
// Surcharge operateur * : produit A*v
vecteur matrice::operator * (const vecteur& v)
{
    assert((this->size1>0)&&(this->size2>0)&&(this->size2==v.dim()));
    vecteur w(this->size1);
    for (int i=0; i<this->size1; i++)
        w(i) = this->matrix[i]*v;
    return w;
}
//-------------------------------------------------------
// Surcharge operateur * : produit v^t * A
vecteur operator * (vecteur& v, const matrice& A)
{
	assert((A.size1>0)&&(A.size2>0)&&(A.size1==v.dim()));
	vecteur w(A.size2);
	matrice C(A.size2,A.size1);
	for (int i=0; i<A.size2; i++)
	     w(i) = v*C.transpose(A)[i];
	return w;
}
//-------------------------------------------------------
// Surcharge operator <<
ostream& operator << (ostream& s, const matrice & mat)
{
  assert((mat.size1>0)&&(mat.size2>0));
  for (int i=0; i<mat.size1; i++)
    s << mat[i];
  s << "\n";
  return s; 
}
//-------------------------------------------------------
// Surcharge operator >>
istream& operator >> (istream& s, const matrice & mat)
{
  assert((mat.size1>0)&&(mat.size2>0));
  cout << "Entrer au clavier une matrice avec " << mat.size1 <<" lignes et " << mat.size2 << " colonnes \n";
  for (int i=0; i<mat.size1; i++)
    s >> mat[i];
  return s; 
}
//---------------------------------------------
// Sauvegarde d'une matrice dans un fichier
void save_matr(const char *Nomfich, const matrice & mat)
{
    ofstream fichier;
    fichier.open(Nomfich); // ouverture du fichier
    
    assert(mat.size1>0 && mat.size2>0);
    fichier << "Dimension de la matrice = " << mat.dim1() << " " << mat.dim2() << "\n";
    for (int i=0; i<mat.size1 ; i++)
        fichier << "i= " << i << "\t" << mat[i];
    fichier.close();           // fermature du fichier
}
//---------------------------------------------
/*
// lecture d'une matrice dans un fichier
// v1 : il faut avoir en memoire une matrice mat de taille suffisamment grande, pas optimal
 void read_matr(const char *Nomfich, matrice & mat)
 {
 ifstream fichier;
 fichier.open(Nomfich); // ouverture du fichier
 
 if (!fichier)
 { cerr << "Probleme d'ouverture de fichier " << endl;
 exit(1);
 }
 
 if (!fichier.eof()) 
 fichier >> mat.dim1() >> mat.dim2();
     
 int i=0;    
 while (!fichier.eof())
 {
 fichier >> mat[i];
 i+=1;
 }
 
 fichier.close();  // fermature du fichier
 cout << "Dimension de la matrice = " << mat.dim1() << " " << mat.dim2() << "\n";
 }
*/

//---------------------------------------------
// lecture d'une matrice dans un fichier
// v2 : on va obtenir une matrice de la bonne taille
matrice read_matr(const char *Nomfich)
{
    ifstream fichier;
    fichier.open(Nomfich); // ouverture du fichier
    
    if (!fichier)
    { cerr << "Probleme d'ouverture de fichier " << endl;           
        exit(1);
    }
    
    int m,n;
    if (!fichier.eof()) 
      fichier >> m >> n;
    matrice M(m,n);
    // lecture d'un vecteur dans un fichier
    for (int i=0;i<M.size1;i++)
    {
        for (int j=0;j<M.size2;j++)
        {
            fichier >> M(i,j);
            if (fichier.eof())
                cout << " Erreur dans le fichier, matrice incomplete \n ";
        }
    }
    
    fichier.close(); // fermeture du fichier
    return M;
}

//--------------------------------------
// factorisation LU d'une matrice carree
// N.B. ici Piv est un vecteur de double, converti ensuite en entier. Sinon utiliser le template vecteur

matrice matrice::lu(vecteur& Piv)
{
    assert(this->size1>0 && this->size2>0 && this->size1==Piv.dim());
    matrice C(this->size1,this->size2);
    if (this->size1!=this->size2)
    {cout << " Matrice non carree, factorisation LU impossible \n"; return C;}
    for (int i=0; i<this->size1;i++)
    {
        Piv(i)=i;
        C.matrix[i]=this->matrix[i];
    };
    
    for (int k=0; k<this->size2-1;k++)
    {
        
        double max = abs(C(int(Piv(k)),k));				  // Factorisation LU	
        int lmax = k;
        
        for (int l=k+1; l<this->size1;l++)
        {
            if (max<abs(C(int(Piv(l)),k)))
            {
                max = abs(C(int(Piv(l)),k));
                lmax = l;
            }     
        };
        
        int aux = Piv(k);
        Piv(k) = Piv(lmax);
        Piv(lmax) = aux;
        
        for (int p=k+1; p<this->size1;p++)
        {
            C(int(Piv(p)),k)=C(int(Piv(p)),k) / C(int(Piv(k)),k);
            
            for (int i=k+1; i<this->size1; i++)
                C(int(Piv(p)),i)=C(int(Piv(p)),i) - C(int(Piv(p)),k) * C(int(Piv(k)),i);
            
        };
    };
    return C;
}               

//------------------------------------------------------
// resolution systemes lineaires triangulaire inf et sup 
// N.B. ici Piv est un vecteur de double, converti ensuite en entier. Sinon utiliser le template vecteur

vecteur matrice::solvelu(vecteur& b, vecteur& Piv)
{
    assert(this->size1>0 && this->size2>0 && this->size1==Piv.dim() && this->size1==b.dim());
    vecteur sol(this->size2);
    int p;
    for (int i=0; i<this->size1; i++)
        sol(i) = b(int(Piv(i)));
    
    for (int k=0; k<this->size1; k++)
    {
        for (int j=0; j<k; j++)
        {p = int(Piv(k)); 
            sol(k) = sol(k) - this->matrix[p](j) * sol(j);}     		  // Résolution LU
    };
    
    for (int h=this->size1-1; h>=0;h--)
    {
        p = int(Piv(h));
        for(int j=h+1;j<this->size2;j++)
            sol(h) = sol(h) - this->matrix[p](j) * sol(j);
        sol(h) = sol(h) / this->matrix[p](h);
    };
    return sol;
} 

double matrice::norme_inf()
{
  double norme, norme_max;
  norme_max = 0;
  for (int i = 0; i < this->size2; i++){
      norme = this->matrix[i].norme_1();
      if (norme > norme_max)
	norme_max = norme;
  }
  return norme_max;
}

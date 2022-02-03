 /**********************************************
*  TP 1 : Creation de la classe vecteur
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
#if !defined(__CSTDLIB_H)
#include <cstdlib>
#endif

using namespace std;

class vecteur
{
private:
    int size;
    double* p;

public:
    vecteur() {size=0;};     // constructeur par default
    vecteur(int);            // constructeur
    vecteur(const vecteur&); // constructeur par copie
    ~vecteur();              // destructeur

    int dim() const {return this->size;}                       // retourne la dimension du vecteur
    double& operator () (int);                           // retourne l'element d'adresse int
    vecteur& operator = (const vecteur&);                      // surcharge de = (affectation) : correcte
    //vecteur operator = (const vecteur&);                     // surcharge de = (affectation) : fausse
    
    vecteur operator / (const double&);                        // surcharge de /
    vecteur operator * (const double&);                        // surcharge de * : res=v*cte
    friend vecteur operator * (const double&, const vecteur&);     // surcharge de * : res=cte*v

    double  operator * (const vecteur&);    // surcharge de * : prod=v1*v2 (produit scalaire asymetrique)
    vecteur operator + (const vecteur&);    // surcharge de +  : asymetrique
    vecteur operator - (const vecteur&);    // surcharge de -  : asymetrique
    int     operator == (const vecteur&);   // surcharge de == : asymetrique

    double norme_inf();                     // norme infinie d'un vecteur
    double norme_1  ();                     // norme 1 d'un vecteur
    double norme_2  ();                     // norme euclidienne d'un vecteur
    double norme_l2 ();                     // norme l2 d'un vecteur
    friend ostream& operator << (ostream&, const vecteur&);       // sortie a l'ecran
    friend istream& operator >> (istream&, const vecteur&);       // lecture au clavier

    friend void save_vect (const char*, const vecteur&);          // ecriture dans un fichier
    friend int  read_vect (const char*, vecteur&);                // lecture dans un fichier
};

//---------------------------------------------
vecteur::vecteur(int s)
{
  assert(s>0);
  this->size=s;
  this->p= new double[s];
  for(int i=0; i<s; i++)
    this->p[i]=0.0;
  //cout << "constructeur avec taille s " << s << " " << this << "\n";
}
//---------------------------------------------
vecteur::vecteur(const vecteur& v)
{
  assert(v.size>0);
  this->size=v.size;
  this->p= new double[v.size];
  for(int i=0; i<size; i++)
    this->p[i]=v.p[i];
  //cout << "constructeur par copie " << this << "\n";
}
//---------------------------------------------
vecteur::~vecteur()
{
  if (this->size !=0) delete this->p;
  this->size=0;
  //cout << "destructeur " << this << "\n";
}
//---------------------------------------------
// Surcharge operator ()
double& vecteur::operator () (int i)
{
  assert((i>=0)&&(i<this->size));
  return this->p[i];
}
//---------------------------------------------
// Surcharge affectation : CORRECTE
vecteur& vecteur::operator = (const vecteur& v)
{
    assert(v.size);
//  if (&v != this) // je peux laisser le test if, mais alors j'ai un warning... mais ca marche aussi
//  {
    assert((this->size==v.size) || (this->size==0));
    if (!size) {
        this->size=v.size;
        this->p=new double[this->size];
    };
    for(int i=0; i<v.size; i++)
        p[i]=v.p[i];
    //cout << "surcharge affectation " << this << "\n";
    return *this;
//  }
}
//---------------------------------------------
/*
// Surcharge affectation : FAUSSE
vecteur vecteur::operator = (const vecteur& v)
{
	assert(v.size);
	vecteur w(v.size); // constructeur avec taille donnee
	for(int i=0; i<v.size; i++)
		w.p[i]=v.p[i];
		cout << "surcharge affectation fausse \n";
		return w;
}
 */
//---------------------------------------------
// Surcharge operator /
vecteur vecteur::operator / (const double& a)
{
  vecteur w(this->size);
  assert(a !=0.0);
  for (int i=0; i<this->size; i++)
    w.p[i]=this->p[i]/a;
  //cout << "surcharge operateur / \n";
  return w;
}
//---------------------------------------------
// Surcharge operator * : res=v*cte
vecteur vecteur::operator * (const double& a)
{
  vecteur w(this->size);
  for (int i=0; i<this->size; i++)
    w.p[i]=this->p[i]*a;
  //cout << "surcharge operateur * \n";
  return w;
}
//---------------------------------------------
// Surcharge operator * : res=cte*v
vecteur operator * (const double& a, const vecteur &v)
{
  vecteur w(v.size);
  for (int i=0; i<v.size; i++)
    w.p[i]=a*v.p[i];
  //cout << "surcharge operateur * v2 \n";
  return w;
}
//---------------------------------------------
// // Surcharge operator * (produit scalaire asymetrique)
double vecteur::operator * (const vecteur &v)
{
	double res=0.0;
	assert((this->size>0) && (this->size==v.size));
	for (int i=0; i<v.size; i++)
		res+=this->p[i]*v.p[i];
    //cout << "surcharge operateur * (produit scalaire) \n";
    return res;
}

//---------------------------------------------
// surcharge de + : asymetrique
vecteur vecteur::operator + (const vecteur &b)
{
  assert((this->size>0) && (this->size==b.size));
  vecteur c(this->size);
  for (int i=0; i<this->size; i++)
    c.p[i]=this->p[i]+b.p[i];
  //cout << "surcharge operateur + \n";
  return c;
}
//---------------------------------------------
// surcharge de - : asymetrique
vecteur vecteur::operator - (const vecteur &b)
{
  assert((this->size>0) && (this->size==b.size));
  vecteur c(this->size);
  for (int i=0; i<this->size; i++)
    c.p[i]=this->p[i]-b.p[i];
  //cout << "surcharge operateur - \n";
  return c;
}
//---------------------------------------------
// surcharge de == : asymetrique
int vecteur::operator == (const vecteur &b)
{
  double epsilon = 1.e-15;
  assert(this->size>0);
  int vrai=(this->size==b.size);
  assert(vrai);
  for (int i=0; i<this->size; i++)
    vrai=pow((this->p[i]-b.p[i]),2)<=epsilon;
  //cout << "surcharge operateur == \n";
  return vrai;
}
//---------------------------------------------
  double vecteur::norme_inf()
{
  double nrm=-1.0;
  assert(this->size>0);
  for (int i=0; i<this->size; i++)
      if (abs(this->p[i])>nrm) nrm = abs(this->p[i]);
  return nrm;
}
//---------------------------------------------
  double vecteur::norme_1()
{
  double nrm=0.0;
  assert(this->size>0);
  for (int i=0; i<this->size; i++)
    nrm += abs(this->p[i]);
  return nrm;
}
//---------------------------------------------
  double vecteur::norme_2()
{
  double nrm=(*this)*(*this);  // surcharge de * (produit scalaire)
  return sqrt(nrm);
}
//---------------------------------------------
double vecteur::norme_l2()
{
    double nrm=(*this)*(*this);  // surcharge de * (produit scalaire)
    return (nrm/this->size);
}
//---------------------------------------------
// Surcharge operator <<
ostream& operator << (ostream& s, const vecteur& v)
{
  assert(v.size>0);
  for (int i=0; i<v.size ; i++)
    s<<v.p[i]<<" ";
  s<<"\n";
  return s;
}
//---------------------------------------------
// Surcharge operator >>
istream& operator >> (istream& s, const vecteur& v)
{
  assert(v.size>0);
  cout << "Entrer au clavier " << v.size << " valeurs \n";
  for (int i=0; i<v.size ; i++)
    s>>v.p[i];
  return s;
}
//---------------------------------------------
// Sauvegarde d'un vecteur dans un fichier
void save_vect(const char *Nomfich, const vecteur& v)
{
    ofstream fichier;
    fichier.open(Nomfich, ios::out); // ouverture du fichier

    assert(v.dim() != 0); // pas assert(v.size) car c'est un champ privé et save_vect est amie de la classe
    fichier << "Dimension du vecteur = " << v.size << endl;
    for (int i=0; i<v.dim() ; i++)
        fichier << "i= " << i << "\t" << v.p[i] << endl;
    fichier.close();           // fermature du fichier
}
//---------------------------------------------
// lecture d'un vecteur dans un fichier
int read_vect (const char *Nomfich, vecteur& v)
{
    ifstream fichier;
    fichier.open(Nomfich, ios::in); // ouverture du fichier

    assert(!fichier.fail()); // verifier ouverture correcte du fichier
    
    int i=0;
    while (!fichier.eof())
    {
        fichier >> v.p[i];
        i++;
        if (i>v.size)
        { cerr << "Probleme taille vecteur pas suffisante pour lire tout le fichier " << endl;
          exit(1);
        }
    }
  v.size=i-1;
    
  fichier.close();           // fermature du fichier
  return v.size;
}

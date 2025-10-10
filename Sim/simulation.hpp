#include "xlife++-libs.h"
#include "domain.hpp"
#include <filesystem>
#include <iostream>
#include <functional>
#include <fstream>
using namespace xlifepp;
using namespace std;
using namespace std::filesystem;
// ============= var ==============
extern Number ord;              // ordre de la methode utilisee

extern Real tf, length;         // dimension spatialle et temporelle du probleme
extern Real CFL;                // coefficient de cfl
extern Real dspace;             // pas d'espace du probleme (dx)

extern Real rho_0,mu_0 ;        // coeffs de l'espace Q- 
extern Real rho_1,mu_1 ;        // coeffs de l'espace Q+
extern Real cm ,cp ;            // vitesse en -/+
extern Real sm ,sp ;            // permeabilité en -/+
extern Real c_max;              // vitesse max (- et +)
extern String mod;              // mode choisi (flat, jump, interface, periode, energie)
extern String cond;             // condition au bord

                                // paramètres nécessaire aux fonctions  
extern Real vit;                // vitesse de l'interface
extern Real vit_abs;            // vitesse de l'interface en val absolue
extern Real epsi;               // paramètre de régularisation de l'interface
extern const Real wave_length;  // periode temporelle du mode periodique
extern const Real space_length; // periode spacial du mode periodique
extern Point k;                 // point fonctionnant comme vecteur de base
extern Real T;                  // periode temporelle du mode energie
extern Real Start_int;          // endroit ou commence l'interface 
extern Real Start_domaine;
extern Real Start_sol;          // endroit ou commence la solution

// definit la manière d'obtenir la valeur de la derivee en t=0
// l'un des deux doit etre vrai, la valeur est privilegier a la derivee
extern bool possede_derive, possede_valeur;

extern bool uni_dir ;           // definit si l'onde se propage dans une seul direction ou non
extern bool error_verb , error_calc;    // a t on des messages d'erreur

extern Real t, dtime, trho;     // paramètres temporelles

extern String fes;              // methode de quadrature
extern FESubType fesub;         // methode de quadrature xlife

extern Number savefile;         // frequence des images sur la simulation
extern Number savepic;          // nombres d'images sur la simulation
extern Number vbl;              // xlife parle t il ?

extern String choix;            // mu et rho en fonction de sigma, c ou l'inverse ? 
extern Number test_nombre;      // nombre de test dans les batteries de test

extern Real R_ref ,T_trans;             // coefficients de Reflexion-Transmission
extern Real tau_1, tau_2, tau_3, tau_4; // coefficients de dilatation/ contraction

extern Number counter;          // compte ce que l'on souhaite, permet de tester le code
extern Number error_count;      // compte des erreurs

// ================= fonctions ===================


void move(const string source, const string destination);
Real give_dt(const Real dspace, const Real max_c);     // rend dt en fonction de max_c, dspace (CFL et ordre sont implicite)

// conditions initiales
Real u_0(const Point&P, Parameters& pa = defaultParameters);    // u(x,0)
Real u_d(const Point&P, Parameters& pa = defaultParameters);    // u'(x,0) ou u(x,dt)
Real u_1(const Point&P, Parameters& pa = defaultParameters);    // u(x,dt)

// parametres physiques
Real rho(const Point& P, Parameters& pa = defaultParameters);   // rho(x,t)
Real mu(const Point& P, Parameters& pa = defaultParameters);    // mu(x,t)
Real beta_0(const Point& P, Parameters& pa = defaultParameters);// beta(x,t)    !!! beta deja pris comme noms de fonction

// solution exacte du probleme d'interface mobile
Real Fm (const Point&P, Parameters& pa = defaultParameters);    
Real Gm (const Point&P, Parameters& pa = defaultParameters);
Real Fp (const Point&P, Parameters& pa = defaultParameters);
Real Gp (const Point&P, Parameters& pa = defaultParameters);
Real u_ex(const Point& P, Parameters& pa = defaultParameters);

// interface / domaines au temps t
Real interface_func(const Point&P, Parameters& pa = defaultParameters);

// 2nd membre / terme source
Real g_1D(const Point& P, Parameters& pa = defaultParameters);
Real h(const Real& t, Parameters& pa = defaultParameters);
Real f(const Point& P, Parameters& pa = defaultParameters);

// fonctions raccords / splin
Real splin_exp(const Real x);
Real splin_lin(const Real x);
Real raccord (const Real x, const Real p1,const  Real p2,const Real a, const Real b);

// lance les simulation du probleme
void simulation(path dossier_trail);
void test_standard(path dossier, const vector<Real> S_s, const vector<Real> S_d, const vector<Real> V);

// definit les domaine du probleme d'energie croissante
Real f_2_temp(const Real x);        // par construction, cette fonction est 2T periodique
Real f_2(const Real x);             // passe d'une fonction definit sur [0,2T] a defini sur R+
Real f_1(const Real x);


class Domain_S
{
    public :
    Real c,sigma,rho,mu,epsi;                               // valeurs physiques dans le domaine
    vector<function<Real(Real,Real)>> bords;                // a valeurs dans {0,1}, si regularisation : [0,1] 
    Domain_S(Real r1, Real r2,String choix,Real _epsi =0)   // declaration d'un domaine
    {
        if(choix =="vitesse"){  c= r1;sigma= r2; rho=sigma/c; mu=sigma*c;}
        else {rho=r1;mu=r2;  c=sqrt(r2/r1);sigma=sqrt(r1*r2);}
        epsi=_epsi;
    }

// ajoute des bords au domaines, 
// si pas de bords, pas de domaines.

    // interface, independant du temps
    void ajout_bords_lin(Real _time, Real _side)    
    {
        function<Real(Real,Real)> bord([_time, _side](Real _x,Real _t){
            Real xi = _t -_time; xi *=_side;
            if(xi<-epsi){return 1.;}
            if(xi<epsi){return raccord(xi,1.,0.,-epsi,epsi);}
            else{return 0.;}
        });
        bords.push_back(bord);
    }

    // interface mobile a vitesse constante
    void ajout_bords_lin(Real _Start, Real _vitesse, Real _side)    
    {
        function<Real(Real,Real)> bord([_Start, _vitesse, _side](Real x,Real t){
            Real xi= x-_vitesse*t -_Start; xi*=_side;
            if(xi<-epsi){return 1.;}
            if(xi<epsi){return raccord(xi,1.,0.,-epsi,epsi);}
            else{return 0.;}
        }); 
        bords.push_back(bord);
    }

    // interface quelconque
    void ajout_bords(function<Real(Real,Real)> bordure)
    {
        bords.push_back(bordure);
    }

    // renvoie une valeur de rho en fonction de si x se trouve dans le domaine
    // si is_inside =1, rho_S = rho, si is_inside =0, rho_S =0, si is_inside in [0,1], rho > rho_S > 0;
    Real rho_S(Real x, Real t)
    {
        Real is_inside=1;
        for(int i=0; i<bords.size(); i++){ is_inside = is_inside*bords[i](x,t); }
        // std::cout<<"x= "<<x<<" "<<is_inside<<" "<<rho<<'\n';
        return rho*is_inside;
    }
    
    // meme dynamique pour mu
    Real mu_S(Real x, Real t)
    {
        Real is_inside=1;
        for(int i=0; i<bords.size(); i++){ is_inside = is_inside*bords[i](x,t); }
        // std::cout<<"x= "<<x<<" "<<is_inside<<" "<<mu<<'\n';
        return mu*is_inside;
    }

};

extern vector<Domain_S> Domaines;
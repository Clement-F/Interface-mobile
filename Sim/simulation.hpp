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
extern Number counter;
extern Number error_count;
extern Number ord;

extern Real tf, length;
extern Real CFL;         // coefficient de cfl
extern Real dspace;

extern Real rho_0,mu_0 ;       // coeffs de l'espace Q- 
extern Real rho_1,mu_1 ;       // coeffs de l'espace Q+
extern Real cm ,cp ;            // vitesse en -/+
extern Real sm ,sp ;            // permeabilité en -/+
extern Real c_max;            // vitesse max (- et +)
extern String mod;            // mode choisi (flat, jump, interface, periode, energie)
extern String cond;           // condition au bord

                              // paramètres nécessaire aux fonctions  
extern Real vit;              // vitesse de l'interface
extern Real vit_abs;          // vitesse de l'interface en val absolue
extern Real epsi;             // paramètre de régularisation de l'interface
extern const Real wave_length;   // periode temporelle du mode periodique
extern const Real space_length;  // periode spacial du mode periodique
extern Point k;   // point fonctionnant comme vecteur de base
extern Real T;                // periode temporelle du mode energie
extern Real Start_int;        // endroit ou commence l'interface 
extern Real Start_domaine;
extern Real Start_sol;        // endroit ou commence la solution

// definit la manière d'obtenir la valeur de la derivee en t=0
// l'un des deux doit etre vrai, la valeur est privilegier a la derivee
extern bool possede_derive, possede_valeur;
extern bool uni_dir ;
extern bool error_verb , error_calc;

extern Real t, dtime, trho;   // paramètres temporelles

extern String fes;
extern FESubType fesub;
extern Number savefile;
extern Number savepic;
extern Number vbl;

extern String choix;
extern bool pre_Start;
extern Number test_nombre;

extern Real R_ref ,T_trans;
extern Real tau_1, tau_2, tau_3, tau_4;


// ================= fonctions ===================


void move(const string source, const string destination);
Real give_dt(const Real dspace, const Real max_c);     // rend dt en fonction de max_c, dspace (CFL et ordre sont implicite)


Real u_0(const Point&P, Parameters& pa = defaultParameters);
Real u_d(const Point&P, Parameters& pa = defaultParameters);
Real u_1(const Point&P, Parameters& pa = defaultParameters);

Real rho(const Point& P, Parameters& pa = defaultParameters);
Real mu(const Point& P, Parameters& pa = defaultParameters);
Real beta_0(const Point& P, Parameters& pa = defaultParameters);

Real Fm (const Point&P, Parameters& pa = defaultParameters);
Real Gm (const Point&P, Parameters& pa = defaultParameters);
Real Fp (const Point&P, Parameters& pa = defaultParameters);
Real Gp (const Point&P, Parameters& pa = defaultParameters);
Real u_ex(const Point& P, Parameters& pa = defaultParameters);
Real interface_func(const Point&P, Parameters& pa = defaultParameters);

Real g_1D(const Point& P, Parameters& pa = defaultParameters);
Real h(const Real& t, Parameters& pa = defaultParameters);
Real f(const Point& P, Parameters& pa = defaultParameters);

Real splin_exp(const Real x);
Real splin_lin(const Real x);
Real raccord (const Real x, const Real p1,const  Real p2,const Real a, const Real b);
void simulation(path dossier_film, path dossier_trail);

Real f_2_temp(const Real x);        // par construction, cette fonction est 2T periodique
Real f_2(const Real x);             // passe d'une fonction definit sur [0,2T] a defini sur R+
Real f_1(const Real x);


class Domain_S
{
    public :
    Real c,sigma,rho,mu;
    vector<function<Real(Real,Real)>> bords;      // a valeurs dans {0,1}, si regularisation : [0,1] 
    Domain_S(Real r1, Real r2,String choix,Real _epsi =0)
    {
        if(choix =="vitesse"){  c= r1;sigma= r2; rho=sigma/c; mu=sigma*c;}
        else {rho=r1;mu=r2;  c=sqrt(r2/r1);sigma=sqrt(r1*r2);}
        epsi=_epsi;
    }

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

    void ajout_bords(function<Real(Real,Real)> bordure)
    {
        bords.push_back(bordure);
    }

    Real rho_S(Real x, Real t)
    {
        Real is_inside=1;
        for(int i=0; i<bords.size(); i++){ is_inside = is_inside*bords[i](x,t); }
        // std::cout<<"x= "<<x<<" "<<is_inside<<" "<<rho<<'\n';
        return rho*is_inside;
    }
    Real mu_S(Real x, Real t)
    {
        Real is_inside=1;
        for(int i=0; i<bords.size(); i++){ is_inside = is_inside*bords[i](x,t); }
        // std::cout<<"x= "<<x<<" "<<is_inside<<" "<<mu<<'\n';
        return mu*is_inside;
    }

};

// void ajout_interface(Domain_S _D1, Domain_S _D2,function<Real(Real,Real)> interface){
//     _D1.ajout_bords(interface); _D2.ajout_bords(interface);
// }

// class Interface_S
    // {
    //     public :
    //     Domain_S D1, D2;
    //     Real vit, Start;
    //     Real reg=0.;

    //     Interface_S(Domain_S _D1, Domain_S _D2, Real _vit, Real _Start, Real _Reg=0): D1(_D1), D2(_D2), vit(_vit),Start(_Start), reg(_Reg) {};

    //     // Real rho(Real x){
    //     //     cout<<"I was here \n";
    //     //     Real x = _P(1) -vit*t;
    //     //     cout<<x<<'\n';
    //     //     cout<<D1.rho<<'\n';
    //     //     if(x<0){return D1.rho;}
    //     //     if(x<reg){return raccord(x,D1.rho,D2.rho,-reg,reg);}
    //     //     else{return D2.rho;}
    //     //     return 1;
    //     // };
    //     // Real mu(Point _P){
    //     //     Real x = _P(1) -vit*t +Start;
    //     //     if(x<-reg){return D1.rho;}
    //     //     if(x<reg){return raccord(x,D1.rho,D2.rho,-reg,reg);}
    //     //     else{return D2.rho;}
    //     // };
    // };


extern vector<Domain_S> Domaines;
// extern vector<Interface_S> Interfaces;
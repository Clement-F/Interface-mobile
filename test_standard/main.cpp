#include "xlife++.h"
#include "simulation.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>

using namespace xlifepp;
using namespace std;
using namespace std::filesystem;
//============================= pre main ===============================
  // ========================== var global ==========================================

    Number counter =0;
    Number error_count =0;
    Number ord;

    Real tf, length;
    Real CFL;         // coefficient de cfl
    Real dspace;
        
    Real rho_0,mu_0 ;       // coeffs de l'espace Q- 
    Real rho_1,mu_1 ;       // coeffs de l'espace Q+
    Real cm ,cp ;            // vitesse en -/+
    Real sm ,sp ;            // permeabilité en -/+
    Real c_max;            // vitesse max (- et +)
    String mod;            // mode choisi (flat, jump, interface, periode, energie)
    String cond;           // condition au bord

                          // paramètres nécessaire aux fonctions  
    Real vit;              // vitesse de l'interface
    Real vit_abs;          // vitesse de l'interface en val absolue
    Real epsi;             // paramètre de régularisation de l'interface
    const Real wave_length=10;   // periode temporelle du mode periodique
    const Real space_length=10;  // periode spacial du mode periodique
    Point k = Point(1.);   // point fonctionnant comme vecteur de base
    Real T;                // periode temporelle du mode energie
    Real Start_int;        // endroit ou commence l'interface 
    Real Start_domaine;
    Real Start_sol;        // endroit ou commence la solution

    // definit la manière d'obtenir la valeur de la derivee en t=0
    // l'un des deux doit etre vrai, la valeur est privilegier a la derivee
    bool possede_derive=false , possede_valeur=true ;
    bool uni_dir ; Number direction;
    bool error_verb=1 , error_calc=1;

    Real t=0., dtime, trho=0.;   // paramètres temporelles

    Number savefile,  savepic,  vbl;
    
    FESubType fesub=_GaussLobattoPoints;
    String choix, fes = "Lobatto";
    Number test_nombre=0;
    
    Real R_ref ,T_trans;
    Real tau_1, tau_2, tau_3, tau_4;
    vector<Domain_S> Domaines;

  // ========================== tool ================================
    void is_stable(const Real dspace, const Real dt, const Real max_c){ // si le max est explicite / donné
        if(max_c * dt > dspace){theCout<<"methode instable \n"; exit(-1); }
    }
    Real give_rho(const Real sigma, const Real c){return sigma/c;}                // rend rho en fonction de sigma, c
    Real give_mu(const Real sigma,const Real c){return sigma*c;}                  // rend mu en fonction de sigma, c


int main(int argc, char** argv){
  init(argc, argv, _lang=en); // mandatory initialization of xlifepp  

  // ========================== parametre batterie de test =============================
    String directorypath = "test_stock";
    path dossier ="test_stock";
    
    if (!exists(dossier)) {create_directory(dossier);}
    ord =1; dspace = 0.002;
    tf = 1; length =2;
    savefile=400;  savepic =500;  verboseLevel(0);
    choix ="vitesse"; cond ="trans"; uni_dir =0; direction=1;
    epsi =0;
    Start_int=0;
    cm= 1.5; cp = 0.6; sm = 3; sp = 7;

    std::cout<<"attribution vecteurs \n";

    Number nb_point=9;
    int nb_test =100;

    vector<Real> A_sm(nb_test,sm), A_sp(nb_test,sp), A_cm(nb_test,cm), A_cp(nb_test,cp); 
    vector<Real> A_epsi(nb_test,epsi), A_dx(nb_test,dspace); vector<Number> A_ord(nb_test,1);
    vector<String> A_mod(nb_test,"interface"); 
    
    std::cout<<"attribution positions \n";

    vector<Real> s_standard = {-0.2  ,0.2  ,-0.2 ,0.1   ,-0.3 ,0.2  ,-0.2 ,0.3  ,-0.2 ,0.2  ,-0.4  ,0.2};
    vector<Real> d_standard = {-1.2  ,-1   ,-0.8 ,-0.8  ,-1.5 ,-0.8 ,-0.4 ,-0.4 ,-0.8 ,-0.4 ,-1.2 ,-0.6};
    vector<Real> V_standard = {0.4,0.4,-0.4,-0.4,-2,-2,2,2,1,-1,-1,1};
    vector<vector<Real>>  A_ss(nb_test,s_standard);
    vector<vector<Real>>  A_sd(nb_test,d_standard);
    vector<vector<Real>>  A_v(nb_test,V_standard);
    vector<Real> zeros(12,0.);

    std::cout<<"vecteurs changeant \n";

    //convergence de dx ->0
    for(int k =0;k<=2;k++){
    for(int i =1;i<=nb_point+1;i++){
      A_dx  [k*(nb_point+1)+i] = 0.05*pow(0.6,i); 
      A_ord [k*(nb_point+1)+i] = Number(k+1);
    }
    }

    // //convergence de epsi -> 0
    // for(int k=1;k<=5;k++){
    //   for(int i =30;i<=39;i++){  
    //   A_dx  [k*(nb_point+1)+i] = 0.05*pow(0.65,i);
    //   A_epsi[k*(nb_point+1)+i]=0.2*pow(0.5,k);}
    //   }
    // A_dx[0] = 0.001; A_epsi[0] = 0.01;

    theCout<<A_dx<<'\n'<<A_ord<<'\n'<<A_epsi<<'\n';
    // //test des derives de rho/mu
    // for(int i=16;i<=19;i++){A_mod[i]="Jump";} 
    // for(int i=20;i<=23;i++){A_mod[i]="interface";A_v[i]= zeros;} 
    // int i=16;
    // A_sm[i]=sqrt(12); A_sp[i] = sqrt(32); A_cm[i]=sqrt(4./3); A_cp[i]=sqrt(1./2); // mm=4=mp; rm =3,rp=8;
    // A_sm[i+1]=sqrt(12); A_sp[i+1] = sqrt(32); A_cm[i+1]=sqrt(1./2); A_cp[i+1]=sqrt(4./3); // mm=4=mp; rm =8,rp=3;
    // A_sm[i+2]=sqrt(12); A_sp[i+2] = sqrt(21); A_cm[i+2]=sqrt(4./3); A_cp[i+2]=sqrt(7./3); // mm=4, mp=7; rm =3=rp;
    // A_sm[i+3]=sqrt(21); A_sp[i+3] = sqrt(12); A_cm[i+3]=sqrt(7./3); A_cp[i+3]=sqrt(4./3); // mm=7, mp=4; rm =3=rp;
    // A_sm[i+4]=sqrt(12); A_sp[i+4] = sqrt(32); A_cm[i+4]=sqrt(4./3); A_cp[i+4]=sqrt(1./2); // mm=4=mp; rm =3,rp=8;
    // A_sm[i+5]=sqrt(12); A_sp[i+5] = sqrt(32); A_cm[i+5]=sqrt(1./2); A_cp[i+5]=sqrt(4./3); // mm=4=mp; rm =8,rp=3;
    // A_sm[i+6]=sqrt(12); A_sp[i+6] = sqrt(21); A_cm[i+6]=sqrt(4./3); A_cp[i+6]=sqrt(7./3); // mm=4, mp=7; rm =3=rp;
    // A_sm[i+7]=sqrt(21); A_sp[i+7] = sqrt(12); A_cm[i+7]=sqrt(7./3); A_cp[i+7]=sqrt(4./3); // mm=7, mp=4; rm =3=rp;

    // vector<Real> dV = {-0.05,-0.05,0.05,0.05,0.05,0.05,-0.05,-0.05,-0.05,0.05,0.05,-0.05};
    // for(int i=29;i<=36;i++){
    //   A_v[i]  = V_standard + (i-28)* dV;
    //   A_v[i+9] = V_standard - (i-28)* dV;
    // }
    
    std::cout<<"test begin \n";
    for(int k=0;k<nb_test;k++){

      theCout<<" \n ========================================================= \n";
      theCout<<" \n                      test "+tostring(k+1)+"                             \n";
      theCout<<" \n ========================================================= \n";

      mod=A_mod[k]; epsi= A_epsi[k]; dspace= A_dx[k];
      cm = A_cm[k]; cp= A_cp[k]; sm= A_sm[k]; sp= A_sp[k];
      ord = A_ord[k];

      mu_0 = give_mu(sm,cm); mu_1 = give_mu(sp,cp);
      rho_0 = give_rho(sm,cm); rho_1 = give_rho(sp,cp); 
      std::cout<<"order "<<ord<<'\n';
      std::cout<<"dx "<<dspace<<"\n";

      path dossier_test = dossier/path("test_"+tostring(k+1));
      if (!exists(dossier_test)) {create_directory(dossier_test);}
      test_standard(dossier_test, A_ss[k], A_sd[k],A_v[k]);
    }
   return 0;
  }


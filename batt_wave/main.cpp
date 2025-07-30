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
    bool possede_derive=false , possede_valeur=true , pre_Start;
    bool uni_dir ; Number direction;
    bool error_verb , error_calc;

    Real t=0., dtime, trho=0.;   // paramètres temporelles

    Number savefile,  savepic,  vbl;
    
    FESubType fesub;
    String choix, fes;
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
  //=========================== recupere les parametres externes==================================
    Options opts;
    // déclaration des options
      opts.addOption("order",1);      opts.addOption("dspace",0.01);
      opts.addOption("subtype","standard"); // standard or lobatto
      opts.addOption("tf",1.);        opts.addOption("length",1.);

      opts.addOption("mu_base",1.);   opts.addOption("rho_base",1.);
      opts.addOption("mu_change",1.); opts.addOption("rho_change",1.);
      opts.addOption("c_minus",1.);   opts.addOption("s_minus",1.);
      opts.addOption("c_plus",1.);    opts.addOption("s_plus",1.);
      
      opts.addOption("choix","vitesse");  opts.addOption("vit", 1.);
      opts.addOption("period_T", 1.); opts.addOption("regularisation", 0.1);
      opts.addOption("model","flat"); opts.addOption("condition","transparent");
      
      opts.addOption("savefile",1);   opts.addOption("savepic",1);    opts.addOption("verbose",0);
      opts.addOption("CFL",0.9);      opts.addOption("preStart",true);
      opts.addOption("Start_int",1.); opts.addOption("Start_sol",0.5); opts.addOption("Start_domaine",1.);

      opts.addOption("error_verb",true);  opts.addOption("error_calc",true);  opts.addOption("uni_dir",true);
      opts.addOption("direction","gauche");


      opts.parse("data.txt", argc, argv);

    //assigner les vars des options
      error_verb = opts("error_verb"); error_calc =opts("error_calc");
      ord=opts("order");  dspace=opts("dspace");
      tf=opts("tf");  length = opts("length");

      String fes_temp=opts("subtype");  fes = fes_temp;   fesub = _standard;
      if(fes!="standard"){ fesub =_GaussLobattoPoints;}

      savefile=opts("savefile");  savepic =opts("savepic");
      vbl=opts("verbose");  verboseLevel(vbl);
      CFL = opts("CFL");

      // paramètres de vitesse / milieux
        String choix_temp = opts("choix");  choix = choix_temp;
        if(choix=="vitesse"){
          cm = opts("c_minus"); cp = opts("c_plus");
          sm = opts("s_minus"); sp = opts("s_plus");
          mu_0 = give_mu(sm,cm); mu_1 = give_mu(sp,cp);
          rho_0 = give_rho(sm,cm); rho_1 = give_rho(sp,cp);}
        else {
          mu_0=opts("mu_base");   rho_0 =opts("rho_base");
          mu_1=opts("mu_change"); rho_1 =opts("rho_change");
          cm = sqrt(mu_0/rho_0); cp = sqrt(mu_1/rho_1);
          sm = sqrt(mu_0*rho_0); sp = sqrt(mu_1*rho_1);
        }

      String mod_temp = opts("model");  mod = mod_temp;
      pre_Start=opts("preStart");
      Start_int=opts("Start_int"); Start_domaine=opts("Start_domaine"); Start_sol=opts("Start_sol");
      cout<<"Start_sol : "<<Start_sol<<" Start_int : "<<Start_int<<" \n";

      vit = opts("vit");  vit_abs = abs(vit);

      uni_dir = opts("uni_dir");
      String direction_temp = opts("direction");  if(direction_temp != "droite"){direction = 1;} else direction =-1;

      epsi = opts("regularisation");  T = opts("period_T"); if(mod=="energie"){T = abs(length/vit);}
      String cond_temp=opts("condition"); cond = cond_temp;

  // ========================== parametre batterie de test =============================
    String directorypath = "batterie_test";
    path dossier ="batterie_test";
    
    if (!exists(dossier)) {create_directory(dossier);}
    ifstream myfile ("data_batterie.txt");
    string data; int line=0;
    bool read=false;
    String name_sim;
    while ( getline (myfile,data) ){
      line +=1;
      std::cout<<line<<"\n";
      Real para3,para4,para5,para6; 
      if(line>10){read=true;}
      if(line %10 ==1 and read)dspace = stod(data);
      if(line %10 ==2 and read)ord = stod(data);
      if(line %10 ==3 and read)para3 = stod(data);
      if(line %10 ==4 and read)para4 = stod(data);
      if(line %10 ==5 and read)para5 = stod(data);
      if(line %10 ==6 and read)para6 = stod(data);
      if(line %10 ==7 and read)vit = stod(data);
      if(line %10 ==8 and read)Start_domaine = stod(data);
      if(line %10 ==9 and read)Start_sol = stod(data);
      // std::cout<<read<<" "<<dspace<<'\n';
      // cout<<"Start_sol : "<<Start_sol<<" Start_int : "<<Start_int<<" \n";
      if(line %10 ==0 and read)
      { 
      if(choix=="vitesse")
      { cm= para3; cp=para4; sm=para5; sp=para6;
        mu_0 = give_mu(sm,cm); mu_1 = give_mu(sp,cp);
        rho_0 = give_rho(sm,cm); rho_1 = give_rho(sp,cp);}

      else
      { rho_0=para3; rho_1=para4;mu_0=para5;mu_1=para6;
        cm = sqrt(mu_0/rho_0); cp = sqrt(mu_1/rho_1);
        sm = sqrt(mu_0*rho_0); sp = sqrt(mu_1*rho_1);}

      name_sim = data;
      std::cout<<name_sim<<'\n';
      test_nombre +=1;

      path dossier_film = dossier/ path("film"+tostring(test_nombre));
      if (!exists(dossier_film)) {create_directory(dossier_film);}

      path dossier_trail = dossier/ path("trail"+tostring(test_nombre));
      if (!exists(dossier_trail)) {create_directory(dossier_trail);}

      vit_abs = abs(vit);
      if(vit_abs<min(cm,cp)){ 
        std::cout<<"vitesse subsonique \n";
        if(Start_sol<0){ tau_1 = (cp-vit)/(cm-vit); tau_2 = (cm +vit)/(cm-vit);  R_ref = (sm-sp)/(sm+sp); T_trans = 2*sm/(sm+sp);}
        if(Start_sol>0){ tau_1 = (cm+vit)/(cp+vit); tau_2 = (cp-vit)/(cp+vit);   R_ref = (sp-sm)/(sm+sp); T_trans = 2*sp/(sm+sp);}
        }
      if(vit_abs>max(cm,cp)){ 
        std::cout<<"vitesse supersonique \n";
        if(Start_sol<0)
        {
          R_ref = (sp-sm)/(2*sp), T_trans = (sp+sm)/(2*sp);
          tau_1 = (vit-cp)/(vit-cm), tau_2 = (cp-vit)/(vit+cm), tau_3 =(vit+cp)/(vit+cm), tau_4=(vit+cp)/(vit-cm);}
        
        if(Start_sol>0)
        {
          R_ref =(sm-sp)/(2*sm), T_trans =  (sp+sm)/(2*sm);
          tau_1 = (-vit-cm)/(-vit-cp), tau_2 = (cm+vit)/(-vit+cp), tau_3 =(-vit+cm)/(-vit+cp), tau_4=(-vit+cm)/(-vit-cp);}
        }

      // if((cp>vit_abs && vit_abs>cm) or (cm>vit_abs && vit_abs>cp)){ std::cout<<"vitesse transsonique \n";error_calc =0; error_verb=0;}
      // else error_calc =1; error_verb=1;

      simulation(dossier_film,dossier_trail);
      }
    }

    return 0;
  }


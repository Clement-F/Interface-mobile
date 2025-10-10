#include "simulation.hpp"

using namespace xlifepp;
using namespace std;
using namespace std::filesystem;


//========================================================================
//  code fait par Clement Fontan lors d'un stage de M1
//  Les variables ne sont pas les meilleurs choisi et les noms ne sont pas les meilleurs mais le programme fonctionne
//  On peut encore voir des vestiges de fonctions qui ont ete faite puis effacer (probablement le meme jour)
//  Malheureusement n'ayant plus acces a la librairie Xlife, je ne peux pas effacer de bout de code sans compromettre son bon fonctionnement
//  Je vous incombe donc la tache de modifier comme bon vous semble ce code
//  Voir le ReadMe pour plus de details sur le fonctionnement du code
// =======================================================================


// ========================== fonction files =============================
// move un fichier d'un point d'origine a un point de destination.
// cette commande peut s'utiliser dans le script cpp
void move(const string source, const string destination)
{
  char cmd[100];

  strcpy(cmd,"mv ");
  strcat(cmd,source.c_str());
  strcat(cmd," ");
  strcat(cmd,destination.c_str());
  system(cmd);
}
// ========================== tool ================================

// donne un pas de temps repondant a une condition CFL simple
Real give_dt(const Real dspace, const Real max_c){return (CFL/ord) * (dspace/max_c);}

// raccort C^infty entre 0,1 sur 0,1
Real splin_exp(const Real x)    
{ 
    return exp(-1/x)*(x>0) / (exp(-1/x)*(x>0)+ exp(-1/(1-x))*(x>0));
}

// raccort C^1 entre 0,1 sur 0,1
Real splin_lin(const Real x)    
{ 
    return x*(x>0)*(x<1) ;
}    

// raccord C infty entre p1 et p2 sur a,b
Real raccord (const Real x, const Real p1,const  Real p2,const Real a, const Real b)
{ 
    // theCout<<"Nealy there \n";
    Real Coeff = (x-a)/(b-a);
    return (1-splin_exp(Coeff))*p1 + splin_exp(Coeff)*p2;
} 

// ========================== condition initiales===================================================

// u(x,0) <- u_0(xi) , xi = x +ct, pa=(c)
Real u_0(const Point&P, Parameters& pa) 
{ 
  Real a =1000; // parametre de concentration de la gaussienne
  Real d = abs(P(1)-Start_sol); 
  Real R= 0.05; // source radius
  Real amp= 1; // source amplitude (constant power)
  return exp(-a*d*d);
}

// u'(x,0) ou u(x,dt) en fonction des possede_
Real u_d(const Point&P, Parameters& pa) 
{
  // rends une valeur exacte de u(x,dt) comme u_0(x-c*dt)
  Real a = 1000,  R= 0.05, amp= 1,d;
  if(uni_dir)d = abs(P(1)-Start_sol - cm*dtime); else d=  abs(P(1)-Start_sol);
  return exp(-a*d*d);

  // rends une approximation de u(x,dt) comme DL 
  // if(mod =="energie")
  // {
  //   possede_derive = true; 
  //   possede_valeur = false;
  //   return cm*2*x*a*exp(-a*(x-length/2)*(x-length/2));
  // }
}



// ========================== fonctions de bords pour le cas energie ====================================

// f_2_temp n est definie que sur [0,2T]
Real f_2_temp(const Real x)
{  
    if((0<=x) and (x<=T/2))     {return cm*T + vit*x;}
    if((T/2 <x) and (x<=3*T/2)) {return cm*T + (vit*T)/2 - vit*(x-T/2);}
    if((3*T/2 <x) and (x<=2*T) ){return cm*T - (vit*T)/2 + vit*(x-3*T/2);} 

    else {return 0;}
  if(error_verb) theCout<<"error !!!!!!!!! \n";
  return 0;
}

// f_2 est definie sur [0,+inf] par periodicite de f_2_temp
Real f_2(const Real x)
{ 
  // evite un overlap avec f_1 a l'origine
  if(0<=x && x<T/2) return cm*T +vit*T-vit*x;    

  Real t_mod = x - floor(x/(2*T)) * 2*T;

  if(t_mod<0 and error_verb) error("free_error"," temps negatif sur f_2 \n");

  return f_2_temp(t_mod);
  if(error_verb) theCout<<"error !!!!!!!!! \n";
  return 0;
}

// f_1 est definie sur [0,+inf] par periodicite de f_2_temp
Real f_1(const Real x)
{
  if(0<=x && x<=T)    return -vit_abs*T + vit_abs*x;
  if(T<x && x<=3*T/2) return vit*T - vit*x;
  else                return f_2(x)-cm*T;

  // une situation n a pas ete pris en compte ==> erreur
  if(error_verb) theCout<<"error !!!!!!!!! \n";
  return 0;
}


// ========================== choix paramètres ================================

// fonction rho(x,t)
// la valeur de t est code comme paramètre globale

Real rho(const Point& P, Parameters& pa)
{   
  // rho est constant sur le domaine
  if(mod=="flat"){return rho_0;}

  // rho change brutalement de valeur en un instant
  if(mod=="Jump")
  {
    if(trho>=tf/4){return  rho_1;} 
    else{return rho_0;}
  }

  // rho change de valeur de maniere periodique
  if(mod=="period"){return rho_0 + rho_1 * sin(wave_length*trho/(2*pi_))*sin(space_length*dot(k,P)/(2*pi_));}

  // probleme d'interface mobile avec fonction raccord si on "regularise" l'interface
  if(mod=="interface"){
    Real x =P(1) -vit*t-Start_int;
      if(x<-epsi){return rho_0;}
      if(x<epsi){return raccord(x,rho_0,rho_1,-epsi,epsi);}
      else{return rho_1;}
  }
  
  // probleme d'interface mobile, cas de l'energie croissante avec fonction raccord si on "regularise" l'interface
  if(mod=="energie"){
    Real x1 =f_1(t), x2= f_2(t);
    if(P(1)<x1-epsi) return rho_1; 
    if(x1-epsi<P(1) && P(1)<x1+epsi) return raccord(P(1),rho_1,rho_0,x1-epsi,x1+epsi);
    if(x1+epsi<P(1) && P(1)<x2-epsi) return rho_0;
    if(x2-epsi<P(1) && P(1)<x2+epsi) return raccord(P(1),rho_0,rho_1,x2-epsi,x2+epsi);
    if(x2+epsi<P(1)) return rho_1;

    else return 100;
  }

  // probleme du "slab"
  // 2 interfaces, l'une commenceant en 0 et l'autre en Start_int
  // possibilite de rajouter une regularisation sur la fonction 
  if(mod=="double"){
    Real x1 =P(1) -vit*t;
    Real x2 =P(1) -vit*t-Start_int;
    if(x1<0 && x2<0) return rho_0;
    if((x1<0 && x2>0) or (x1>0 && x2<0)) return rho_1;
    if(x1>0 && x2>0) return rho_0;
    else return 0.;
  } 
  

  // cas general de l'interface mobile
  // on definit un domaine i a l'aide de fonction f_i,j
  // on dit que (x,t) est dans le domaine i si f_i,j(x,t)>0 pour tout j
  // plus de details dans le cpp
  // si (x,t) n'appartient pas a un domaine, rho(x,t) = 0
  // d'ou l'importance de definir tout le demi plan (x,t)
  // si (x,t) appartient a plus d'un domaine (i,j), rho(x,t) = rho_i + rho_j
  if(mod=="domain")
  {
    Real RHO=0;
    for(int i=0;i<Domaines.size();i++){
      // cout<<i<<" "<<Domaines[i].rho_S(P(1),t)<<"\n";
      RHO += Domaines[i].rho_S(P(1),t);
    }
    // std::cout<<"rho="<<RHO<<"\n";
    return RHO;
  }
  

  if(error_verb)  theCout<<"error !!!!!!!!! rho \n";
  return 0;
}

// se construit comme la fonction rho
Real mu(const Point& P, Parameters& pa)
{
  if(mod=="flat"){return mu_0;}
  if(mod=="Jump")
  {
    if(t>= tf/4){return  mu_1;}
    else{return mu_0;}
  }


  if(mod=="interface"){
      Real x =P(1) -vit*t-Start_int;
      if(x<-epsi){return mu_0;}
      if(x<epsi){return raccord(x,mu_0,mu_1,-epsi,epsi);}
      else{return mu_1;}
      }

  if(mod=="period"){return mu_0 + mu_1 * sin(wave_length*t/(2*pi_))*sin(space_length*dot(k,P)/(2*pi_));}

  if(mod=="energie"){
    Real x1 =f_1(t), x2= f_2(t);
    if(P(1)<x1-epsi) return mu_1; 
    if(x1-epsi<P(1) && P(1)<x1+epsi) return raccord(P(1),mu_1,mu_0,x1-epsi,x1+epsi);
    if(x1+epsi<P(1) && P(1)<x2-epsi) return mu_0;
    if(x2-epsi<P(1) && P(1)<x2+epsi) return raccord(P(1),mu_0,mu_1,x2-epsi,x2+epsi);
    if(x2+epsi<P(1)) return mu_1;

    else return 100;
    }

  if(mod=="double"){
    Real x1 =P(1) -vit*t;
    Real x2 =P(1) -vit*t-Start_int;
    if(x1<0 && x2<0) return mu_0;
    if((x1<0 && x2>0) or (x1>0 && x2<0))  return mu_1;
    if(x1>0 && x2>0) return mu_0;
    else return 0.;
  }

  if(mod=="domain")
  {
    Real MU=0;
    for(int i=0;i<Domaines.size();i++){
      // cout<<i<<" "<<Domaines[i].mu_S(P(1),t)<<"\n";
      MU += Domaines[i].mu_S(P(1),t);
    }
    // std::cout<<"mu="<<MU<<"\n";
    return MU;
  }
  
    if(error_verb) theCout<<"error !!!!!!!!! mu \n";
  return 0;
}

// est defini sur le bord
// beta est deja pris comme fonction. d'ou beta_0 
Real beta_0(const Point& P, Parameters& pa)
{
  // condition transparente, la solution sort du domaine, pas de reflexion
  if(cond=="trans"){return sqrt(rho(P)*mu(P));} 

  // condition de Neumann, la solution rebondit sur le bord, reflexion totale.
  if(cond=="neumann"){return 0;}
  
  if(error_verb) theCout<<"error !!!!!!!!! beta \n";
  return 0;
}

// ========================== solution exacte =================================================== 

// rends les solutions exactes du probleme de l'interface mobile
// il est possible de les generalisee pour le probleme de "slab" ou de l'energie croissante

// due aux pertes de symmetries des regimes supersonique et transoniques, 
// il est necessaire de preciser si la solution commence a droite ou a gauche de l'interface.

// si il y a un probleme avec la solution exacte (call sur un probleme du regime transsonique ou mauvaise definition)
// la fonction rends 10 et envoie un message d'erreur pouvant indiquer ou se trouve l'erreur.

Real Fm (const Point&P, Parameters& pa)
{
  Real xi = P(1)-cm*t;
  Point P_xi = Point(xi);
  if(xi<=Start_int){return 0.5*u_0(P_xi);}
  else{
    if(mod=="flat"){return 0.5*u_0(P_xi);}

    if(mod=="interface"){
      if(vit_abs<min(cm,cp)){if(error_verb) theCout<<"ERROR_int_neg ";return 10;} 

      if(vit_abs>max(cm,cp)){
        if(Start_sol<Start_int){
          if(vit>=0){return 0;}
          if(vit<0){if(error_verb) theCout<<"ERROR_neg_fm_n";return 10;}}
        if(Start_sol>Start_int){
          if(vit<0){if(error_verb) theCout<<"ERROR_neg_fm_p";return 10;}
          if(vit>=0){return 0.5*(T_trans * u_0(P_xi/tau_3 + Start_int*(cp-cm)/(-vit-cm))+ R_ref*u_0(P_xi/tau_4 + Start_int*(cm+cp)/(cp+vit)));}
          }
        else if(error_verb) theCout<<"ERROR_nochoice";return 10;

        }
      if(cm<vit_abs<cp){
        if(vit<0){if(error_verb) theCout<<"ERROR_int_neg ";return 10;}
        if(vit>=0){1+1;}} //rajouter l'expression de m

      if(cp<vit_abs<cm){if(error_verb) theCout<<"ERROR_int_neg ";return 10;}

    }
    else {if(error_verb) theCout<<"ERROR_neg_fm_all";return 10;}
  }
  if(error_verb) theCout<<"ERROR_Fm";return 10;
}
Real Gm (const Point&P, Parameters& pa)
{
  Real xi = P(1)+cm*t;
  Point P_xi = Point(xi);

  if(xi<=Start_int){return 0.5*u_0(P_xi);}
  else{
    if(mod=="flat"){return 0.5*u_0(P_xi);}
    if(mod=="interface"){
      // if(abs((-P_xi+Start_int)(1)/tau_2 -Start_sol)<0.1) theCout<<P(1)<<" GM "<<'\n';
      if(vit_abs<min(cm,cp)){
        if(Start_sol<Start_int){return R_ref/2 * u_0(-P_xi/tau_2);}
        if(Start_sol>Start_int){return T_trans/2 * u_0(P_xi/tau_1);}}

      if(vit_abs>max(cm,cp)){
        if(Start_sol<Start_int){
          if(vit>=0){return 0;}
          if(vit<0){if(error_verb) theCout<<"ERROR_neg_gm_n";return 10;}}
        if(Start_sol>Start_int){
            if(vit<0){if(error_verb) theCout<<"ERROR_neg_gm_p";return 10;}
            if(vit>=0){return 0.5*(T_trans*u_0(P_xi/tau_1 + Start_int*(cm-cp)/(cm-vit)) +R_ref*u_0(-P_xi/tau_2 + Start_int*(cm+cp)/(cp-vit)));}}
        else if(error_verb) theCout<<"ERROR_nochoice";return 10;
        }

      if(cm<vit_abs<cp){1+1;}

      if(cp<vit_abs<cm){1+1;}
    }
    else {if(error_verb) theCout<<"ERROR_neg_gm_all";return 10;}
  }
  

  if(error_verb) theCout<<"ERROR_Gm";return 10;
}
Real Fp (const Point&P, Parameters& pa)
{
  Real xi = P(1)-cp*t;
  Point P_xi = Point(xi);
  // theCout<<tau_2<<" "<<tau_1<<" \n";
    if(xi>=Start_int){return 0.5*u_0(P_xi);}
    else{
      if(mod=="flat"){return 0.5*u_0(P_xi);}

      if(mod=="interface" ){
        // if(abs((P_xi+Start_int)(1)/tau_1 -Start_sol)<0.1) theCout<<P(1)<<" FP "<<'\n';
        if(vit_abs<min(cm,cp)){
          if(Start_sol<Start_int){return T_trans/2 * u_0((P_xi)/tau_1);}
          if(Start_sol>Start_int){return R_ref/2 *u_0(-P_xi/tau_2);}
          }

        if(vit_abs>max(cm,cp)){ 
          if(Start_sol<Start_int){
            if(vit>0){if(error_verb) theCout<<"ERROR_neg_fp_p";return 10;}
            if(vit<=0){return 0.5*(T_trans*u_0(P_xi/tau_1 + Start_int*(cm-cp)/(vit-cp)) +R_ref*u_0(-P_xi/tau_2 + Start_int*(cm+cp)/(cp-vit)));}}
          if(Start_sol>Start_int){
            if(vit<=0){return 0;}
            if(vit>0){if(error_verb) theCout<<"ERROR_neg_fp_p";return 10;}}
        }
        if(cm<vit_abs<cp){1+1;}

        if(cp<vit_abs<cm){1+1;}
      }
      else {if(error_verb) theCout<<"ERROR_pos";return 10;}
    }
  if(error_verb) theCout<<"ERROR_Fp";return 10;
}
Real Gp (const Point&P, Parameters& pa)
{
  Real xi = P(1)+cp*t;
  Point P_xi = Point(xi);
  if(xi>=Start_int){return 0.5*u_0(P_xi);}
  else{
    if(mod=="flat"){return 0.5*u_0(P_xi);}

    if(mod=="interface"){
      if(vit_abs<min(cm,cp)){if(error_verb) theCout<<"ERROR_int_pos";return 10;}

      if(vit_abs>max(cm,cp)){
        if(Start_sol<Start_int){
          if(vit>0){if(error_verb) theCout<<"ERROR_neg_gp_p";return 10;}
          if(vit<=0){return 0.5*(T_trans * u_0(P_xi/tau_3 +Start_int*(cp-cm)/(vit+cp))+ R_ref*u_0(P_xi/tau_4 + Start_int*(cm+cp)/(cp+vit)));}}
        if(Start_sol>Start_int){
          if(vit<=0){return 0;}
          if(vit>0){if(error_verb) theCout<<"ERROR_neg_gp_p";return 10;}}
      }
      else if(error_verb) theCout<<"ERROR_nochoice";return 10;
      
      if(cm<vit_abs<cp){1+1;}

      if(cp<vit_abs<cm){1+1;}
    }
    else {if(error_verb) theCout<<"ERROR_pos";return 10;}
  }
  if(error_verb) theCout<<"ERROR_Gp";return 10;
}

// combine les fonctions precedentes en une seule fonction donnant la solution exacte

Real u_ex(const Point& P, Parameters& pa)
{
  if(mod == "flat"){ return Fm(P)+Gm(P);}
  if(mod =="interface"){
  if(P(1)-vit*t<Start_int) return (Fm(P)+Gm(P));
  if(P(1)-vit*t>Start_int) return (Fp(P)+Gp(P));

  // dans le cas improbable que l'on demande la solution exactement a l'interface, rends le point milieu
  if(P(1)-vit*t==Start_int) return 0.5*(Fm(P)+Gm(P)) + 0.5 * (Fp(P)+Gp(P));}
  if(error_verb) theCout<<"ERROR_u_ex";return 10;
} 

// rends le domaine sur lequel on se trouve dans le cas de l'interface mobile

Real interface_func(const Point&P, Parameters& pa)
{
  if(mu(P)==mu_0 and rho(P)==rho_0) return 0;
  if(mu(P)==mu_1 and rho(P)==rho_1) return 1;
  else return 0.5;
}

// ========================== 2nd membre / terme source ===================================================
// tres peu tester, peut etre source de bug
Real g_1D(const Point& P, Parameters& pa){return 0;}
Real h(const Real& t, Parameters& pa){return 0;}
Real f(const Point& P, Parameters& pa){return 0;}

// fait tout la simulation en une seule fonction
// necessaire de definir tout les parametre globaux avant de call la fonction
// si un parametre est mal definit, la fonction va gueuler. 

void simulation(path dossier_trail)
{ 

  //=========================== check data ==============================================
  // Liste non exhaustive de probleme de parametre
    if(error_verb or error_calc){

      // le temps final n'est pas avant le debut de la simulation
      if(tf<0) error("temps de simulation négatif \n");

      // le domaine possede une longueur positive
      if(length<0) error("domaine de simulation negatif \n");

      // le pas d'espace est positive et inferieur a la longueur du domaine
      if(dspace<0 or dspace>length) error("pas d'espace mal defini \n");

      // l'ordre de la methode est positive
      if(ord<0) error("ordre mal defini \n");

      // la methode n est pas defini si l'on utilise une methode autre que Gauss Legendre (standard) ou Gauss Lobatto
      // !!!! il faut refaire le code pour faire avec la methode Gauss Legendre, pour le moment seul Gauss Lobatto est codee -> ligne 608
      if(fes!="standard" and fes!="Lobatto") error("subtype mal defini \n");

      // on sauvegarde tout les 'images positive'
      if(savefile<0) error("savefile mal defini \n");

      // on sauvegarde au plus 'images' positive
      if(savepic<0) error("savepic mal defini \n");

      // les parametre physiques sont positif
      if(mu_0<0 or mu_1<0)  error("mu negatif \n");
      if(rho_0<0 or rho_1<0) error("rho negatif \n");

      if(cm<0 or cp<0) error("celerite negatif \n");
      if(sm<0 or sp<0) error("permeabilite negatif \n");

      // une periode temporelle est positive
      if(T<0) error("period negative \n");

      // une regularisation est positive
      if(epsi<0) error("regularisation negative \n");

      // on definit les parametre physique soit par rapport aux vitesses :
      // mu,rho en fonction de sigma et c
      // ou en fonction des "parametres"
      // sigma et c en fonction de mu,rho 
      if(choix !="vitesse" and choix!="para") error("choix non defini \n");

      // determine le choix des fonctions rho et mu
      // si l'on rajoute une methode, ajouter le check
      if(mod !="flat" and mod !="Jump" and mod!="interface" and mod!="period" and mod!="energie" and mod!="double" and mod !="domain") error("model non supportee");

      // determine les conditions au bord
      if(cond !="neumann" and cond!="trans") error("condition non supportee \n");

      // le coefficient CFL est dans [0,1]
      if(0>CFL or CFL>1) error("CFL mal defini");  }
    
  //=========================== declaration ===================================================
    
  // a defaut de mieux, voici les coefficients CFL empirique en fonction de l'ordre.
    if(false){
    if(ord ==1){CFL = 0.95;}
    if(ord ==2){CFL = 0.8;}
    if(ord ==3){CFL = 0.6;}
    if(ord >=4){CFL = 0.5;}
  }
    theCout<<"la CFL est à "<<CFL<<'\n';

    // xlife things
      // mesh
      // nombre de nodes/noeux
      Number nodes = (int)length/dspace;

      // definit le domaine sur lequel se pose la solution
      // on peut le remplacer par un domaine 2D simplement
      Segment seg(_xmin=Start_domaine, _xmax=Start_domaine+length, _nnodes=nodes, _domain_name="Omega", _side_names="Gamma");

      // definit le maillage
      Mesh mesh1d(seg, _generator=structured, _name="P1-mesh");

      // definit le domaine de solution a partir du mesh
      // Omega le domaine interieur, Gamma le domaine du bord
      Domain omega=mesh1d.domain("Omega"), gamma=mesh1d.domain("Gamma");

      // interpolation
      // espace de solution / Galerkin
      Space V(_domain=omega, _FE_type=Lagrange, _FE_subtype=fesub, _order=ord ,_name="V");

      // inconnue et fonction test
      // important dans Xlife
      Unknown u(V, _name="u");  TestFunction v(u, _name="v");

      theCout<<"========================= Mesh & space created ========================"<<'\n';

    // time parameters
    
      // Grand ennoncee des parametres du probleme
      // Si il y a un probleme avec le code, premier endroit ou regarder
      // Aucun parametre ne devrait etre 0
      c_max = max(cm,cp);
      theCout<<" model : "<<mod<<" \n";
      if(mod=="interface"){
        if(vit_abs<min(cm,cp)) theCout<<" régime subsonique \n";
        else if(vit_abs>max(cm,cp)) theCout<<" régime supersonique \n";
        else if(vit_abs>cp) theCout<<" régime transonique +\n";
        else if(vit_abs>cm) theCout<<" régime transonique -\n";
      }
      theCout<<" célérité maximal : "<<c_max;
      dtime = give_dt(dspace,c_max);
      Number nbt = (int)(tf/dtime);
      theCout<<" \n nombre de pas de temps :"<<nbt<<" pas de temps :"<<dtime<<"\n";
      theCout<<" les bords du domaines sont : "<<Start_domaine<<" et "<<Start_domaine+length<<'\n';
      theCout<<" la solution commence en :"<<Start_sol<<'\n';
      theCout<<" nombre de noeuds : "<<length/dspace *ord<<"\n";
      theCout<<" ordre de la méthode : "<<ord<<" \n";
      theCout<<" paramètre de régularisation :"<<epsi<<" \n";

      // Nous evite de recalculer ces constantes, fait gagner du temps de calcul
      Real dt2=dtime*dtime, dt_2= dtime/2, dt_4d3 =1/(4*dtime*dtime*dtime);

      // le nombre de frame ou l'on souhaite avoir des infos
      // savepic est une frequence de frame ou l'on souhaite avoir des infos
      savefile = min(savefile,Number(ceil(Real(nbt)/savepic)));
      theCout<<" savefile ="<<savefile<<" \n";

    // storage creation
      // Stock la solution approchee u en t
      TermVectors U(nbt, _name="U_sol");  

      // Stock la solution u exacte en t
      TermVectors U_exacts((int)(nbt/savefile)+1, _name="U_ex");

      // Stock l'etat des domaines en t
      TermVectors Interface((int)(nbt/savefile)+1, _name="Inter");

      // Stock la valeur de Rho et Mu en t sur le domaine
      TermVectors RHO((int)(nbt/savefile)+1, _name="RHO");
      TermVectors MU((int)(nbt/savefile)+1, _name="MU");

      // Stock l'energie au temps t
      std::vector<Real> E((int)(nbt/savefile)+1); 

      // Stock les erreur L2 et H1 au temps t
      std::vector<Real> ErrorL2((int)(nbt/savefile)+1), ErrorL2_relat((int)(nbt/savefile)+1);
      std::vector<Real> ErrorH1((int)(nbt/savefile)+1), ErrorH1_relat((int)(nbt/savefile)+1);
      theCout<<"========================= Storage created ========================"<<'\n';

    // cond initial

      // vecteur de bases definit sur le domaine omega
      TermVector zeros(u, omega, 0.), ones(u,omega,1.);

      // Conditions initiales
      TermVector U0(u,omega,u_0,_name="U_sol"),  U1(u,omega,u_d,_name="U_sol"); 
    
      // max de l'erreur sur la simulation
      Real max_Error_L2=0., max_Error_relat_L2=0.;
      Real max_Error_H1=0., max_Error_relat_H1=0.;
    
      // initialise U
      U(1)= U0;

      // theCout<<"normes des premiers termes "<<norminfty(U0)<<" "<<norminfty(U1)<<endl;
      theCout<<"========================= Cond init set ========================"<<'\n';

      // definition des formes
      t=0;
      trho = t+dt_2; 

      // !!!! les formes ne sont definit que pour GaussLobatto
      // pour avoir du Gauss Lagrange (standard), il faut redefinir cette partie
      BilinearForm a = intg(omega, mu*grad(u)|grad(v), _quad = GaussLobatto, _order =2*ord-1);
      BilinearForm m = intg(omega, rho*u*v, _quad = GaussLobatto, _order =(2*ord-1));
      BilinearForm b = intg(gamma, beta_0*u*v, _quad = GaussLobatto, _order =(2*ord-1));
      LinearForm   l = intg(omega, f*v, _quad = GaussLobatto, _order =2*ord-1);

      // declaration des differentes matrices EF
      TermMatrix Masse = TermMatrix(intg(omega, u*v, _quad = GaussLobatto, _order =(2*ord-1)));   
      TermMatrix Rigid = TermMatrix(intg(omega, grad(u)|grad(v) , _quad = GaussLobatto, _order =(2*ord-1)));
      TermMatrix An =   TermMatrix(a,_name="An_"+tostring(t)); 
      TermMatrix Mnp1 = TermMatrix(m,_name="M_"+tostring(t)); 

      TermMatrix L;     ldltFactorize(Mnp1, L);

      // approximation a l'ordre 2 de la solution au temps dt
      TermVector U2 =   factSolve(L, dt2*U0);

      // si la solution se dirige dans une seul direction ( la derivee est non nulle ), l'approximation est fausse -> solution exacte par u_d
      // si la solution ne se dirige pas dans une seul direction ( la derivee est nulle ), l'approximation est d'ordre 2 -> on l'utilise
      // oui le noms des variables sont mal choisi
      if(uni_dir){U(2)=U1; theCout<<"hiiiiiiiiii"<<'\n';}
      else{U(2) = U0 - 0.5* U2 ;theCout<<"heyyyyyyyyyy"<<'\n';}
      
      // j'etais curieux de la taille de l'approximation
      saveToFile("U1",U(1)-U(2),_format=raw);

      theCout<<"========================= Form defined ========================"<<'\n';
    
    // derniere preparations
    
      TermMatrix Mn(m,_name="M_t"+tostring(t)); //Mn-1/2
      TermMatrix Anm1(a, _name="A"+tostring(t));

      Real tmat=0., tsmb=0., tinv=0., tsol=0., tproc=0.; //temps de calculs

      // derniere declaration du probleme
      // si il y a un probleme, c'est le dernier moment pour le voir
      // aucun parametre ne devrait etre 0
      theCout<<"paramètres utilisée (rho, mu, c, sigma) : \n";

      theCout<<"milieu 0: "<<rho_0<<" ;"<<mu_0<<" ;"<<cm<<" ;"<<sm<<" \n";
      theCout<<"milieu 1: "<<rho_1<<" ;"<<mu_1<<" ;"<<cp<<" ;"<<sp<<" \n";

      if(mod=="interface") theCout<<"vitesse = "<<vit<<" \n";
      if(mod=="period") theCout<<"period spacial = "<<space_length<<" period temporel = "<<wave_length<<" \n";
      if(mod=="energie") theCout<<"vitesse = "<<vit<<", period temporel = "<<T<<" \n";
        theCout<<"les coeffs de reflexion, transmissions sont : R ="<<R_ref<<" T="<<T_trans<<" \n";
        theCout<<"les coeffs de contraction dilatation sont : tau_1="<<tau_1<<" tau_2="<<tau_2<<" tau_3="<<tau_3<<" tau_4="<<tau_4<<" \n";

      // declaration des dernier matrices que l'on va utiliser.
      TermMatrix M1_temp, M2_temp, M_moy;
      TermMatrix Bn, MB;
      TermVector Fn, G, Err, U_1, U_2, MBv;    
  // ========================== loop ============================================================== 
    theCout<<"========================= loop ========================"<<'\n';
    for (Number n=2; n<nbt; n++, t+=dtime)
    {    
      // déclaration des matrices en temps n

      trho = t+dt_2; 
      elapsedTime();
      An=   TermMatrix(a,_name="An_"+tostring(t)); 
      Mnp1= TermMatrix(m,_name="M_"+tostring(t)); 

      trho=t;
      Bn=   TermMatrix(b,_name="Bn_"+tostring(t));
      tmat+=elapsedTime();

      // inversion de la matrice MB
      MB = Mnp1 + dt_2*Bn ;
      if(fesub==_standard) ldltFactorize(MB, L);  else L=inverse(MB);

      // tentative "d'accelerer" le calcul
      // il est possible que ce soit exactement le meme calcul qui est fait en arriere plan lorsque l'on demande inverse(MB)

      // if(MB.isDiagonal())
      // {
      //   MBv = MB*ones;
      //   MBv = TermVector(MBv, 1./x_1);
      //   L = TermMatrix(MBv,_name="MB"+tostring(t));
      // }
      // else {ldltFactorize(MB, L);}

      L.name("L");
      tinv+=elapsedTime();

      // schema numerique
      Fn = TermVector(l,_name="Fn_"+tostring(t));
      G = - Mn*U(n-1) + dt_2*Bn*U(n-1) + (Mnp1+Mn)*U(n) + dt2*(Fn-An*U(n));
      tsmb+=elapsedTime();

      // resolution du systeme implicite
      U(n+1)=L*G;
      tsol+=elapsedTime();


      if((n)%savefile ==0){      //traitement de donnees
        // calcul de l'energie au temps t
        E[(int)n/savefile] = abs(0.5/dt2 * (Mnp1 * U(n+1)-Mnp1*U(n))|(U(n+1)-U(n))) +0.5 * abs(An*U(n+1)|U(n)); 

        theCout<<"========= "<<test_nombre<< " =============== "<<n<<" =============== " <<(Real(n)/nbt)*tf<<" ========================= \n";

        // regarde quelques parametres interessant nous permettant de juger de la stabilite de la solution rapidement
        // cela ne garantie aucunement de la regularite de la solution et encore moins de sa validite
        // juste qu'elle n'explose pas
        Real normU = norminfty(U(n));
        theCout<<"amplitude numerique : "<<normU<<" \n";
        theCout<<"energie             : "<<E[(int)n/savefile]<<'\n';
        theCout<<floor(Real(n)/nbt*10000)/100<<"%"<<" du calcul fini: "<<n<<"/"<<nbt<<" \n";
        }
        tproc+=elapsedTime();

        // les matrices n+1 au temps n deviennent les matrices n au temps n+1
        // evite de les recalculee 
        Mn=Mnp1;
        Anm1 = An;
    }

  // ========================== Out of the loop============================================================

    // quelques parametres utiles pour voir si il y a eu des problemes
    theCout<<"le counter est à : "<<counter<<"\n";
    theCout<<"il y a eu : "<<error_count<<" erreurs \n";

    // max_Error_L2 n'est actuellement pas update dans la boucle, il faudrait l'update dans la partie traitement de donnees
    theCout<<"l'erreur maximal est : "<<max_Error_L2<<" \n"; 
    
    if(savefile>0)
    {
      elapsedTime();
      Real tsave=0.;
      theCout<<"========================= save to file ========================"<<'\n';      
      for (Number n=savefile, indice =1; n<=nbt; n+=savefile, indice+=1)  // la boucle fait un lien implicite entre savefile t et n
      { 
        t = n*dtime -dtime-dt_2  ; trho = t+dt_2;
        theCout<<"t ="<<t<<" trho ="<<trho<<" n="<<n<<" indice ="<<indice<<'\n';

        // calcule la solution exacte au temps t 
        U_exacts(indice) = TermVector(u,omega,u_ex);

      // calcule l'erreur au temps t
        Err =U(n)-U_exacts(indice);                                                               // U-U_ex

        // calcules des erreurs L2
        ErrorL2[indice] = abs((Masse*Err|Err));                                                   // int(M*(U-U_ex)*(U-U_ex))
        ErrorL2_relat[indice] = ErrorL2[indice]/abs((Masse*U_exacts(indice)|U_exacts(indice)));   // int(M*(U-U_ex)*(U-U_ex))/ int(M * U_ex*U_ex)
        if(ErrorL2[indice]>max_Error_L2) max_Error_L2 = ErrorL2[indice];                          // max_t int(M*(U-U_ex)*(U-U_ex))
        if(ErrorL2_relat[indice]>max_Error_relat_L2) max_Error_relat_L2 = ErrorL2_relat[indice];  // max_t int(M*(U-U_ex)*(U-U_ex))/ int(M * U_ex*U_ex)
        
        // calcules des erreurs H1
        ErrorH1[indice] = abs((Rigid*Err|Err));                                                   // int(K*(U-U_ex)*(U-U_ex))
        ErrorH1_relat[indice] = ErrorH1[indice]/abs((Rigid*U_exacts(indice)|U_exacts(indice)));   // int(K*(U-U_ex)*(U-U_ex))/ int(K * U_ex*U_ex)
        if(ErrorH1[indice]>max_Error_H1) max_Error_H1 = ErrorH1[indice];                          // max_t int(K*(U-U_ex)*(U-U_ex))
        if(ErrorH1_relat[indice]>max_Error_relat_H1) max_Error_relat_H1 = ErrorH1_relat[indice];  // max_t int(K*(U-U_ex)*(U-U_ex))/ int(K * U_ex*U_ex)


        // calcules des domaines / valeur des parametres au temps t
        Interface(indice) = TermVector(u,omega,interface_func);
        RHO(indice) = TermVector(u,omega,rho);
        MU(indice)  = TermVector(u,omega,mu);

        TermVector U_tosave = U(n); 


        // sauvegarde les solutions au temps t
        // raw ne prends que les valeurs et .vtu est un format particulier permettant une exportation vers paraview
        // U_exacts(indice) = TermVector(u,omega,u_ex);
        U_tosave.name("Usol"); Interface(indice).name("Int"); U_exacts(indice).name("Uex");
        RHO(indice).name("rho"); MU(indice).name("mu");
        saveToFile("U"+tostring(indice)+"_trail",U(n),_format=raw);
        saveToFile("Uex"+tostring(indice)+"_trail",U_exacts(indice),_format=raw);
        saveToFile("Int"+tostring(indice)+"_trail",Interface(indice),_format=raw);
        saveToFile("Rho"+tostring(indice)+"_trail",RHO(indice),_format=raw);
        saveToFile("Mu"+tostring(indice)+"_trail",MU(indice),_format=raw);


        // place les fichiers dans un dossier pre defini

        // move("U"+tostring(indice)+"_film.vtu", dossier_film);
        // move("Interface"+tostring(indice)+"_film.vtu", dossier_film);
        move("U"+tostring(indice)+"_trail.txt", dossier_trail);
        move("Uex"+tostring(indice)+"_trail.txt", dossier_trail);
        move("Int"+tostring(indice)+"_trail.txt", dossier_trail);
        move("Rho"+tostring(indice)+"_trail.txt", dossier_trail);
        move("Mu"+tostring(indice)+"_trail.txt", dossier_trail);
        theCout<<"------ post traitement-------"<<n<<"/"<<nbt<<" \n";
      } 

      // Sauvegarde et place divers fichier avec parametres
      std::ofstream myfile;
      myfile.open("E_trail.txt");
      for (Number i=0; i<=((int)nbt/savefile - 1) ;i++ ){myfile<<E[i]<<"\n";}
      myfile.close();         
      move("E_trail.txt", dossier_trail);

      myfile.open("ErrorL2_trail.txt");
      for (Number i=0; i<=((int)nbt/savefile - 1) ;i++ ){myfile<<ErrorL2[i]<<" "<<ErrorL2_relat[i]<<"\n";}
      myfile.close();         
      move("ErrorL2_trail.txt", dossier_trail);

      myfile.open("ErrorH1_trail.txt");
      for (Number i=0; i<=((int)nbt/savefile - 1) ;i++ ){myfile<<ErrorH1[i]<<" "<<ErrorH1_relat[i]<<"\n";}
      myfile.close();         
      move("ErrorH1_trail.txt", dossier_trail);


      // compiles divers parametres du probleme
      // donnees physiques, numeriques et temps de calcul pour une analyse simple 
      myfile.open("data_trail.txt");
      myfile<<"ordre: "<<ord<<'\n'<<"dx: "<<dspace<<'\n';
      myfile<<"c_m: "<<cm<<" \n"<<"c_p: "<<cp<<" \n"<<"v: "<<vit<<" \n"<<"s_m: "<<sm<<" \n"<<"s_p: "<<sp<<" \n";
      myfile<<"length: "<<length<<" \n"<<"tf: "<<tf<<'\n';
      myfile<<"start_interface: "<<Start_int<<" \n"<<"start_solution: "<<Start_sol<<'\n';

      tsave+=elapsedTime();
      myfile<<"mat: "<<tmat<<"\ninv: "<<tinv<<"\nsmb: "<<tsmb<<"\nsol: "<<tsol<<"\ntproc: "<<tproc<<"\ntsave: "<<tsave<<"\ntotal: "<<tmat+tsmb+tinv+tsol+tproc+tsave<<'\n';
      myfile<<"erreur_maximal_L2: "<<max_Error_L2<<"\nerreur_maximal_relative_L2: "<<max_Error_relat_L2<<" \n";
      myfile<<"erreur_maximal_H1: "<<max_Error_H1<<"\nerreur_maximal_relative_H1: "<<max_Error_relat_H1<<" \n";
      myfile<<"start_domaine: "<<Start_domaine<<'\n'<<"regularisation: "<<epsi<<'\n';;
      myfile.close();
      move("data_trail.txt", dossier_trail);

      theCout<<"===================================== \n";
      theCout<<"time leapfrog : mat = "<<tmat<<" inv = "<<tinv<<" smb = "<<tsmb<<" sol = "<<tsol<<" tproc="<<tproc<<" tsave ="<<tsave<<" total = "<<tmat+tsmb+tinv+tsol+tproc+tsave<<'\n';
      theCout<<"erreur maximal L2 : "<<max_Error_L2<<" erreur maximal relative L2 :"<<max_Error_relat_L2<<" \n";
      theCout<<"erreur maximal H1: "<<max_Error_H1<<"  erreur_maximal relative_H1: "<<max_Error_relat_H1<<" \n";
    }
}



// Pour faire un test sur plusieurs parametres, ou un test sur un parametres avec plusieurs ordre de grandeurs (convergence)
// le test standard prends divers positions initiales de l'interface S_d et de la solution S_s et divers vitesses de l'interface V
// et rends un test sur les 12 cas test de l'interface mobile a vitesse constante. 

// si le mod n'est pas interface, il prends les meme positions et vitesses et test un seul probleme.

void test_standard(path dossier,const vector<Real> S_s, const vector<Real> S_d,const vector<Real> V){
  

  if(mod=="interface")
  {
  string regime;
  for(int k=0; k<12;k++){
    test_nombre +=1;
    vit=V[k]; vit_abs = abs(vit);
    Start_domaine=S_d[k]; Start_sol= S_s[k];
    
    if(vit_abs<min(cm,cp)){ 
        if(Start_sol<0) {tau_1 = (cp-vit)/(cm-vit), tau_2 = (cm +vit)/(cm-vit);  R_ref = (sm-sp)/(sm+sp), T_trans = 2*sm/(sm+sp);}
        if(Start_sol>0) {tau_1 = (cm+vit)/(cp+vit), tau_2 = (cp-vit)/(cp+vit);   R_ref = (sp-sm)/(sm+sp), T_trans = 2*sp/(sm+sp);}
        }
    if(vit_abs>max(cm,cp)){ 
        if(Start_sol<0)
        {
          R_ref = (sp-sm)/(2*sp), T_trans = (sp+sm)/(2*sp);
          tau_1 = (vit-cp)/(vit-cm), tau_2 = (cp-vit)/(vit+cm), tau_3 =(vit+cp)/(vit+cm), tau_4=(vit+cp)/(vit-cm);}
        
        if(Start_sol>0)
        {
          R_ref =(sm-sp)/(2*sm), T_trans =  (sp+sm)/(2*sm);
          tau_1 = (-vit-cm)/(-vit-cp), tau_2 = (cm+vit)/(-vit+cp), tau_3 =(-vit+cm)/(-vit+cp), tau_4=(-vit+cm)/(-vit-cp);}
      }
    
    if(k<4){regime ="subsonique";}
    else if(k<8){regime = "supersonique";}
    else {regime ="sonique";}
    if(k<8){error_calc=1; error_verb=1;}
    if(k>=8){error_calc=0; error_verb=0;}
    path dossier_test = dossier/path("trail_" + tostring(k+1));
    if (!exists(dossier_test)) {create_directory(dossier_test);}
    simulation(dossier_test);
  }
  }

  else
  {
  path dossier_test = dossier/path("trail_"+mod);
  if (!exists(dossier_test)) {create_directory(dossier_test);}
  simulation(dossier_test);
  }
  
}
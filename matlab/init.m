function [data_phys, data_num, data_time]=init()
file_data=["data_trail.txt"];
fileID=fopen(file_data,'r');
data = textscan(fileID, '%s %f');
fclose(fileID);

data_string= data{1,1}; data_numb=data{1,2};

#================== récupère les données de data ==============

ord = data_numb(1); dx= data_numb(2);c_m= data_numb(3);c_p= data_numb(4);
v= data_numb(5);s_m= data_numb(6);s_p= data_numb(7);length= data_numb(8);
tf= data_numb(9);start_interface= data_numb(10);start_solution= data_numb(11);
mat= data_numb(12);inv= data_numb(13);smb= data_numb(14);sol= data_numb(15);
tproc= data_numb(16);tsave= data_numb(17);total= data_numb(18);
erreur_maximal= data_numb(19);erreur_maximal_relative= data_numb(20);
start_domaine=data_numb(21); epsi=data_numb(22);

#==============================================================
CFL=0.9;
if(ord==1) CFL=0.95; endif;
if(ord==2)CFL = 0.8; endif;
if(ord==3)CFL = 0.6; endif;
if(ord>=4)CFL = 0.5; endif;
c_max= max(c_m,c_p);
dt = (CFL/ord) * (dx/max(c_m,c_p));
nbt = round(tf/dt);

#-------------------------

file_data=["E_trail.txt"];
fileID=fopen(file_data,'r');
Eneg= textscan(fileID, '%f');
nbt_file = size(Eneg{:, 1},1);
fclose(fileID);

#-------------------------

filename = ["int1_trail.txt"];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, "%f");
fclose(fileID);
nbx=size(dataArray{:, 1},1);

#-------------------------
data_phys=[[v,c_max];[c_m,c_p];[s_m,s_p];[length,tf];[epsi,0]];
data_num=[[CFL,ord,0];[dx,dt,0];[nbx,nbt,nbt_file];[start_interface,start_solution,start_domaine];[erreur_maximal,erreur_maximal_relative,0]];
data_time=[ mat,inv,smb,sol,tproc,tsave,total];


function plot_energ_err(data_phys,data_num)

tf = data_phys(4,2);
nbt= data_num(3,2);nbt_file= data_num(3,3);
dt=data_num(2,2);
delimiter = {''};
formatSpec = '%f';

filename = ['E_trail.txt'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec);
fclose(fileID);

figure;
formatSpec = '%f %f';
savefile=round(tf/(dt*(nbt_file-1)));
ts= linspace(0,tf,nbt_file);
y_en=dataArray{:, 1};
plot(ts,y_en);
title(["graphe de l'Energie en fonction du temps"]);

hold off

figure;
filename = ['Error_trail.txt'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec);
fclose(fileID);

y_err=dataArray{:, 1};
y_err_relat=dataArray{:, 2};
hold on
plotyy(ts,y_err_relat, ts, y_err);
title(["graphe de l'erreur en fonction du temps"]);
hold off

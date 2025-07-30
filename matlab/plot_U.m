function plot_U(data_phys,data_num)

length = data_phys(4,1); tf = data_phys(4,2);
dx = data_num(2,1);     dt=data_num(2,2);
nbx= data_num(3,1); nbt= data_num(3,2);nbt_file= data_num(3,3);
start_int= data_num(4,1); start_sol = data_num(4,2); start_dom = data_num(4,3);
vit = data_phys(1,1);
val=zeros(1,nbx);

end_int = start_int + vit*tf;


delimiter = {''};
formatSpec = '%f';
filename = ['U1_trail.txt'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec);
fclose(fileID);
val(1,:) = dataArray{1, :};

for k=2:nbt_file
 filename = ['U' num2str(k) '_trail'  '.txt'];
 fileID = fopen(filename,'r');
 dataArray = textscan(fileID, formatSpec);
 fclose(fileID);
 val=[val;dataArray{1, :}'];
end

#T_end = 0.5*(1 /(5*1.)) ;

max_val = max(max(val));
min_val = min(min(val));

end_dom= start_dom+length-dx;
ts= linspace(0,tf,nbt_file);
xs=start_dom:length*dx:end_dom;

[X,T]=meshgrid(xs,ts);

figure('color','white');

surf(X,T,val);
line([start_int end_int],[0 tf],[max_val max_val], 'Color', 'm', 'Linewidth',.5);

shading interp
view(2)
axis image
colormap jet;
caxis([-max_val ,max_val]);
title(['graphe des caractéristiques']);
ylim([0,tf]); xlim([start_dom,end_dom]);

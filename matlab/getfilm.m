clear all


[data_phys, data_num, data_time]=init();

length = data_phys(4,1); tf = data_phys(4,2);
dx = data_num(2,1);     dt=data_num(2,2);
nbx= data_num(3,1); nbt= data_num(3,2);nbt_file= data_num(3,3);
start_int= data_num(4,1); start_sol = data_num(4,2); start_dom = data_num(4,3);
vit = data_phys(1,1); epsi=data_phys(5,1);
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
endfor

#T_end = 0.5*(1 /(5*1.)) ;

max_val = max(max(val));
min_val = min(min(val));

end_dom= start_dom+length-dx;
ts= linspace(0,tf,nbt_file);
xs=start_dom:length*dx:end_dom;

out_dir = "temp_img";
mkdir (out_dir);

for k=1:nbt_file
  t= k*(nbt/nbt_file)*dt;
  plot(xs,val(k,:), 'k-');
  line([vit*t+epsi,vit*t+epsi],[min_val,max_val]);
  line([vit*t-epsi,vit*t-epsi],[min_val,max_val]);
  title(["temps =" ,num2str(t)]);
  caxis([-max_val ,max_val]);
  ylim([0,tf]); xlim([start_dom,end_dom]);

  fname = fullfile (out_dir, sprintf ("img%03d.png", k));
  imwrite (getframe (gcf).cdata, fname);
endfor
cmd = sprintf ("ffmpeg -framerate 20 -i ./%s/img%%03d.png -vf scale=1080:-1 Example1.avi", out_dir)
system (cmd)


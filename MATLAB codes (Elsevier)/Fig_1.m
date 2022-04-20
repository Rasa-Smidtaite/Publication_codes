function Fig_1

clear all; close all;
fig1=figure('Units','normalized','Position',[0.05 0.05 0.55 0.9],'Color',[1 1 1]);%[0.03 0.25 0.4 0.45]); %[left bottom width height];   

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time

epsilon = 0.02; 
a       = 0.8; 
b       = 0.03;
 
N =[5 10 15 20 30]./dt; % iterations (time moments) to be visualized  
T_max = 3000; % the last time moment to be visualized  
N_max = T_max/dt; % the last iteration to be visualized
     
kkk = 0; 
%- Initial conditions -
lambda11_0 = zeros(m); lambda11_0(:,ceil(m/2)+1:m) = 0.9;  % U field: lambda11_0*D1+lambda12_0*D2
lambda12_0 = zeros(m); lambda12_0(:,ceil(m/2)+1:m) = 0.9;  % U field 

lambda21_0 = zeros(m); lambda21_0(1:ceil(m/2),:)   = 0.9;  % V field: lambda21_0*D1+lambda22_0*D2
lambda22_0 = zeros(m); lambda22_0(1:ceil(m/2),:)   = 0.9;  % V field


%- Visualization of initial conditions -
kkk = kkk+1;
sub(kkk)=subplot(2,4,kkk); 
imagesc(lambda11_0);   
colormap(bluewhitered(256));    
set(gca,'XTickLabel',[]); 
set(gca,'YTickLabel',[]);

kkk = kkk+1;
sub(kkk)=subplot(2,4,kkk); 
imagesc(lambda21_0);   
colormap(bluewhitered(256));    
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
    
% - Computations -
for i=1:N_max  
    
   lambda11_1 = lambda11_0 + (f(lambda11_0,lambda21_0,epsilon,a,b)+Lattice(lambda11_0,m))*dt; % u:lambda11_0*D1+lambda12_0*D2
   lambda12_1 = lambda12_0 + (f(lambda12_0,lambda22_0,epsilon,a,b)+Lattice(lambda12_0,m))*dt; 
                                 % u:lambda11_0*D1+lambda12_0*D2

   lambda21_1 = lambda21_0 + (g(lambda11_0,lambda21_0) + D*Lattice(lambda21_0,m))*dt; % v
   lambda22_1 = lambda22_0 + (g(lambda12_0,lambda22_0) + D*Lattice(lambda22_0,m))*dt;
                                % v: lambda21_0*D1+lambda22_0*D2

   t = t + dt;
   lambda11_0 = lambda11_1; lambda12_0 = lambda12_1;    
   lambda21_0 = lambda21_1; lambda22_0 = lambda22_1;

% - Vizualization of the dynamics of lambda21_0 - 
figure(fig1);
if  ismember(i,N) || i==N_max
    kkk =kkk+1;    
    sub(kkk)=subplot(2,4,kkk); 
    imagesc(lambda21_0);   % v
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
end

end

%- Unified colorbar -----
Color_Lim = sub(1).CLim;
for iii=2:kkk
    Color_Lim=cat(2,Color_Lim,sub(iii).CLim);
end
Color_max = max(Color_Lim);
Color_min = min(Color_Lim);

for iii=1:kkk
    sub(iii).CLim = [Color_min Color_max];
end

%- Vizualization cosmetics -
%- Subplot positions: fig1 -
pos1 = get(sub(5),'Position'); set(sub(5),'Position',[pos1(1)-0.11 pos1(2)-0.03 pos1(3)-0.01 pos1(4)+0.01]); %[left buttom ]
pos1 = get(sub(1),'Position'); set(sub(1),'Position',[pos1(1)-0.11 pos1(2)-0.06 pos1(3)-0.01 pos1(4)+0.01]); %[left buttom ]
pos1 = get(sub(6),'Position'); set(sub(6),'Position',[pos1(1)-0.14 pos1(2)-0.03 pos1(3)-0.01 pos1(4)+0.01]); %[left buttom ]
pos1 = get(sub(2),'Position'); set(sub(2),'Position',[pos1(1)-0.14 pos1(2)-0.06 pos1(3)-0.01 pos1(4)+0.01]);
pos1 = get(sub(7),'Position'); set(sub(7),'Position',[pos1(1)-0.17 pos1(2)-0.03 pos1(3)-0.01 pos1(4)+0.01]); %[left buttom ]
pos1 = get(sub(3),'Position'); set(sub(3),'Position',[pos1(1)-0.17 pos1(2)-0.06 pos1(3)-0.01 pos1(4)+0.01]);
pos1 = get(sub(4),'Position'); set(sub(4),'Position',[pos1(1)-0.2  pos1(2)-0.06 pos1(3)-0.01 pos1(4)+0.01]);
pos2 = get(sub(8),'Position'); set(sub(8),'Position',[pos2(1)-0.2  pos2(2)-0.03 pos1(3)-0.01 pos1(4)+0.01]);

cbh=colorbar('FontName','Times New Roman','FontSize',22);
cbh.Position=[pos2(1)-0.025  pos2(2)+0.14  0.013 pos2(4)+0.18]; 
cbh.Ticks = linspace(0, 0.9, 10) ; %Create 8 ticks from zero to 1

%- Subplot annotations: fig1 -
annotation(fig1,'textbox',[0.08 0.025 0.055 0.055],'String',{'(e)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08 0.465 0.055 0.055],'String',{'(a)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.255 0.025 0.055 0.055],'String',{'(f)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.255 0.465 0.055 0.055],'String',{'(b)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
       
annotation(fig1,'textbox',[0.43 0.025 0.055 0.055],'String',{'(g)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.43 0.465 0.055 0.055],'String',{'(c)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.61 0.025 0.055 0.055],'String',{'(h)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.61  0.46 0.055 0.055],'String',{'(d)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

%- Functions used in computations -
function [out_data] = Lattice(in_data,in_m)       
 data_temp = Shift(in_data,in_m);
 out_data    = (data_temp.k +  data_temp.d +  data_temp.v +  data_temp.a) - 4* in_data;   
end
 
function [outM] = Shift(inM,in_m)
outM.k = zeros(in_m); % left neighbor 
outM.d = zeros(in_m); % right neighbor
outM.v = zeros(in_m); % upper neighbor
outM.a = zeros(in_m); % bottom neighbor

% no-flux boundary conditions
outM.k = circshift(inM,[ 0  1]); outM.k(:,1)=outM.k(:,2); 
outM.d = circshift(inM,[ 0 -1]); outM.d(:,end)=outM.d(:,end-1); 
outM.v = circshift(inM,[ 1  0]); outM.v(1,:)=outM.v(2,:); 
outM.a = circshift(inM,[-1  0]); outM.a(end,:)=outM.a(end-1,:);  
end 

 % - f(u,v) lambda11 ir  lambda12 -
 function [out_lambda1] = f(in_lambda1,in_lambda2,in_epsilon,in_a,in_b) %  in_x atitinka "u", in_y atitinka "v" (2) formule is M. Vaidelys
 out_lambda1 = (1/in_epsilon)*in_lambda1.*(1-in_lambda1).*(in_lambda1-(1/in_a)*(in_lambda2+in_b));
 end 
 
% -g(u,v) lambda21 ir lambda22 - 
 function [out_data] = g(in_data1,in_data2) 
 out_data = in_data1-in_data2;
 end 

 end




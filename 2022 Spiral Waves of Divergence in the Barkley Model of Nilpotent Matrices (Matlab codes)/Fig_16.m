function Fig_16

clear all;
close all;

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time

epsilon = 0.02; 
a       = 0.8; 
b       = 0.005; 

N =[20 20.5]./dt; % iterations to be visualized for a spiral wave 
T_max = 21; % the last time moment to be visualized  
N_max = T_max/dt; % the last iteration to be visualized


fig1=figure('Units','normalized','Position',[0.2 0.1 0.41 0.9],'Color',[1 1 1]);
kkk = 0; 
  
%- Initial conditions -
lambda1_0 = zeros(m); lambda1_0(:,ceil(m/2)+1:m) = 0.9; % U field
lambda2_0 = zeros(m); lambda2_0(1:ceil(m/2),:)   = 0.9; % V field
mu1_0 = ones(m); mu2_0 = ones(m);  % U and V field accordingly

% - Computations -  
for i=1:N_max 
    
   mu1_1 = mu1_0 + (f_mu1(mu1_0,mu2_0,lambda1_0,lambda2_0,epsilon,a,b)+Lattice(mu1_0,m))*dt; % v
   mu2_1 = mu2_0 + (g(mu1_0,mu2_0) + D*Lattice(mu2_0,m))*dt;  % v

   lambda1_1 = lambda1_0 + (f_lambda1(lambda1_0,lambda2_0,epsilon,a,b)+Lattice(lambda1_0,m))*dt; % u
   lambda2_1 = lambda2_0 + (g(lambda1_0,lambda2_0) + D* Lattice(lambda2_0,m))*dt;  %u
   
   t = t + dt;
   mu1_0 = mu1_1; mu2_0 = mu2_1;   
   lambda1_0 =lambda1_1; lambda2_0 = lambda2_1; 

% - Vizualization of the dynamics of mu2 - 
if  ismember(i,N) || i==N_max
    figure(fig1);    
    kkk =kkk+1;    
    sub(kkk)=subplot(2,3,kkk); 
    imagesc(mu2_0);   

     % - Symetric colorbar in each subplot -
     Color_Lim = sub(kkk).CLim; 
     max_abs_color = max(abs(Color_Lim));
     sub(kkk).CLim = [-max_abs_color max_abs_color];

     colormap(bluewhitered(256));     
     colorbar('FontName','Times New Roman','FontSize',22); 
     set(gca,'XTickLabel',[]);
     set(gca,'YTickLabel',[]);
end


end

%- Moire gratings -
raiska = 35; 
h_x = 2*pi/(raiska-1);
h_y = 2*pi/(raiska-1);

[xx,yy] = meshgrid(0:h_x:30*pi,0:h_y:30*pi);
lambda = 10; 
k = 2*pi/lambda;

theta_laipsniais = 20; 
theta = theta_laipsniais*2*pi/360;

M1 = (1/2)*(1+cos(k*xx));
xx_rot = xx*cos(theta) - yy*sin(theta); 
M2 = (1/2)*(1+cos(k*xx_rot));

M_superpozicija = (M1+M2)/2;

kkk=kkk+1;
sub(kkk)=subplot(2,3,kkk); 
imagesc(M1);   
colormap(sub(kkk),'gray'); 
cbh=colorbar('FontName','Times New Roman','FontSize',22); 
set(cbh,'XTick',[0.0001 0.2 0.4 0.6 0.8 1]);
set(cbh,'XTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

kkk=kkk+1;
sub(kkk)=subplot(2,3,kkk); 
imagesc(M2);   
colormap(sub(kkk),'gray'); 
cbh=colorbar('FontName','Times New Roman','FontSize',22);

set(cbh,'XTick',[0.0001 0.2 0.4 0.6 0.8 1]);
set(cbh,'XTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

kkk=kkk+1;
sub(kkk)=subplot(2,3,kkk); 
imagesc(M_superpozicija);   
colormap(sub(kkk),'gray'); 
cbh=colorbar('FontName','Times New Roman','FontSize',22);

set(cbh,'XTick',[0.0001 0.2 0.4 0.6 0.8 1]);
set(cbh,'XTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

%- Vizualization cosmetics -
%- Subplot positions: fig1 -
pos1 = get(sub(4),'Position'); set(sub(4),'Position',[pos1(1)-0.11 pos1(2)-0.02 pos1(3) pos1(4)]); %[left buttom ]
pos2 = get(sub(5),'Position'); set(sub(5),'Position',[pos2(1)-0.14 pos2(2)-0.02 pos1(3) pos2(4)]); %[left buttom ]
pos3 = get(sub(6),'Position'); set(sub(6),'Position',[pos3(1)-0.17 pos3(2)-0.02 pos1(3) pos3(4)]); %[left buttom ]
pos1 = get(sub(1),'Position'); set(sub(1),'Position',[pos1(1)-0.11 pos1(2)-0.06 pos1(3) pos1(4)]);
pos2 = get(sub(2),'Position'); set(sub(2),'Position',[pos2(1)-0.14 pos2(2)-0.06 pos1(3) pos1(4)]); %[left buttom ]
pos2 = get(sub(3),'Position'); set(sub(3),'Position',[pos2(1)-0.17 pos2(2)-0.06 pos1(3) pos2(4)]);


%- Subplot annotations: fig1 -
 annotation(fig1,'textbox',[0.095 0.46 0.055 0.055],'String',{'(a)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
 
 annotation(fig1,'textbox',[0.344 0.46 0.055 0.055],'String',{'(b)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
       
 annotation(fig1,'textbox',[0.6 0.46 0.055 0.055],'String',{'(c)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

 annotation(fig1,'textbox',[0.095 0.028 0.055 0.055],'String',{'(d)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
 
 annotation(fig1,'textbox',[0.344 0.028 0.055 0.055],'String',{'(e)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
       
 annotation(fig1,'textbox',[0.6 0.028 0.055 0.055],'String',{'(f)'},'FontName',...
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

%- f(u,v) mu1 -
 function [out_mu1] = f_mu1(in_mu1,in_mu2,in_lambda1,in_lambda2,in_epsilon,in_a,in_b) %  in_x atitinka "u", in_y atitinka "v" (2) formule is M. Vaidelys
 out_mu1 = (1/in_epsilon)*((1-2*in_lambda1).*(in_lambda1-(1/in_a)*(in_lambda2+in_b)).*in_mu1+in_lambda1.*(1-in_lambda1).*(in_mu1-(1/in_a)*in_mu2));
 end 
 
  %- f(u,v) lambda1 - 
 function [out_lambda1] = f_lambda1(in_lambda1,in_lambda2,in_epsilon,in_a,in_b) 
 out_lambda1 = (1/in_epsilon)*in_lambda1.*(1-in_lambda1).*(in_lambda1-(1/in_a)*(in_lambda2+in_b));
 end 
 
  %- g(u,v) lambda2 and mu2 - 
 function [out_data] = g(in_data1,in_data2) 
 out_data = in_data1-in_data2;
 end 

 end




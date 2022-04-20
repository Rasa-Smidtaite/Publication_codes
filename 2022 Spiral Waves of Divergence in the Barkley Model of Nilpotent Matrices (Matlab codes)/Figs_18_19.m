function Figs_18_19

clear all;
close all;
fig2=figure('Units','normalized','Position',[0.03 0.25 0.24 0.4],'Color',[1 1 1]); 
fig1=figure('Units','normalized','Position',[0.05 0.05 0.55 0.9],'Color',[1 1 1]);

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time

epsilon = 0.02; 
a       = 0.8; 
b       = 0.03;

N =[5 10 15 20 50 100 300]./dt; % iterations (time moments) to be visualized  
T_max = 3000; % the last time moment to be visualized  
N_max = T_max/dt; % the last iteration to be visualized

kkk = 0; 
mu1_0 = ones(m);  mu2_0 = ones(m);   % mu1 and mu2 initial conditions

%- Reproduce the exact initial conditions use rng('Seed_Figs_18_19.mat') -
Data = load('Seed_Figs_18_19.mat'); 
rng(Data.seed); 
 
%-- lambda01 initial conditions (with a diffusion type distribution) --
lambda1_0 = zeros(m); 
noise_begin = round(0.40*m); 
kk = 2; % the number of colors
 
for i = 1:kk-1 aa(i) = 10*i/(kk*(m-1)); bb(i) = -aa(i);  end 
for ix = 1:m
    for iy = noise_begin:m
    for i=1:kk-1 p(i) = aa(i)*(iy-round(noise_begin)+1) + bb(i);  end
     s = rand;
     if s<p(1); lambda1_0(ix,iy) = 0.9; end
    end 
end

%-- lambda02 initial conditions (with a diffusion type distribution) --
 lambda2_0 = zeros(m); 
 for ix = 1:m
    for iy = noise_begin:m
    for i=1:kk-1 p(i) = aa(i)*(iy-round(noise_begin)+1) + bb(i);  end
     s = rand;
     if s<p(1); lambda2_0(iy,ix) = 0.9;      end
    end 
end
lambda2_0 = flip(lambda2_0);

% - Vizualization of initial conditions: fig2 -
figure(fig2);  
sub2(1)=subplot(1,2,1); 
imagesc(lambda1_0);    
colormap(bluewhitered(256));     
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

sub2(2)=subplot(1,2,2); 
imagesc(lambda2_0);   
colormap(bluewhitered(256));     
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]); 

%- Subplot annotations: fig2 -
annotation(fig2,'textbox',[0.15 0.1 0.055 0.055],'String',{'(a)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
annotation(fig2,'textbox',[0.14+0.4  0.1 0.055 0.055],'String',{'(b)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

%- Subplot positions: fig2 -
pos1 = get(sub2(1),'Position'); set(sub2(1),'Position',[pos1(1)-0.1 pos1(2)+0.06 pos1(3)-0.01 pos1(4)-0.11]); %[left buttom ]
pos2 = get(sub2(2),'Position'); set(sub2(2),'Position',[pos2(1)-0.17 pos2(2)+0.06 pos2(3)-0.01 pos2(4)-0.11]); %[left buttom ]

%- Unified colorbar and colorbar position: fig2 -
Color_Lim = sub2(1).CLim;
for iii=2:2 Color_Lim=cat(2,Color_Lim,sub2(iii).CLim); end
Color_max = max(Color_Lim); Color_min = min(Color_Lim);
for iii=1:2 sub2(iii).CLim = [Color_min Color_max]; end

cbh=colorbar('FontName','Times New Roman','FontSize',20,'XTick', 0:0.1:0.9);
cbh.Position=[pos2(1)+0.22  pos2(2)+0.06  0.019 pos2(4)-0.11]; 

% - Computations -
for i=1:N_max 
    
   mu1_1 = mu1_0 + (f_mu1(mu1_0,mu2_0,lambda1_0,lambda2_0,epsilon,a,b)+Lattice(mu1_0,m))*dt; % v
   mu2_1 = mu2_0 + (g(mu1_0,mu2_0) + D*Lattice(mu2_0,m))*dt;  % v

   lambda1_1 = lambda1_0 + (f_lambda1(lambda1_0,lambda2_0,epsilon,a,b)+Lattice(lambda1_0,m))*dt; % u
   lambda2_1 = lambda2_0 + (g(lambda1_0,lambda2_0) + D* Lattice(lambda2_0,m))*dt;  %u
  
   t = t + dt;
   mu1_0 = mu1_1; mu2_0 = mu2_1;   
   lambda1_0 =lambda1_1; lambda2_0 = lambda2_1; 

Max(i) = max(max(abs(mu2_0)));

% - Vizualization of the dynamics of mu2 - 
figure(fig1);
if  ismember(i,N) || i==N_max
    kkk =kkk+1;    
    sub1(kkk)=subplot(2,4,kkk); 
    imagesc(mu2_0);   

    % - Symetric colorbar in each subplot -
      Color_Lim = sub1(kkk).CLim;
      max_abs_color = max(abs(Color_Lim));
      sub1(kkk).CLim = [-max_abs_color max_abs_color];

      colormap(bluewhitered(256));     
      colorbar('FontName','Times New Roman','FontSize',20);
      set(gca,'XTickLabel',[]);
      set(gca,'YTickLabel',[]);
end

end

%- Vizualization cosmetics -
%- Subplot positions: fig1 -
pos1 = get(sub1(5),'Position'); set(sub1(5),'Position',[pos1(1)-0.11 pos1(2)-0.03 pos1(3) pos1(4)]);   
pos1 = get(sub1(1),'Position'); set(sub1(1),'Position',[pos1(1)-0.11 pos1(2)-0.07 pos1(3) pos1(4)]);  
pos1 = get(sub1(6),'Position'); set(sub1(6),'Position',[pos1(1)-0.14 pos1(2)-0.03 pos1(3) pos1(4)]);  
pos1 = get(sub1(2),'Position'); set(sub1(2),'Position',[pos1(1)-0.14 pos1(2)-0.07 pos1(3) pos1(4)]);
pos1 = get(sub1(7),'Position'); set(sub1(7),'Position',[pos1(1)-0.17 pos1(2)-0.03 pos1(3) pos1(4)]);  
pos1 = get(sub1(3),'Position'); set(sub1(3),'Position',[pos1(1)-0.17 pos1(2)-0.07 pos1(3) pos1(4)]);
pos1 = get(sub1(4),'Position'); set(sub1(4),'Position',[pos1(1)-0.2  pos1(2)-0.07 pos1(3) pos1(4)]);
pos2 = get(sub1(8),'Position'); set(sub1(8),'Position',[pos2(1)-0.2  pos2(2)-0.03 pos1(3) pos1(4)]);
clear pos1 pos2;

%- Subplot annotations: fig1 -
annotation(fig1,'textbox',[0.07 0.015 0.055 0.055],'String',{'(e)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.07 0.45 0.055 0.055],'String',{'(a)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.07+0.18 0.015 0.055 0.055],'String',{'(f)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.07+0.18 0.45 0.055 0.055],'String',{'(b)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
       
annotation(fig1,'textbox',[0.07+0.18+0.172 0.015 0.055 0.055],'String',{'(g)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.07+0.18+0.172 0.45 0.055 0.055],'String',{'(c)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08+0.18+0.172+0.172 0.015 0.055 0.055],'String',{'(h)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08+0.18+0.172+0.172 0.45 0.055 0.055],'String',{'(d)'},'FontName',...
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




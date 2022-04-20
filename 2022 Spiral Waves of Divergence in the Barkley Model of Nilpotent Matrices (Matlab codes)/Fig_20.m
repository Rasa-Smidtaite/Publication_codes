function  Fig_20

clear all;
close all;
fig1=figure('Units','normalized','Position',[0.03 0.25 0.5 0.45],'Color',[1 1 1]);

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time

epsilon = 0.02; 
a       = 0.8; 
b       = 0.03;

N =[5 10 15]./dt; % iterations (time moments) to be visualized  
T_max = 20; % the last time moment to be visualized  
N_max = T_max/dt; % the last iteration to be visualized

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

Max(i) = max(max(abs(mu2_0)));

% - Perturbation is added randomly max or min at T=10 for nodes  45:55 and 70:80 -
rng(3); % to reproduce the exact noise use rng function
a1 = 45; b1 = 55; a2 = 75;  b2 = 85; tasku_kiekis = 20;
r1 = randi([a1 b1],tasku_kiekis,1);
r2 = randi([a2 b2],tasku_kiekis,1);
A = [1,-1]; r3 = A(randi([1,2],tasku_kiekis,1)); r3=r3';
if i==1000
 for ii=1:tasku_kiekis; mu2_0(r1(ii),r2(ii)) = r3(ii)*Max(i); end    
end 

% - Vizualization of the dynamics of mu2 - 
figure(fig1);
if  ismember(i,N) || i==N_max
    kkk =kkk+1;    
    sub(kkk)=subplot(1,4,kkk); 
    imagesc(mu2_0);    

    % - Symetric colorbar in each subplot -
    Color_Lim = sub(kkk).CLim; % paimam colorbar ribas
    max_abs_color = max(abs(Color_Lim));
    sub(kkk).CLim = [-max_abs_color max_abs_color];

    colormap(bluewhitered(256));     
    colorbar('FontName','Times New Roman','FontSize',20);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end

end
 
%- Vizualization cosmetics -
%- Subplot positions: fig1 -
 pos1 = get(sub(1),'Position'); set(sub(1),'Position',[pos1(1)-0.11 pos1(2)+0.04 pos1(3)+0.02 pos1(4)-0.18]); %[left buttom ]
 pos2 = get(sub(2),'Position'); set(sub(2),'Position',[pos2(1)-0.11 pos2(2)+0.04 pos2(3)+0.02 pos2(4)-0.18]); %[left buttom ]
 pos3 = get(sub(3),'Position'); set(sub(3),'Position',[pos3(1)-0.11 pos3(2)+0.04 pos3(3)+0.02 pos3(4)-0.18]); %[left buttom ]
 pos4 = get(sub(4),'Position'); set(sub(4),'Position',[pos4(1)-0.11 pos4(2)+0.04 pos3(3)+0.02 pos3(4)-0.18]); %[left buttom ]
 clear pos1 pos2 pos3 pos4;

%- Subplot annotations -
annotation(fig1,'textbox',[0.08 0.1 0.055 0.055],'String',{'(a)'},'FontName',...
           'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08+0.205 0.1 0.055 0.055],'String',{'(b)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08+0.205+0.205 0.1 0.055 0.055],'String',{'(c)'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');

annotation(fig1,'textbox',[0.08+0.205+0.205+0.205 0.1 0.055 0.055],'String',{'(d)'},'FontName',...
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
 function [out_mu1] = f_mu1(in_mu1,in_mu2,in_lambda1,in_lambda2,in_epsilon,in_a,in_b) 
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




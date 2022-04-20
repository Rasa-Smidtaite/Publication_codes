function Fig_8

clear all; close all;
fig3=figure('Units','normalized','Position',[0.5 0.1 0.25 0.5],'Color',[1 1 1]);%[0.03 0.25 0.4 0.45]); %[left bottom width height];   

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time

epsilon = 0.02; 
a       = 0.8; 
b       = 0.005;
   
T = 100; % the last time moment to be visualized  
N = T/dt; % the last iteration to be visualized

%- Initial conditions -
lambda1_0 = zeros(m); lambda1_0(:,ceil(m/2)+1:m) = 0.9; % U field
lambda2_0 = zeros(m); lambda2_0(1:ceil(m/2),:)   = 0.9; % V field
mu1_0 = ones(m); mu2_0 = ones(m);  % U and V field accordingly

% - Computations -
for i=1:N
    
   mu1_1 = mu1_0 + (f_mu1(mu1_0,mu2_0,lambda1_0,lambda2_0,epsilon,a,b)+Lattice(mu1_0,m))*dt; % v
   mu2_1 = mu2_0 + g(mu1_0,mu2_0)*dt;  % v

   lambda1_1 = lambda1_0 + (f_lambda1(lambda1_0,lambda2_0,epsilon,a,b)+Lattice(lambda1_0,m))*dt; % u
   lambda2_1 = lambda2_0 + g(lambda1_0,lambda2_0)*dt;  % u

   t = t + dt;
   mu1_0 = mu1_1; mu2_0 = mu2_1;   
   lambda1_0 =lambda1_1; lambda2_0 = lambda2_1; 

Max(i) = max(max(abs(mu2_0)));

end

% - Vizualization of the dynamics of mu2 - 
k = 1;
Mu2_max_vaizdavimui(1) = 1; 
for k=2:(1+T/(20*dt)) 
    Mu2_max_vaizdavimui(k) = max(Max(20*(k-2)+1:20*(k-2)+20)); 
end
[maximum,i_max]= max(Mu2_max_vaizdavimui); 
maximum = round(maximum,2);

set(gca, 'TickLabelInterpreter', 'latex');
xticks_array = [0+1 25+1 50+1  75+1 100+1  250+1  500];
xticks_label_array = {'0' '5' '10' '15' '20' '50' '100'};

sub1_1 =subplot(1,1,1); 
yyaxis left
plot(Mu2_max_vaizdavimui,'LineWidth',1.1,'LineStyle','-' ,'Color',[0 0 0]);
ax = gca;
ax.YColor = [0 0 0];
y_lim = ylim;  y_lim(2) = maximum;
hold on; plot([i_max i_max],y_lim,'LineStyle','--' ,'Color',[0 0 0])

xlim([0 k]);
x_lim = xlim;
xticks(xticks_array);
xticklabels(xticks_label_array);

%- Logarithmic scale -
yyaxis right
plot(sign(Mu2_max_vaizdavimui).*(log10(1+abs(Mu2_max_vaizdavimui))),'LineWidth',1.1,'LineStyle','-' ,'Color',[0 0 0.6]);
ax = gca;
ax.YColor = [0 0 0.6];

annotation(fig3,'textbox',[0.01 0.9 0.22 0.055],'String',{'max |{\it\mu}_2^{(\it{t})}|'},'FontName',...
           'Times New Roman','FontSize',22,'FitBoxToText','off','LineStyle','none'); %[left bottom width height];

annotation(fig3,'textbox',[0.86 0.03 0.07 0.07],'String',{'\it T  '},'FontName',...
           'Times New Roman','FontSize',22,'FitBoxToText','off','LineStyle','none'); %[left bottom width height];

annotation(fig3,'textbox',[0.71 0.9 0.25 0.055],'String',{'log_{10} (1+max |{\it\mu}^{(\it{t})}_2|)'},'FontName',...
           'Times New Roman','FontSize',22,'Color',[0 0 0.6],'FitBoxToText','off','LineStyle','none');

set(gca,'FontName','Times New Roman','FontSize',19);     
pos1 = get(sub1_1,'Position'); set(sub1_1,'Position',[pos1(1)-0.06 pos1(2)-0.02 pos1(3) pos1(4)-0.1]); %[left buttom width height]

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




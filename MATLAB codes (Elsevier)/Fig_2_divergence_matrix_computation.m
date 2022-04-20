function Divergence_matrix=Fig_2_divergence_matrix_computation

m = 100; % domain size m x m
dt = 0.01; %time step - note as h in the article

D = 0; % ratio of diffusion coeffcients
t = 0; % initial time
epsilon = 0.02; 

h_a = 0.001;
h_b = 0.0005;
a0 = 0.2; a1 = 1; % parameter a limits
b0 = 0  ; b1 = 0.14; % parameter b limits
 
T =100;
N = T/dt; 

Divergence_matrix = (-1)*ones(1+floor((a1-a0)/h_a),1+floor((b1-b0)/h_b));  % Divergence array 
 
i=0;
for a = a0:h_a:a1
    i=i+1;
    disp(i);
    j=0;
    for b = b0:h_b:b1
        j=j+1;
        mu1_0 = ones(m);  mu2_0 = ones(m);  
                
       lambda1_0 = zeros(m); lambda1_0(:,ceil(m/2)+1:m) = 0.9;
       lambda2_0 = zeros(m); lambda2_0(1:ceil(m/2),:)   = 0.9;
           
  for n = 1:N
       
   mu1_1 = mu1_0 + (f_mu1(mu1_0,mu2_0,lambda1_0,lambda2_0,epsilon,a,b)+Lattice(mu1_0,m))*dt;
   mu2_1 = mu2_0 + (g(mu1_0,mu2_0) + D* Lattice(mu2_0,m))*dt;
 
   lambda1_1 = lambda1_0 + (f_lambda1(lambda1_0,lambda2_0,epsilon,a,b)+Lattice(lambda1_0,m))*dt;
   lambda2_1 = lambda2_0 + (g(lambda1_0,lambda2_0) + D* Lattice(lambda2_0,m))*dt;
 
 
t = t + dt;
mu1_0 = mu1_1; mu2_0 = mu2_1;   
lambda1_0 =lambda1_1; lambda2_0 = lambda2_1; 
 
Max(n+1) = max(max(abs(mu2_0)));

end
 
maximum = max(Max);
max_mu2_0 = max(max(abs(mu2_0)));
if      ((max_mu2_0>=1000)==1) || (sum(sum(isnan(mu2_0)))>=1) 
        Divergence_matrix(i,j) = 1; % type (i)
elseif  ((max_mu2_0 < 1000) == 1) 
      if  ((maximum>=1000)==1)
           Divergence_matrix(i,j) = 0.75; % type (ii)
      elseif  ((maximum>=100)==1) && ((maximum<1000)==1) 
           Divergence_matrix(i,j) = 0.5; % type (iii)
      elseif ((maximum>=10) ==1) && ((maximum<100)==1)   
          Divergence_matrix(i,j) = 0.25; % type (iv)
      elseif ((maximum>=0)==1) && ((maximum<10)==1)  
          Divergence_matrix(i,j) = 0; % type (v)
      end
end
 
end
end

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

function [P_e,eff,emissivity] = BUFED_TPV_main(t1,t2,t3)

load('Ge_revised.mat')
load('InAs_revised.mat')

h_bar = 1.054571817e-34;        %reduced Planck's constant
k_b = 1.3807e-23;               %Boltzmann's constant
c = 2.99792e8;                  %speed of light

l = 5;
d = [0 0 t1 t1+t2 t1+t2+t3];
T = 1100;

gamma = 273e11;
w_p = 13.69e15;

eps0 = zeros(1,length(W));
eps1 = zeros(1,length(W));
eps2 = zeros(1,length(W));
eps3 = zeros(1,length(W));
eps4 = zeros(1,length(W));

for i=1:length(W)
   eps0(i) = 1 - w_p^2/(W(i)*(W(i)+1i*gamma));
   eps1(i) = Eps_Ge_r(i) + 1i*Eps_Ge_i(i);
   eps2(i) = eps_Si(W(i),T);
   eps3(i) = Eps_Ge_r(i) + 1i*Eps_Ge_i(i);
   eps4(i) = 1;
end

Eps(:,1) = W;
Eps(:,2) = eps0;
Eps(:,3) = eps1;
Eps(:,4) = eps2;
Eps(:,5) = eps3;
Eps(:,6) = eps4;

q_w = BUFED_v5(l,1,d,T,Eps) + BUFED_v5(l,2,d,T,Eps) + BUFED_v5(l,3,d,T,Eps) + BUFED_v5(l,4,d,T,Eps);
theta = h_bar*W./(exp(h_bar.*W/(k_b.*T))-1);
q_b_w = theta.*W.^2/(4*pi^2*c^2);
emissivity = q_w./q_b_w;
[P_e,eff,~] = InAs_TPV(T,emissivity,InAs_IQE);

end

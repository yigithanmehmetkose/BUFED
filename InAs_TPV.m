function [P_e,eff,e_tot] = InAs_TPV(T_e,E_e,InAs_IQE)

%Parameters:
h_bar = 1.054571817e-34;    % reduced Planck constant, J-s
k_b = 1.380649e-23;         % Boltzmann constant, J/K
sigma = 5.670374419e-8;     % Stefan-Boltzmann constant, W/m^2-K^4
q = 1.60217663e-19;         % electron charge, C
c = 2.99792458e8;           % speed of light, m/s

%InAs Material Properties
E_g = 0.354*q;              % bandgap, J
w_g = E_g/h_bar;            % bandgap frequency, rad/s
N_i = 1e15;                 % intrinsic carrier concentration, cm^-3
N_a = 1e18;                 % acceptor concentration, cm^-3
D_e = 1e3;                  % electron diffusion coefficient, cm^2/s
tau_e = 30e-9;              % minority electron lifetime, s
N_d = 1e19;                 % donor concentration, cm^-3
D_h = 1.2;                  % hole diffusion coefficient, cm^2/s
tau_h = 3e-6;               % minority hole lifetime, s

%Inputs
T_c = 300;                  % cell temperature, K
N = 1083;                   % number of integration steps

%Preprocessing
dw = 1e12;           % integration step size
W = dw*(1:3*N+1);           % angular frequency
I_j = zeros(1,N);
I_p = zeros(1,N);

%Integration
for p=1:N
    %Preprocessing
    w = [W(3*p-2) W(3*p-1) W(3*p) W(3*p+1)];
    IQE = [InAs_IQE(3*p-2) InAs_IQE(3*p-1) InAs_IQE(3*p) InAs_IQE(3*p+1)];
    e = [E_e(3*p-2) E_e(3*p-1) E_e(3*p) E_e(3*p+1)];
    P_w = zeros(1,4);
    N_w = zeros(1,4);
    %Simpson's 3/8 Rule
    for i=1:4
        P_w(i) = e(i)*h_bar*w(i)^3/(4*pi^2*c^2*(exp(h_bar*w(i)/(k_b*T_e))-1));
        N_w(i) = IQE(i)*e(i)*w(i)^2/(4*pi^2*c^2*(exp(h_bar*w(i)/(k_b*T_e))-1));
    end
    I_p(p) = 3*dw/8*(P_w(1)+3*P_w(2)+3*P_w(3)+P_w(4));
    I_j(p) = 3*dw/8*(N_w(1)+3*N_w(2)+3*N_w(3)+N_w(4));
end

%Postprocessing
P_r = sum(I_p);                                                             % radiative power, W/m^2
J_ph = q*sum(I_j);                                                          % photocurrent density, A/m^2
J_s = q*N_i^2*(sqrt(D_e/tau_e)/N_a + sqrt(D_h/tau_h)/N_d);                  % saturation current density, A/m^2
V_oc = k_b*T_c/q*log(J_ph/J_s + 1);                                         % open circuit voltage, V
P_e = J_ph*V_oc*(1-1/log(J_ph/J_s))*(1-log(log(J_ph/J_s))/log(J_ph/J_s));   % electrical power, W/m^2
eff = P_e/P_r;                                                              % conversion efficiency
e_tot = P_r/(sigma*T_e^4);                                                  % emitter total emissivity

end
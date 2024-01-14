function epsilon = eps_Si(w,T)

h_bar = 1.054571817e-34;        %reduced Planck's constant
q = 1.60217663e-19;             %electron charge

E_h = w*h_bar/q;
T0 = 293;
Tdel = (T-T0)/T0;

a_0 = [4.87e-3 0.7722];
a_1 = [-8.936e-4 5.984e-3];
a_2 = [7.854e-4 5.586e-4];

omega_0 = [0.1289 0.3129];
omega_1 = [-4.571e-3 6.405e-4];
omega_2 = [9.421e-4 -6.527e-4];

gamma_0 = [1.875e-2 9.742e-2];
gamma_1 = [9.274e-4 4.814e-2];
gamma_2 = [1.651e-3 -1.518e-2];

gammaprime_0 = [0.1387 9.505e-2];
gammaprime_1 = [1.161e-2 5.607e-2];
gammaprime_2 = [-1.543e-2 -1.948e-2];

a = a_0 + a_1*Tdel + a_2*Tdel^2;
omega = omega_0 + omega_1*Tdel + omega_2*Tdel^2;
gamma = gamma_0 + gamma_1*Tdel + gamma_2*Tdel^2;
gammaprime = gammaprime_0 + gammaprime_1*Tdel + gammaprime_2*Tdel^2;

F = @(w) a.*(omega.^2-1i*gammaprime*w)./(omega.^2-w^2-1i*w*gamma);
f = sum(F(E_h));
epsilon = (1+2*f)/(1-f);

end


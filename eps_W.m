function eps = eps_W(w)

h_bar = 1.054571817e-34;        %reduced Planck's constant
q = 1.60217663e-19;             %electron charge

E = w*h_bar/q;

omega_p = 13.22;

f0 = 0.206;
gamma0 = 0.064;

f1 = 0.054;
gamma1 = 0.530;
omega1 = 1.004;

f2 = 0.166;
gamma2 = 1.281;
omega2 = 1.917;

f3 = 0.706;
gamma3 = 3.332;
omega3 = 3.580;

f4 = 2.590;
gamma4 = 5.836;
omega4 = 7.498;

Omega_p = sqrt(f0)*omega_p;

eps1 = 1-Omega_p^2/(E*(E+1i*gamma0));
eps2 =  f1*omega_p^2 / ((omega1^2-E^2)-1i*E*gamma1);
eps3 = f2*omega_p^2 / ((omega2^2-E^2)-1i*E*gamma2);
eps4 = f3*omega_p^2 / ((omega3^2-E^2)-1i*E*gamma3);
eps5 = f4*omega_p^2 / ((omega4^2-E^2)-1i*E*gamma4);

eps = eps1 + eps2 + eps3 + eps4 + eps5;

end
function [A_TE,B_TE,C_TE,D_TE,A_TM,B_TM,C_TM,D_TM] = S_matrix_method(w,kp,Eps,Z,s,l)

%% Inputs:
c = 2.99792e8;                  %speed of light
N = length(Z);                  %number of layers except half spaces

%% Initialization:
s11_TE_0_l = zeros(1,N);
s12_TE_0_l = zeros(1,N);
s21_TE_0_l = zeros(1,N);
s22_TE_0_l = zeros(1,N);
s11_TM_0_l = zeros(1,N);
s12_TM_0_l = zeros(1,N);
s21_TM_0_l = zeros(1,N);
s22_TM_0_l = zeros(1,N);

s11_TE_s_n = zeros(1,N);
s12_TE_s_n = zeros(1,N);
s21_TE_s_n = zeros(1,N);
s22_TE_s_n = zeros(1,N);
s11_TM_s_n = zeros(1,N);
s12_TM_s_n = zeros(1,N);
s21_TM_s_n = zeros(1,N);
s22_TM_s_n = zeros(1,N);

s11_TE_0_l(1) = 1;
s12_TE_0_l(1) = 0;
s21_TE_0_l(1) = 0;
s22_TE_0_l(1) = 1;

s11_TM_0_l(1) = 1;
s12_TM_0_l(1) = 0;
s21_TM_0_l(1) = 0;
s22_TM_0_l(1) = 1;

s11_TE_s_n(s) = 1;
s12_TE_s_n(s) = 0;
s21_TE_s_n(s) = 0;
s22_TE_s_n(s) = 1;

s11_TM_s_n(s) = 1;
s12_TM_s_n(s) = 0;
s21_TM_s_n(s) = 0;
s22_TM_s_n(s) = 1;

%% S-matrix between layer 0 and l:

for i=1:N-1
    eps1 = Eps(i);
    eps2 = Eps(i+1);
    kz1 = sqrt(eps1*(w/c)^2-kp^2);                      %z-components of wavevectors
    kz2 = sqrt(eps2*(w/c)^2-kp^2);

    r_s = (kz1-kz2)/(kz1+kz2);                          %reflection coefficient in s polarization
    t_s = 1+ (kz1-kz2)/(kz1+kz2);                       %tra+nsmission coefficient in s polarization
    r_p = (eps2*kz1-eps1*kz2)/(eps2*kz1+eps1*kz2);      %reflection coefficient in p polarization
    t_p = 2*sqrt(eps1*eps2)*kz1/(eps2*kz1+eps1*kz2);    %transmission coefficient in p polarization    

    s11_TE_0_l(i+1) = s11_TE_0_l(i)*t_s*exp(1i*kz1*(Z(i+1)-Z(i)))/(1-s12_TE_0_l(i)*r_s*exp(2i*kz1*(Z(i+1)-Z(i))));
    s12_TE_0_l(i+1) = (s12_TE_0_l(i)*exp(2i*kz1*(Z(i+1)-Z(i)))-r_s)/(1-s12_TE_0_l(i)*r_s*exp(2i*kz1*(Z(i+1)-Z(i))));
    s21_TE_0_l(i+1) = s11_TE_0_l(i+1)*s22_TE_0_l(i)*r_s*exp(1i*kz1*(Z(i+1)-Z(i)))/t_s + s21_TE_0_l(i);
    s22_TE_0_l(i+1) = s22_TE_0_l(i)*(r_s*s12_TE_0_l(i+1)+1)*exp(1i*kz1*(Z(i+1)-Z(i)))/t_s;

    s11_TM_0_l(i+1) = s11_TM_0_l(i)*t_p*exp(1i*kz1*(Z(i+1)-Z(i)))/(1-s12_TM_0_l(i)*r_p*exp(2i*kz1*(Z(i+1)-Z(i))));
    s12_TM_0_l(i+1) = (s12_TM_0_l(i)*exp(2*1i*kz1*(Z(i+1)-Z(i)))-r_p)/(1-s12_TM_0_l(i)*r_p*exp(2i*kz1*(Z(i+1)-Z(i))));
    s21_TM_0_l(i+1) = (s11_TM_0_l(i+1)*s22_TM_0_l(i)*r_p*exp(1i*kz1*(Z(i+1)-Z(i))))/t_p + s21_TM_0_l(i);
    s22_TM_0_l(i+1) = s22_TM_0_l(i)*(r_p*s12_TM_0_l(i+1)+1)*exp(1i*kz1*(Z(i+1)-Z(i)))/t_p;
end

%% S-matrix between layer s and N:

for i=s:N-1
    eps1 = Eps(i);
    eps2 = Eps(i+1);    
    kz1 = sqrt(eps1*(w/c)^2-kp^2);
    kz2 = sqrt(eps2*(w/c)^2-kp^2);    
    
    r_s = (kz1-kz2)/(kz1+kz2);
    t_s = 1+ (kz1-kz2)/(kz1+kz2);
    r_p = (eps2*kz1-eps1*kz2)/(eps2*kz1+eps1*kz2);
    t_p = 2*sqrt(eps1*eps2)*kz1/(eps2*kz1+eps1*kz2);    

    s11_TE_s_n(i+1) = s11_TE_s_n(i)*t_s*exp(1i*kz1*(Z(i+1)-Z(i)))/(1-s12_TE_s_n(i)*r_s*exp(2*1i*kz1*(Z(i+1)-Z(i))));
    s12_TE_s_n(i+1) = (s12_TE_s_n(i)*exp(2*1i*kz1*(Z(i+1)-Z(i)))-r_s)/(1-s12_TE_s_n(i)*r_s*exp(2*1i*kz1*(Z(i+1)-Z(i))));
    s21_TE_s_n(i+1) = (s11_TE_s_n(i+1)*s22_TE_s_n(i)*r_s*exp(1i*kz1*(Z(i+1)-Z(i))))/t_s + s21_TE_s_n(i);
    s22_TE_s_n(i+1) = s22_TE_s_n(i)*(r_s*s12_TE_s_n(i+1)+1)*exp(1i*kz1*(Z(i+1)-Z(i)))/t_s;

    s11_TM_s_n(i+1) = s11_TM_s_n(i)*t_p*exp(1i*kz1*(Z(i+1)-Z(i)))/(1-s12_TM_s_n(i)*r_p*exp(2*1i*kz1*(Z(i+1)-Z(i))));
    s12_TM_s_n(i+1) = (s12_TM_s_n(i)*exp(2*1i*kz1*(Z(i+1)-Z(i)))-r_p)/(1-s12_TM_s_n(i)*r_p*exp(2*1i*kz1*(Z(i+1)-Z(i))));
    s21_TM_s_n(i+1) = (s11_TM_s_n(i+1)*s22_TM_s_n(i)*r_p*exp(1i*kz1*(Z(i+1)-Z(i))))/t_p + s21_TM_s_n(i);
    s22_TM_s_n(i+1) = s22_TM_s_n(i)*(r_p*s12_TM_s_n(i+1)+1)*exp(1i*kz1*(Z(i+1)-Z(i)))/t_p;    
end

%% Finding the coefficients:

eps_s = Eps(s);
kzs = sqrt(eps_s*(w/c)^2-kp^2);

S_plus = exp(1i*kzs*Z(s));
S_minus = exp(1i*kzs*Z(s));

B_TE_s = s21_TE_s_n(end)*S_plus/(1-s21_TE_s_n(end)*s12_TE_0_l(s));
B_TE_0 = s22_TE_0_l(s)*B_TE_s;
A_TE_s = s12_TE_0_l(s)*B_TE_s;
C_TE_s = s12_TE_0_l(s)*S_minus/(1-s12_TE_0_l(s)*s21_TE_s_n(end));
D_TE_s = s21_TE_s_n(end)*C_TE_s;
D_TE_0 = s22_TE_0_l(s)*(D_TE_s+S_minus);

B_TM_s = s21_TM_s_n(end)*S_plus/(1-s21_TM_s_n(end)*s12_TM_0_l(s));
B_TM_0 = s22_TM_0_l(s)*B_TM_s;
A_TM_s = s12_TM_0_l(s)*B_TM_s;
C_TM_s = s12_TM_0_l(s)*S_minus/(1-s12_TM_0_l(s)*s21_TM_s_n(end));
D_TM_s = s21_TM_s_n(end)*C_TM_s;
D_TM_0 = s22_TM_0_l(s)*(D_TM_s+S_minus);

if s==1
    B_TE = (s21_TE_0_l(end)-s21_TE_0_l(l))/s22_TE_0_l(l);
    A_TE = s11_TE_0_l(l)+s12_TE_0_l(l)*B_TE;
    C_TE = 0;
    D_TE = 0;

    B_TM = (s21_TM_0_l(end)-s21_TM_0_l(l))/s22_TM_0_l(l);
    A_TM = s11_TM_0_l(l)+s12_TM_0_l(l)*B_TM;
    C_TM = 0;
    D_TM = 0;
elseif s~=1 && l<s
    B_TE = B_TE_0/s22_TE_0_l(l);
    A_TE = s12_TE_0_l(l)*B_TE;
    D_TE = D_TE_0/s22_TE_0_l(l);
    C_TE = s12_TE_0_l(l)*D_TE;
    
    B_TM = B_TM_0/s22_TM_0_l(l);
    A_TM = s12_TM_0_l(l)*B_TM;
    D_TM = D_TM_0/s22_TM_0_l(l);
    C_TM = s12_TM_0_l(l)*D_TM;    
else
    B_TE = (B_TE_s-s21_TE_s_n(l)*(A_TE_s+S_plus))/s22_TE_s_n(l);
    A_TE = s11_TE_s_n(l)*(A_TE_s+S_plus)+s12_TE_s_n(l)*B_TE;
    D_TE = (D_TE_s-s21_TE_s_n(l)*C_TE_s)/s22_TE_s_n(l);
    C_TE = s11_TE_s_n(l)*C_TE_s+s12_TE_s_n(l)*D_TE;

    B_TM = (B_TM_s-s21_TM_s_n(l)*(A_TM_s+S_plus))/s22_TM_s_n(l);
    A_TM = s11_TM_s_n(l)*(A_TM_s+S_plus)+s12_TM_s_n(l)*B_TM;
    D_TM = (D_TM_s-s21_TM_s_n(l)*C_TM_s)/s22_TM_s_n(l);
    C_TM = s11_TM_s_n(l)*C_TM_s+s12_TM_s_n(l)*D_TM;
end

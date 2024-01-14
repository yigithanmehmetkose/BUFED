function q_w = BUFED_v5(l,s,d,T,Eps)

%Constants
h_bar = 1.054571817e-34;        %reduced Planck's constant
k_b = 1.3807e-23;               %Boltzmann's constant
c = 2.99792e8;                  %speed of light

d = d - d(s);

N_kr = 301;
ts = d(s+1)-d(s);

W = Eps(:,1);
Eps_l = Eps(:,l+1);
Eps_s = Eps(:,s+1);
Eps_r = Eps(:,(2:end));
q_w = zeros(1,length(W));

parfor p=1:length(W)
    w = W(p);
    kv = w/c;
    theta = h_bar*w/(exp(h_bar*w/(k_b*T))-1);

    eps_l = Eps_l(p);
    eps_s = Eps_s(p);
    eps_r = Eps_r(p,:);

    kr_u = 0.999*kv;
    kr_l = 0;
    dkr = (kr_u-kr_l)/(3*N_kr);
    k_r = kr_l:dkr:kr_u;
    kr = zeros(1,4);
    I_kr = zeros(1,N_kr);

    for t=1:N_kr
        kr(1) = k_r(3*t-2);                     %discretized parallel wavevector for simpson integration
        kr(2) = k_r(3*t-1);
        kr(3) = k_r(3*t);
        kr(4) = k_r(3*t+1);
        F = zeros(1,4);

        for i=1:4
            ks = sqrt(eps_s)*kv;                %wavevector in emitting layer
            kl = sqrt(eps_l)*kv;                %wavevector in absorbing layer
            kzs = sqrt(ks^2-kr(i)^2);           %z component of wavevector in emitting layer
            kzl = sqrt(kl^2-kr(i)^2);           %z component of wavevector in absorbing layer

            if s==1
                [A_TE,~,~,~,A_TM,~,~,~] = S_matrix_method(w,kr(i),eps_r,d,s,l);

                g_E_0_rr = 1i*kzl/(2*kl*ks)*A_TM;
                g_E_0_rz = 1i*kzl*kr(i)/(2*kzs*ks*kl)*(-A_TM);
                g_E_0_tt = 1i/(2*kzs)*A_TE;

                g_H_0_tr = kl/(2*ks)*(-A_TM);
                g_H_0_tz = kl*kr(i)/(2*ks*kzs)*A_TM;
                g_H_0_rt = kzl/(2*kzs)*A_TE;

                F(i) =  kr(i)/(2*imag(kzs))*(g_E_0_rr*conj(g_H_0_tr)+g_E_0_rz*conj(g_H_0_tz)-g_E_0_tt*conj(g_H_0_rt));
            else
                [A_TE,B_TE,C_TE,D_TE,A_TM,B_TM,C_TM,D_TM] = S_matrix_method(w,kr(i),eps_r,d,s,l);

                g_E_rr_g_H_tr_star_plus_g_E_rz_g_H_tz_star = 1i*conj(kl)*kzl/(8*real(kzs)*imag(kzs)*kl*abs(ks)^2*abs(kzs)^2) * (real(kzs)*(exp(2*imag(kzs)*ts)-1)*(abs(kzs)^2+kr(i)^2)*...
                    (-abs(A_TM)^2 - A_TM*conj(B_TM) + conj(A_TM)*B_TM + abs(B_TM)^2) + 1i*imag(kzs)*(exp(-2i*real(kzs)*ts)-1)*(abs(kzs)^2-kr(i)^2)*(A_TM*conj(C_TM) + A_TM*conj(D_TM) - B_TM*conj(C_TM) - B_TM*conj(D_TM)) + ...
                    1i*imag(kzs)*(1-exp(2i*real(kzs)*ts))*(abs(kzs)^2-kr(i)^2)*(conj(A_TM)*C_TM + conj(B_TM)*C_TM - conj(A_TM)*D_TM - conj(B_TM)*D_TM) + real(kzs)*(1-exp(-2*imag(kzs)*ts))*(abs(kzs)^2+kr(i)^2)*(-abs(C_TM)^2 - ...
                    C_TM*conj(D_TM) + conj(C_TM)*D_TM + abs(D_TM)^2));

                g_E_tt_g_H_rt_star = 1i*conj(kzl)/(8*real(kzs)*imag(kzs)*abs(kzs)^2) * (real(kzs)*(exp(2*imag(kzs)*ts)-1)*(abs(A_TE)^2 - A_TE*conj(B_TE) + conj(A_TE)*B_TE - abs(B_TE)^2) + ...
                    1i*imag(kzs)*(exp(-2i*real(kzs)*ts)-1)*(A_TE*conj(C_TE) - A_TE*conj(D_TE) + B_TE*conj(C_TE) - B_TE*conj(D_TE)) + 1i*imag(kzs)*(1-exp(2i*real(kzs)*ts))*(conj(A_TE)*C_TE - conj(B_TE)*C_TE + conj(A_TE)*D_TE - ...
                    conj(B_TE)*D_TE) + real(kzs)*(1-exp(-2*imag(kzs)*ts))*(abs(C_TE)^2 - C_TE*conj(D_TE) + conj(C_TE)*D_TE - abs(D_TE)^2));

                F(i) = kr(i)*(g_E_rr_g_H_tr_star_plus_g_E_rz_g_H_tz_star - g_E_tt_g_H_rt_star);
            end
            if imag(eps_s) <= 0
                F(i) = 0;
            end
        end
        I_kr(t) = 3*dkr/8*(F(1)+3*F(2)+3*F(3)+F(4));
    end
    q_w(p) = (kv/pi)^2*theta*real(1i*imag(eps_s)*sum(I_kr));
end

end

PROGRAM BUFED
    IMPLICIT NONE

    REAL(8):: KR_U, KR_L, DKR, EPS_INF, W_LO, W_TO, GAM, H_BAR, K_B, C, DW, THT, PI, W, KV, KR_J, T, T_S, A
    REAL(8), ALLOCATABLE, DIMENSION(:):: D, Q_W, KR, K_R
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: EPS
    COMPLEX(8):: EPS_L, EPS_S, I, A_TE, B_TE, C_TE, D_TE, A_TM, B_TM, C_TM, D_TM, G_E_0_RR, G_E_0_RZ, G_E_0_TT, G_H_0_TR, &
        G_H_0_TZ, G_H_0_RT, KS, KL, KZS, KZL, G_E_RR_G_H_TR_STAR_PLUS_G_E_RZ_G_H_TZ_STAR, G_E_TT_G_H_RT_STAR
    COMPLEX(8), ALLOCATABLE, DIMENSION(:):: INT_KR, F
    INTEGER:: L, N_L, I_W, I_KR, S, N_KR, I1, I2
   
    !CONSTANTS
    H_BAR = 1.054571817D-34
    K_B = 1.3807D-23
    C = 2.99792D8
    PI = 4.0D0*ATAN(1.0D0)
    I = (0.0D0,1.0D0)

    !CBN PARAMETERS
    EPS_INF = 4.46D0
    W_LO = 2.451D14
    W_TO = 1.985D14
    GAM = 9.934D11

    !INPUTS
    N_L = 4
    S = 2
    D = 1.00D-9*(/0.0D0, 0.0D0, 10.0D0, 110.0D0/)
    L = 4
    T = 300.0D0
    DW = 1.0D12
    N_KR = 1000

    !ALLOCATION
    ALLOCATE(EPS(4,N_L))
    ALLOCATE(KR(N_KR))
    ALLOCATE(INT_KR(N_KR))
    ALLOCATE(K_R(4))
    ALLOCATE(Q_W(260))
    ALLOCATE(F(4))
   
    EPS(1,:) = (/1.0D0, EPS_INF, 1.0D0, EPS_INF/)
    EPS(2,:) = (/0.0D0, W_LO, 0.0D0, W_LO/)
    EPS(3,:) = (/0.0D0, W_TO, 0.0D0, W_TO/)
    EPS(4,:) = (/0.0D0, GAM, 0.0D0, GAM/)

    T_S = D(S+1)-D(S)
   
    WLOOP: DO I_W=1,260
        W = DW*I_W
        KV = W/C
        THT = H_BAR*W/(EXP(H_BAR*W/(K_B*T))-1.0D0)
        EPS_L = EPS_INF*(W**2.0D0-W_LO**2.0D0+I*GAM*W)/(W**2.0D0-W_TO**2.0D0+I*GAM*W)
        EPS_S = EPS_INF*(W**2.0D0-W_LO**2.0D0+I*GAM*W)/(W**2.0D0-W_TO**2.0D0+I*GAM*W)
       
        KR_U = 100.0D0*KV
        KR_L = 0.1D-3*KV
        DKR = (KR_U-KR_L)/(3.0*N_KR)
        KR = [(KR_L + DKR*(I1-1), I1=1,3*N_KR+1)]
           
        KLOOP: DO I_KR=1,N_KR
            K_R(1) = KR(3*I_KR-2)
            K_R(2) = KR(3*I_KR-1)
            K_R(3) = KR(3*I_KR)
            K_R(4) = KR(3*I_KR+1)
           
            KR_SIMPS_LOOP: DO I2=1,4
                KR_J = K_R(I2)
                KS = SQRT(EPS_S)*KV
                KL = SQRT(EPS_L)*KV
                KZS = SQRT(KS**2.0D0-KR_J**2.0D0)
                KZL = SQRT(KL**2.0D0-KR_J**2.0D0)
               
                IF (S .EQ. 1) THEN
               
                    CALL S_MATRIX_METHOD(A_TE,B_TE,C_TE,D_TE, &
                        A_TM,B_TM,C_TM,D_TM,W,KR_J, &
                        EPS,D,S,L)
                   
                    G_E_0_RR = I*KZL/(2.0D0*KL*KS)*A_TM
                    G_E_0_RZ = I*KZL*KR_J/(2.0D0*KZS*KS*KL)*(-A_TM)
                    G_E_0_TT = I/(2.0D0*KZS)*A_TE
               
                    G_H_0_TR = KL/(2.0D0*KS)*(-A_TM)
                    G_H_0_TZ = KL*KR_J/(2.0D0*KS*KZS)*A_TM
                    G_H_0_RT = KZL/(2.0D0*KZS)*A_TE
                       
                    F(I2) = KR_J/IMAG(KZS)*(G_E_0_RR*CONJG(G_H_0_TR)+G_E_0_RZ*CONJG(G_H_0_TZ)-G_E_0_TT*CONJG(G_H_0_RT))
                ELSE
                    CALL S_MATRIX_METHOD(A_TE,B_TE,C_TE,D_TE, &
                        A_TM,B_TM,C_TM,D_TM,W,KR_J, &
                        EPS,D,S,L)
                       
                    G_E_RR_G_H_TR_STAR_PLUS_G_E_RZ_G_H_TZ_STAR = I* &
                        CONJG(KL)*KZL/(8.0D0*REAL(KZS)*IMAG(KZS)*KL*ABS(KS)**2.0D0*ABS(KZS)**2.0D0) * &
                        (REAL(KZS)*(EXP(2.0D0*IMAG(KZS)*T_S)-1.0D0)*(ABS(KZS)**2.0D0+KR_J**2.0D0)* &
                        (-ABS(A_TM)**2.0D0 - A_TM*CONJG(B_TM) + CONJG(A_TM)*B_TM + ABS(B_TM)**2.0D0) + &
                        I*IMAG(KZS)*(EXP(-2.0D0*I*REAL(KZS)*T_S)-1.0D0)*(ABS(KZS)**2.0D0-KR_J**2.0D0)*(A_TM*CONJG(C_TM) + A_TM* &
                        CONJG(D_TM) - B_TM*CONJG(C_TM) - B_TM*CONJG(D_TM)) + I* &
                        IMAG(KZS)*(1.0D0-EXP(2.0D0*I*REAL(KZS)*T_S))*(ABS(KZS)**2.0D0-KR_J**2.0D0)*(CONJG(A_TM)*C_TM + CONJG(B_TM)*&
                        C_TM - CONJG(A_TM)*D_TM - CONJG(B_TM)*D_TM) + REAL(KZS)* &
                        (1.0D0-EXP(-2.0D0*IMAG(KZS)*T_S))*(ABS(KZS)**2.0D0+KR_J**2.0D0)*(-ABS(C_TM)**2.0D0 - C_TM*CONJG(D_TM) + &
                        CONJG(C_TM)*D_TM + ABS(D_TM)**2.0D0))
                       
                    G_E_TT_G_H_RT_STAR = I*CONJG(KZL)/(8.0D0*REAL(KZS)*IMAG(KZS)*ABS(KZS)**2.0D0) * &
                        (REAL(KZS)*(EXP(2.0D0*IMAG(KZS)*T_S)-1.0D0)*(ABS(A_TE)**2.0D0 - A_TE*CONJG(B_TE) + CONJG(A_TE)*B_TE - &
                        ABS(B_TE)**2.0D0) + I*IMAG(KZS)*(EXP(-2.0D0*I*REAL(KZS)*T_S)-1.0D0)*(A_TE*CONJG(C_TE) - A_TE*CONJG(D_TE) + &
                        B_TE*CONJG(C_TE) - B_TE*CONJG(D_TE)) + I*IMAG(KZS)*(1.0D0-EXP(2.0D0*I*REAL(KZS)*T_S))*(CONJG(A_TE)*C_TE - &
                        CONJG(B_TE)*C_TE + CONJG(A_TE)*D_TE - CONJG(B_TE)*D_TE) + REAL(KZS)*(1.0D0-EXP(-2.0D0*IMAG(KZS)*T_S))* &
                        (ABS(C_TE)**2.0D0 - C_TE*CONJG(D_TE) + CONJG(C_TE)*D_TE - ABS(D_TE)**2.0D0))
                   
                    IF (IMAG(EPS_S) .EQ. 0) THEN
                        G_E_RR_G_H_TR_STAR_PLUS_G_E_RZ_G_H_TZ_STAR = 0.0D0
                        G_E_TT_G_H_RT_STAR = 0.0D0

                    END IF
                    F(I2) = KR_J*(G_E_RR_G_H_TR_STAR_PLUS_G_E_RZ_G_H_TZ_STAR-G_E_TT_G_H_RT_STAR)
                 END IF            
                 INT_KR(I_KR) = 3.0*DKR/8.0*(F(1)+3.0*F(2)+3.0*F(3)+F(4))
            END DO KR_SIMPS_LOOP        
        END DO KLOOP    
        Q_W(I_W) = (KV/PI)**2.0D0*THT*REAL(I*IMAG(EPS_S)*SUM(INT_KR))  
        PRINT *, Q_W(I_W)
    END DO WLOOP

CONTAINS

    SUBROUTINE S_MATRIX_METHOD(A_TE,B_TE,C_TE,D_TE,A_TM,B_TM,C_TM,D_TM,W,KR,EPS, &
                D,S,L)

        REAL(8):: W, KR, C
        REAL(8), ALLOCATABLE, DIMENSION(:):: D, EPS_INF, W_LO, W_TO, GAM
        REAL(8), ALLOCATABLE, DIMENSION(:,:):: EPS
        COMPLEX(8):: I, EPS1, EPS2, KZ1, KZ2, RS, TS, RP, TP, EPS_S, KZS, S_PLUS, &
                S_MINUS, A_TE, B_TE, C_TE, D_TE, A_TM, B_TM, C_TM, &
                D_TM, B_TE_S, B_TE_0, A_TE_S, C_TE_S, D_TE_S, D_TE_0, &
                B_TM_S, B_TM_0, A_TM_S, C_TM_S, D_TM_S, D_TM_0
        COMPLEX(8), ALLOCATABLE, DIMENSION(:):: S11_TE_0_L, S12_TE_0_L, S21_TE_0_L, &
                S22_TE_0_L, S11_TE_S_N, S12_TE_S_N, S21_TE_S_N, &
                S22_TE_S_N, S11_TM_0_L, S12_TM_0_L, S21_TM_0_L, &
                S22_TM_0_L, S11_TM_S_N, S12_TM_S_N, S21_TM_S_N, S22_TM_S_N  
        INTEGER:: S,L,N,J

        !CONSTANTS
        C = 2.99792D8
        I = (0.0D0,1.0D0)

        N = SIZE(D)
        EPS_INF = EPS(1,:)
        W_LO = EPS(2,:)
        W_TO = EPS(3,:)
        GAM = EPS(4,:)
               
        ALLOCATE(S11_TE_0_L(N))
        ALLOCATE(S12_TE_0_L(N))
        ALLOCATE(S21_TE_0_L(N))
        ALLOCATE(S22_TE_0_L(N))
        ALLOCATE(S11_TE_S_N(N))
        ALLOCATE(S12_TE_S_N(N))
        ALLOCATE(S21_TE_S_N(N))
        ALLOCATE(S22_TE_S_N(N))
        ALLOCATE(S11_TM_0_L(N))
        ALLOCATE(S12_TM_0_L(N))
        ALLOCATE(S21_TM_0_L(N))
        ALLOCATE(S22_TM_0_L(N))
        ALLOCATE(S11_TM_S_N(N))
        ALLOCATE(S12_TM_S_N(N))
        ALLOCATE(S21_TM_S_N(N))
        ALLOCATE(S22_TM_S_N(N))

        S11_TE_0_L(1) = 1.0D0
        S12_TE_0_L(1) = 0.0D0
        S21_TE_0_L(1) = 0.0D0
        S22_TE_0_L(1) = 1.0D0

        S11_TM_0_L(1) = 1.0D0
        S12_TM_0_L(1) = 0.0D0
        S21_TM_0_L(1) = 0.0D0
        S22_TM_0_L(1) = 1.0D0

        S11_TE_S_N(S) = 1.0D0
        S12_TE_S_N(S) = 0.0D0
        S21_TE_S_N(S) = 0.0D0
        S22_TE_S_N(S) = 1.0D0

        S11_TM_S_N(S) = 1.0D0
        S12_TM_S_N(S) = 0.0D0
        S21_TM_S_N(S) = 0.0D0
        S22_TM_S_N(S) = 1.0D0

        DO J=1,N-1
            EPS1 = EPS_INF(J)*(W**2.0D0-W_LO(J)**2.0D0+I*GAM(J)*W)/(W**2.0D0-W_TO(J)**2.0D0+I*GAM(J)*W)
            EPS2 = EPS_INF(J+1)*(W**2.0D0-W_LO(J+1)**2.0D0+I*GAM(J+1)*W)/(W**2.0D0-W_TO(J+1)**2.0D0+I*GAM(J+1)*W)
            KZ1 = SQRT(EPS1*(W/C)**2.0D0-KR**2.0D0)
            KZ2 = SQRT(EPS2*(W/C)**2.0D0-KR**2.0D0)
           
            RS = (KZ1-KZ2)/(KZ1+KZ2)
            TS = 1.0D0+(KZ1-KZ2)/(KZ1+KZ2)
            RP = (EPS2*KZ1-EPS1*KZ2)/(EPS2*KZ1+EPS1*KZ2)
            TP = 2.0D0*SQRT(EPS1*EPS2)*KZ1/(EPS2*KZ1+EPS1*KZ2)
           
            S11_TE_0_L(J+1) = S11_TE_0_L(J)*TS*EXP(I*KZ1*(D(J+1)-D(J)))/(1.0D0-S12_TE_0_L(J)*RS*EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S12_TE_0_L(J+1) = (S12_TE_0_L(J)*EXP(2.0D0*I*KZ1*(D(J+1)-D(J)))-RS)/(1.0D0-S12_TE_0_L(J)*RS* &
                EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S21_TE_0_L(J+1) = (S11_TE_0_L(J+1)*S22_TE_0_L(J)*RS*EXP(I*KZ1*(D(J+1)-D(J))))/TS + S21_TE_0_L(J)
            S22_TE_0_L(J+1) = S22_TE_0_L(J)*(RS*S12_TE_0_L(J+1)+1.0D0)*EXP(I*KZ1*(D(J+1)-D(J)))/TS

            S11_TM_0_L(J+1) = S11_TM_0_L(J)*TP*EXP(I*KZ1*(D(J+1)-D(J)))/(1.0D0-S12_TM_0_L(J)*RP*EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S12_TM_0_L(J+1) = (S12_TM_0_L(J)*EXP(2.0D0*I*KZ1*(D(J+1)-D(J)))-RP)/(1.0D0-S12_TM_0_L(J)*RP*&
                EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S21_TM_0_L(J+1) = (S11_TM_0_L(J+1)*S22_TM_0_L(J)*RP*EXP(I*KZ1*(D(J+1)-D(J))))/TP + S21_TM_0_L(J)
            S22_TM_0_L(J+1) = S22_TM_0_L(J)*(RP*S12_TM_0_L(J+1)+1.0D0)*EXP(I*KZ1*(D(J+1)-D(J)))/TP
        END DO

        DO J=S,N-1
            EPS1 = EPS_INF(J)*(W**2.0D0-W_LO(J)**2.0D0+I*GAM(J)*W)/(W**2.0D0-W_TO(J)**2.0D0+I*GAM(J)*W)
            EPS2 = EPS_INF(J+1)*(W**2.0D0-W_LO(J+1)**2.0D0+I*GAM(J+1)*W)/(W**2.0D0-W_TO(J+1)**2.0D0+I*GAM(J+1)*W)
            KZ1 = SQRT(EPS1*(W/C)**2.0D0-KR**2.0D0)
            KZ2 = SQRT(EPS2*(W/C)**2.0D0-KR**2.0D0)
           
            RS = (KZ1-KZ2)/(KZ1+KZ2)
            TS = 1.0D0+(KZ1-KZ2)/(KZ1+KZ2)
            RP = (EPS2*KZ1-EPS1*KZ2)/(EPS2*KZ1+EPS1*KZ2)
            TP = 2.0D0*SQRT(EPS1*EPS2)*KZ1/(EPS2*KZ1+EPS1*KZ2)
           
            S11_TE_S_N(J+1) = S11_TE_S_N(J)*TS*EXP(I*KZ1*(D(J+1)-D(J)))/(1.0D0-S12_TE_S_N(J)*RS*EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S12_TE_S_N(J+1) = (S12_TE_S_N(J)*EXP(2.0D0*I*KZ1*(D(J+1)-D(J)))-RS)/(1.0D0-S12_TE_S_N(J)*RS*&
                EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S21_TE_S_N(J+1) = (S11_TE_S_N(J+1)*S22_TE_S_N(J)*RS*EXP(I*KZ1*(D(J+1)-D(J))))/TS + S21_TE_S_N(J)
            S22_TE_S_N(J+1) = S22_TE_S_N(J)*(RS*S12_TE_S_N(J+1)+1.0D0)*EXP(I*KZ1*(D(J+1)-D(J)))/TS

            S11_TM_S_N(J+1) = S11_TM_S_N(J)*TP*EXP(I*KZ1*(D(J+1)-D(J)))/(1.0D0-S12_TM_S_N(J)*RP*EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S12_TM_S_N(J+1) = (S12_TM_S_N(J)*EXP(2.0D0*I*KZ1*(D(J+1)-D(J)))-RP)/(1.0D0-S12_TM_S_N(J)*RP*&
                EXP(2.0D0*I*KZ1*(D(J+1)-D(J))))
            S21_TM_S_N(J+1) = (S11_TM_S_N(J+1)*S22_TM_S_N(J)*RP*EXP(I*KZ1*(D(J+1)-D(J))))/TP + S21_TM_S_N(J)
            S22_TM_S_N(J+1) = S22_TM_S_N(J)*(RP*S12_TM_S_N(J+1)+1.0D0)*EXP(I*KZ1*(D(J+1)-D(J)))/TP
        END DO

        EPS_S = EPS_INF(S)*(W**2.0D0-W_LO(S)**2.0D0+I*GAM(S)*W)/(W**2.0D0-W_TO(S)**2.0D0+I*GAM(S)*W)
        KZS = SQRT(EPS_S*(W/C)**2.0D0-KR**2.0D0)

        S_PLUS = EXP(I*KZS*D(S))
        S_MINUS = EXP(I*KZS*(D(S)))

        B_TE_S = S21_TE_S_N(N)*S_PLUS/(1.0D0-S21_TE_S_N(N)*S12_TE_0_L(S))
        B_TE_0 = S22_TE_0_L(S)*B_TE_S
        A_TE_S = S12_TE_0_L(S)*B_TE_S
        C_TE_S = S12_TE_0_L(S)*S_MINUS/(1.0D0-S12_TE_0_L(S)*S21_TE_S_N(N))
        D_TE_S = S21_TE_S_N(N)*C_TE_S
        D_TE_0 = S22_TE_0_L(S)*(D_TE_S+S_MINUS)

        B_TM_S = S21_TM_S_N(N)*S_PLUS/(1.0D0-S21_TM_S_N(N)*S12_TM_0_L(S))
        B_TM_0 = S22_TM_0_L(S)*B_TM_S
        A_TM_S = S12_TM_0_L(S)*B_TM_S
        C_TM_S = S12_TM_0_L(S)*S_MINUS/(1.0D0-S12_TM_0_L(S)*S21_TM_S_N(N))
        D_TM_S = S21_TM_S_N(N)*C_TM_S
        D_TM_0 = S22_TM_0_L(S)*(D_TM_S+S_MINUS)

        IF (S .EQ. 1) THEN
            B_TE = (S21_TE_0_L(N)-S21_TE_0_L(L))/S22_TE_0_L(L)
            A_TE = S11_TE_0_L(L)+S12_TE_0_L(L)*B_TE
            C_TE = 0.0D0
            D_TE = 0.0D0
           
            B_TM = (S21_TM_0_L(N)-S21_TM_0_L(L))/S22_TM_0_L(L)
            A_TM = S11_TM_0_L(L)+S12_TM_0_L(L)*B_TM
            C_TM = 0.0D0
            D_TM = 0.0D0            
        ELSEIF (S .NE. 1 .AND. L .LT. S) THEN
            B_TE = B_TE_0/S22_TE_0_L(L)
            A_TE = S12_TE_0_L(L)*B_TE
            D_TE = D_TE_0/S22_TE_0_L(L)
            C_TE = S12_TE_0_L(L)*D_TE
           
            B_TM = B_TM_0/S22_TM_0_L(L)
            A_TM = S12_TM_0_L(L)*B_TM
            D_TM = D_TM_0/S22_TM_0_L(L)
            C_TM = S12_TM_0_L(L)*D_TM
        ELSE
            B_TE = (B_TE_S-S21_TE_S_N(L)*(A_TE_S+S_PLUS))/S22_TE_S_N(L)
            A_TE = S11_TE_S_N(L)*(A_TE_S+S_PLUS)+S12_TE_S_N(L)*B_TE
            D_TE = (D_TE_S-S21_TE_S_N(L)*C_TE_S)/S22_TE_S_N(L)
            C_TE = S11_TE_S_N(L)*C_TE_S+S12_TE_S_N(L)*D_TE

            B_TM = (B_TM_S-S21_TM_S_N(L)*(A_TM_S+S_PLUS))/S22_TM_S_N(L)
            A_TM = S11_TM_S_N(L)*(A_TM_S+S_PLUS)+S12_TM_S_N(L)*B_TM
            D_TM = (D_TM_S-S21_TM_S_N(L)*C_TM_S)/S22_TM_S_N(L)
            C_TM = S11_TM_S_N(L)*C_TM_S+S12_TM_S_N(L)*D_TM
        END IF
    END SUBROUTINE S_MATRIX_METHOD
END PROGRAM BUFED
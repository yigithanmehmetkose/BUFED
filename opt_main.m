clc
clear all
close all

load('InAs_revised.mat')
load('E_ideal_eff.mat')
load('Ge_revised.mat')

if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

%psoptions = optimoptions('particleswarm','UseParallel', true, 'UseVectorized', false);

lb = 1e-9*[5 5 5];
ub = 1e-9*[250 250 250];

%x0 = 1e-9*[245 55.1 77.7 35.3 5.5];

x = particleswarm(@obj_v2,nvars,lb,ub)
%options = optimoptions('patternsearch','MeshTolerance',1e-9,'StepTolerance',1e-8);

%x = patternsearch(@obj_v2,x0,[],[],[],[],lb,ub,[],options)

[P_e,eff,emissivity] = BUFED_TPV_main(x(1),x(2),x(3));

save('AgGeSiGe.mat')

%plot(W,emissivity)
%hold on
%plot(W,InAs_IQE)
%hold on
%plot(W,E_ideal)

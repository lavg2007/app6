%% APP6
clearvars
clc
close all


identification
conditions_initiales

gamma_ref = deg2rad(-40); %% angle constant pour la RAA

h = [h_ini:-1:0];

p0 = 0.02
p = p0*exp(-h/hs);

%% RAA
v = v_ini*exp(0.5*B*hs*(p-p(1))/sin(gamma_ref));

figure
plot(h,v)
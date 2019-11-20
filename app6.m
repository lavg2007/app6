%% APP6
clearvars
clc
close all


identification

v_fin = v_fin1
% v_fin = v_fin2


h = [h_ini:-1:v_fin];
%p = p0*exp(-h/hs);

%% RAA

gamma_ref = deg2rad(-90); %% angle constant pour la RAA
v = v_ini*exp(0.5*B*hs*(p-p(1))/sin(gamma_ref));



figure
subplot(2,1,1)
title('Vitesse vs Hauteur')
hold on
plot(pos_mes,-v)
plot(pos_mes,vel_mes)
legend('Théorique','Expérimentale')
subplot(2,1,2)
title('Erreur')
plot(pos_mes,v+vel_mes)



r = h+R_mars;
delta_v_aero_ini = v_fin - sqrt(v_ini.^2 + 2*mu_mars*((1/r_fin)-(1./r_ini)));
% 
% figure
% plot(h, delta_v_aero)
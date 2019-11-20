%% APP6
clearvars
clc
close all


identification

v_fin = v_fin1;
% v_fin = v_fin2

rho = @(h) p0*exp(-h/hs);
h = [h_ini:-1:v_fin];


%% RAA

% gamma_ref = deg2rad(-90); %% angle constant pour la RAA

RAA = @(v_i, gamma, p) v_i*exp(0.5*B*hs*(p-p(1))/sin(-90));

% RAA a angle de -90 pour comparaison avec les russes
v = RAA(v_ini, -90, p);

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


%% loi de guidage (RAA + Grav)
p_fin = rho(h_fin);

r = h+R_mars;

deltaVaero = @(v_fin, v, h) v_fin - sqrt(v.^2 + 2*mu_mars*((1/r_fin)-(1./(h+R_mars))));
gammaRef = @(dva,h,v) asin(0.5*B*hs*(p_fin-rho(h))/log(1+dva/v));

d_v_a1 = deltaVaero(v_fin1, v_ini, h_fin);
d_v_a2 = deltaVaero(v_fin2, v_ini, h_fin);

gamma_ref1 = gammaRef(d_v_a1, h_ini, v_ini);
gamma_ref2 = gammaRef(d_v_a2, h_ini, v_ini);

disp('--Angles de commande--')
disp(' |    300 m/s   |    250 m/s   | ')  
disp([' | ' num2str(rad2deg(gamma_ref1)) ' deg | ' num2str(rad2deg(gamma_ref2)) ' deg |'])
disp(newline)
% 
% figure
% plot(h, delta_v_aero)


%détermination des deux angles
% delta_v_aero


%% Newton Raphson

tol = 1e-8;
% func1 = @(




function [x i] = NP(a,b,c,tol)
    i = 0;
    x = a;
    while 1
        if abs(func1) <= tol
           break
        else
            i = i+1;
            x = x - func1/dfunc1
        end
    end
end
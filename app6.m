%% APP6
clearvars
clc
close all


identification

v_fin = v_fin1;
% v_fin = v_fin2

rho = @(h) p0*exp(-h/hs);
rho_inv = @(p) hs*(log(p0) - log(p));
h = [h_ini:-1:v_fin];
p_ini = rho(h_ini);

%% RAA

% gamma_ref = deg2rad(-90); %% angle constant pour la RAA

RAA = @(v_i, gamma, p) v_i*exp(0.5*B*hs*(p-p_ini)/sin(gamma));

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


% 
% figure
% plot(h, delta_v_aero)



%% Newton Raphson


tol = 1e-8;

% fonctions angle nominal
func1 = @(pe) pe.*(v_ini.*exp(0.5.*B.*hs.*(pe-p_ini)./sin(gamma_ref1))).^2 - 4000/(S*C_D0);
dfunc1 = @(pe) (v_ini^2*exp(-B*hs*(pe-p_ini)/sin(gamma_ref1))*(sin(gamma_ref1) + B*hs*pe))/(sin(gamma_ref1));

% fonctions angle ideal
func2 = @(pe) pe.*(v_ini.*exp(0.5.*B.*hs.*(pe-p_ini)./sin(gamma_ref2))).^2 - 4000/(S*C_D0);
dfunc2 = @(pe) (v_ini^2*exp(-B*hs*(pe-p_ini)/sin(gamma_ref2))*(sin(gamma_ref2) + B*hs*pe))/(sin(gamma_ref2));

h_test = [0:0.1:120000];
p_test = rho(h_test);
func_test = func1(p_test);

figure
plot(h_test, func_test)


%% t_lim pour angle nominal

[p_max1 i1] = newtonRaphson(0,func1, dfunc1, tol);
h_max1 = rho_inv(p_max1);
v_max1 = RAA(v_ini,gamma_ref1,p_max1);

[p_max2 i2] = newtonRaphson(0.005,func1, dfunc1, tol);
h_max2 = rho_inv(p_max2);
v_max2 = RAA(v_ini,gamma_ref1,p_max2);

h_test = [120000:-1:0];


figure
hold on
plot(h_test, RAA(v_ini, gamma_ref1, rho(h_test)))
plot(h_test, RAA(v_ini, gamma_ref2, rho(h_test)))

v_mean = mean([v_max1 v_max2]);
t_lim_nom = (h_max1 - h_max2)/v_mean;

% figure(5)
% hold on
% plot(h_test, RAA(v_ini,gamma_ref1,rho(h_test)))
% plot(h_max1, RAA(v_ini,gamma_ref1,rho(h_max1)), 's')
% plot(h_max2, RAA(v_ini,gamma_ref1,rho(h_max2)), 's')

%% t_lim pour angle idéal

[p_max1 i1] = newtonRaphson(0,func2, dfunc2, tol);
h_max1 = rho_inv(p_max1);
v_max1 = RAA(v_ini,gamma_ref2,p_max1);

[p_max2 i2] = newtonRaphson(0.005,func2, dfunc2, tol);
h_max2 = rho_inv(p_max2);
v_max2 = RAA(v_ini,gamma_ref2,p_max2);

h_test = [120000:-1:0];

v_mean = mean([v_max1 v_max2]);
t_lim_ide = (h_max1 - h_max2)/v_mean;

% figure(5)
% hold on
% plot(h_test, RAA(v_ini,gamma_ref2,rho(h_test)))
% plot(h_max1, RAA(v_ini,gamma_ref2,rho(h_max1)), 's')
% plot(h_max2, RAA(v_ini,gamma_ref2,rho(h_max2)), 's')

%% Commande en tranlation

% retroaction linearisante, forme standard Kp/(s+Kp) => Kp = 1/tau

tau = 0.2;
Kp = 1/tau;


%% commande en rotation

wn = 20;
zeta = 0.7;

Kp = wn^2;
Kd = 2*zeta*wn;
%% Console


disp('Vitesses finales |      300 m/s |      250 m/s | ')  
disp('--------------------------------------------------')
disp(['Angle constant   | ' num2str(rad2deg(gamma_ref1)) ' deg | ' num2str(rad2deg(gamma_ref2)) ' deg |'])
disp(['Tlim D > 2000N   |  ' num2str(t_lim_nom) ' s   |  ' num2str(t_lim_ide) ' s   |'])
disp(newline)

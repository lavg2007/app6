m = 50; %masse de la capsule
j = 1.5; %inertie de la capsule
R_mars = 3397e03; %rayon de mars
mu_mars = 42830e09; %param grav mars
S = 0.8; % surface aerodynamique capsule
d = 0.05; % dimension aero capsule
C_D0 = 1.2; % coeff trainée
C_La = 0.8; % coeff portance
C_Ma = -0.07; % coeff couple
C_Mq = -0.05; % coeff amortissement
C_Md = 0.1 ;% coeff volet aero

B = S*C_D0/m;

%% conditions initiales
v_ini = 6100;
gamma_ini = deg2rad(-20.5);
h_ini  = 120000;
r_ini = h_ini + R_mars
s_ini = 0;
theta_ini = deg2rad(-80);
q_ini = 0;

v_fin1 = 300;
v_fin2 = 250;
h_fin = 10000;
r_fin = h_fin + R_mars;

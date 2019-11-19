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
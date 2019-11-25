function f = eqnynominale(t,x)
%% importation des constantes
cte;
hs = 1.100353042442160e+04;
p0 = 0.021571989401399;
pfin = p0 .* exp(-10000./hs);
pinit = p0 * exp(-120000/hs);
%% calcul de yref
dvaero = vfin.nominale - sqrt(x(1)^2 +2*um*(1/rfin - 1/(r+x(3))));
yref = asin(B*0.5*hs*(pfin - p0 * exp(-x(3)/hs))/log(1+(dvaero/x(1))));
%% calcul des variables utilisees
a = x(5) - x(2);                                                            % angle d'attaque
Pdyn =  0.5 * p0 * exp(-x(3)/hs)* x(1)^2;                                   % Pression aero dynamique
rf = r+x(3);                                                                % Rayon mars + altitude
g = um/rf^2;                                                                % force gravtitationnelle
% valeurs desire
gamma_des = yref;
q_des = 0; % je suis vrm pas certain de celle la....
% caracteristiques des compensateurs
kpy = 1/0.2;
kpd =20^2;
kdd = 2*0.7*20;
%% commande de controle theta
fy = -(Pdyn*S*Cla*x(2))/(x(1)*m) + (x(1)^2/rf-g)*(cos(x(2))/x(1));
gy = (Pdyn*S*Cla)/(x(1)*m);
tta_cmd = -fy/gy + kpy/gy*(gamma_des-x(2));
if abs(tta_cmd) > deg2rad(60)
    tta_cmd = deg2rad(60) * sign(tta_cmd);
end
%% commande de controle delta
const = 1/J * Pdyn * S * d;
fq = const * Cma * (x(5) - x(2)) + const * d/(2*x(1))*Cmq*x(6);
gq = const * Cmd;
dta_cmd = -fq/gq + kpd/gq*(tta_cmd - x(5)) + kdd/gq*(q_des - x(6));
% dta_cmd = 0;
%% calcul de l'aerodynamique
D = Pdyn * S * Cdo;                                                         % Daero
L = Pdyn * S * Cla * (x(5)-x(2));                                           % Laero
M = Pdyn * S * d * (Cma*(x(5)-x(2)) + d/(2*x(1))*Cmq*x(6)+Cmd*dta_cmd);     % Maero  
%% calcul de la dynamique
f(1) = -D/m - g*sin(x(2));
f(2) = 1/x(1) * (L/m +(x(1)^2/rf-g)*cos(x(2)));
f(3) = x(1)*sin(x(2));
f(4) = x(1)/rf * cos(x(2));
f(5) = x(6);
f(6) = 1/J *M;
f(7) = D>2000;
  f = f(:);
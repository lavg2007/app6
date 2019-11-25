function f = eqny_no(t,x)
cte;
hs = 1.100353042442160e+04;
p0 = 0.021571989401399;
pinit = p0 * exp(120000/hs);

deltacmd = 0 ;                                                              % commande d'angle
a = x(5) - x(2);                                                            % angle d'attaque
Pdyn =  0.5 * p0 * exp(-x(3)/hs)* x(1)^2;  
D = Pdyn * S * Cdo;                                                         % Daero
L = Pdyn * S * Cla *a;                                                      % Laero
M = Pdyn * S * d * (Cma*(x(5) - x(2)) + d/(2*x(1))*Cmq*x(6)+Cmd*deltacmd);  % Maero
rf = r+x(3);                                                                % Rayon mars + altitude
g = um/rf^2;                                                                % force gravtitationnelle



% Eqn de dynamique
f(1) = -D/m - g*sind(x(2));
f(2) = 1/x(1) * (L/m +(x(1)^2/rf-g)*cosd(x(2)));
f(3) = x(1)*sind(x(2));
f(4) = x(1)/rf * cosd(x(2));
f(5) = x(6);
f(6) =1/J *M;
f(7) = D>2000;

  f = f(:);
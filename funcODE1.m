function f = funcODE1(t,y)
% y = [v y h s theta q]
% f = dy

params
r = R_mars + y(3);
g = mu_mars/r^2;
p = p0*exp(-y(3)/hs);
P_dyn = 0.5*p*y(1)^2;
a = y(5) - y(2);
D_aero = P_dyn*S*C_D0;
L_aero = P_dyn*S*C_La*a;
M_aero = P_dyn*S*d*(C_Ma*a + 0.5*d*C_Mq*y(6)/y(1));


f(1) = -D_aero/m - g*sin(y(2));
f(2) = (L_aero/m + (y(1)^2/r - g)*cos(y(2)))/y(1);
f(3) = y(1)*sin(y(2));
f(4) = y(1)*cos(y(2))/r;
f(5) = y(6);
f(6) = M_aero/J;
f = f';
end


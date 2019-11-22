function f = funcODE1(t,y)
% y = [v y h s theta q]
% f = dy


params
r = R_mars + y(3);
g = mu_mars/r^2;
p = p0*exp(-y(3)/hs);
P_dyn = 0.5*p*y(1)^2;
a = y(5) - y(2);



deltaVaero = v_fin2 - sqrt(y(1)^2 + 2*mu_mars*(1/r_fin - 1/r));
gamma_ref = asin(0.5*B*hs*(p_fin-p)/log(1+deltaVaero/(y(1))));

Kp_theta = 5;
f_theta = ((y(1)^2/r - g)*cos(y(2)) - P_dyn*S*C_La*y(2)/m)/y(1);
g_theta = P_dyn*S*C_La/(m*y(1));
theta_eq = -f_theta/g_theta;
theta_cmd = theta_eq + Kp_theta*(gamma_ref-y(2));
if abs(theta_cmd) > 1.0472
    theta_cmd = 1.0472*sign(theta_cmd);
end

Kp_delta = 400;
Kd_delta = 28;
f_delta = (P_dyn*S*d*(C_Ma*a + d*C_Mq*y(6)/(2*y(1))))/J;
g_delta = P_dyn*S*d*C_Md;
delta_eq = -f_delta/g_delta;
delta_cmd = delta_eq + Kp_delta*(theta_cmd-y(5)) + Kd_delta*(0-y(6));

D_aero = P_dyn*S*C_D0;
a_cmd = theta_cmd-y(2);
L_aero = P_dyn*S*C_La*a_cmd;
M_aero = P_dyn*S*d*(C_Ma*a + 0.5*d*C_Mq*y(6)/y(1) + C_Md*delta_cmd);

f(1) = -D_aero/m - g*sin(y(2));
f(2) = (L_aero/m + (y(1)^2/r - g)*cos(y(2)))/y(1);
f(3) = y(1)*sin(y(2));
f(4) = y(1)*cos(y(2))/r;
f(5) = y(6);
f(6) = M_aero/J;
f = f';
end


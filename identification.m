%% APP6
clc
clear


load('Accelero_Data_from_Moscow.mat')
params


%% intégration numérique des mesures

vel_mes(1) = -6100;
for i = 2:numel(acc_mes)
    vel_mes(i) = trapz(t(1:i),acc_mes(1:i)) + vel_mes(1);
end
vel_mes = vel_mes';
err_vel = integError(t,acc_mes)

pos_mes(1) = 120000;
for i = 2:numel(acc_mes)
    pos_mes(i) = trapz(t(1:i),vel_mes(1:i)) + pos_mes(1);
end
pos_mes = pos_mes';
err_vel = integError(t,vel_mes)


%%
D_aero = acc_mes*m;
P_dyn = D_aero/(S*C_D0);
p = 2*P_dyn./(vel_mes.^2);

[b m] = moindreCarre(pos_mes, log(p));
p0 = exp(b);
hs = -1/m;

p_lisse = p0*exp(-pos_mes/hs);

figure
hold on
plot(pos_mes, p)
plot(pos_mes, p_lisse)


%%

disp('--Identification des paramètres--')
disp([' H_s = ' num2str(hs) newline ' p_0 = ' num2str(p0)])

% figure
% subplot(3,1,1)
% plot(t,acc_mes)
% subplot(3,1,2)
% plot(t,vel_mes)
% subplot(3,1,3)
% plot(t,pos_mes)


clearvars variables acc_mes vel_mes pos_mes t D_aero P_dyn b m p p_lisse i 

function [b, m] = moindreCarre(xn, yn)
    A = [numel(xn) sum(xn); sum(xn) sum(xn.^2)];
    B = [sum(yn); sum(xn.*yn)];
    rep = inv(A)*B;
    b = rep(1);
    m = rep(2);
end

function F = trapezoidale(xn,yn)
    h = mean(diff(xn));
    if numel(xn) == 1
       F = 0;
    elseif numel(xn) == 2;
       F = (yn(1)+yn(2))*h/2;
    else
       F = 0.5*h*(yn(1) + yn(end) + 2*sum(yn(2:end-1)));
    end
end

function err = integError(xn,yn)
    h = mean(diff(xn));
    dy = diff(yn);
    err = h^2*(dy(end)-dy(1))/12;
end
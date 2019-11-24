%% APP6
clc
clear


load('Accelero_Data_from_Moscow.mat')
params


%% intégration numérique des mesures

vel_mes(1) = -6100;
for i = 2:numel(acc_mes)
    vel_mes(i) = trapezoidale(t(1:i),acc_mes(1:i)) + vel_mes(1);
end
vel_mes = vel_mes';
err_vel = integError(t,acc_mes);

pos_mes(1) = 120000;
for i = 2:numel(acc_mes)
    pos_mes(i) = trapezoidale(t(1:i),vel_mes(1:i)) + pos_mes(1);
end
pos_mes = pos_mes';
err_pos = integError(t,vel_mes);

disp('--Erreurs d''intégration--')
disp(['  vitesse  = ' num2str(err_vel) ' m/s'])
disp(['  position = ' num2str(err_pos) ' m'])
disp(newline)




%%
D_aero = acc_mes*m;
P_dyn = D_aero/(S*C_D0);
p_mes = 2*P_dyn./(vel_mes.^2);

[b M] = moindreCarre(pos_mes, log(p_mes));
p0 = exp(b);
hs = -1/M;

p_lisse = p0*exp(-pos_mes/hs);

figure
hold on
plot(pos_mes, p_mes)
plot(pos_mes, p_lisse)

%% Quality control
RMS = sqrt(sum((p_lisse-p_mes).^2)/numel(p_mes))

y_ = mean(p_mes);
R2 = sum((p_lisse-y_).^2)/sum((p_mes-y_).^2)


a_approx = 0.5*p_lisse.*vel_mes.^2*S*C_D0/m;

RMS_acc_abs = sqrt(sum((acc_mes-a_approx).^2)/numel(acc_mes))
RMS_acc_rel = sqrt(sum(((acc_mes-a_approx)./acc_mes).^2)/numel(acc_mes))

figure
hold on
plot(t, acc_mes)
plot(t, a_approx)
legend('Mesurée', 'Approximée')
xlabel('Temps (s)')
ylabel('Accélération (m/s^2)')
grid on
saveas(gcf, [pwd '\Figures\acc_mes_vs_approx.png'])



%%

disp('--Identification des paramètres--')
disp([' H_s = ' num2str(hs) newline ' p_0 = ' num2str(p0)])
disp(newline)

figure
subplot(3,1,1)
plot(t,acc_mes)
xlabel('Temps (s)')
ylabel('Accélération (m/s^2)')
grid on
subplot(3,1,2)
plot(t,vel_mes)
xlabel('Temps (s)')
ylabel('Vitesse (m/s)')
grid on
subplot(3,1,3)
plot(t,pos_mes)
xlabel('Temps (s)')
ylabel('Altitude (m)')
grid on
saveas(gcf, [pwd '\Figures\acc_vel_pos.png'])

function [b, m] = moindreCarre(xn, yn)
    A = [numel(xn) sum(xn); sum(xn) sum(xn.^2)];
    B = [sum(yn); sum(xn.*yn)];
    rep = inv(A)*B;
    b = rep(1);
    m = rep(2);
end

function F = trapezoidale(xn,yn)
    h = 0.5;
    if numel(xn) == 1
       F = 0;
    elseif numel(xn) == 2;
       F = (yn(1)+yn(2))*h/2;
    else
       F = 0.5*h*(yn(1) + yn(end) + 2*sum(yn(2:end-1)));
    end
end

function err = integError(xn,yn)
    h = 0.5;
    dy = diff(yn);
    err = abs(h^2*(dy(end)-dy(1))/12);
end
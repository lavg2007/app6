%% simulation de APP6
clear all
close all
clc
cte;
hs = 1.100353042442160e+04;
p0 = 0.021571989401399;
%% simulation comparaison avec les Russes
% tspan = [0, 100];
% z0_1 = [6100, -90,120000,0,-90,0,0]';
% reltol1 = 1e-08;
% options = odeset('abstol' ,1e-06, 'reltol', reltol1);
% [TOUT,YOUT] = ode45('eqny_no',tspan,z0_1,options);
% figure()
% plot(TOUT,YOUT(:,3))
% legend('v','y','h','s','theta','q')
%% simulation avec compensation pour yref nominale
tspan = [0, 200];
z0_1 = [6100, deg2rad(-20.5),120000,0,deg2rad(-80),0,0]';
reltol1 = 1e-08;
options = odeset('abstol' ,1e-06, 'reltol', reltol1);
[TOUT,YOUT] = ode45('eqnynominale',tspan,z0_1,options);
% %%
leg = ['v     ';'y     ';'h     ';'s     ';'\theta';'q     ';'daero '];

figure()
subplot(3,2,1)
plot(TOUT,YOUT(:,2))
legend('\gamma')
xlabel('Temps')
ylabel('gamma')

subplot(3,2,2)
plot(YOUT(:,3),YOUT(:,1))
legend('v')
xlabel('h')
ylabel('v')

subplot(3,2,3)
hold on
plot(TOUT,YOUT(:,5))
plot(TOUT,YOUT(:,5)-YOUT(:,2))
legend('\theta','\alpha')
xlabel('h')
ylabel('v')

subplot(3,2,4)
plot(TOUT,YOUT(:,6))
legend('\dot{\theta}')
xlabel('Temps')
ylabel('vitesse angulaire')

Pdyn =  0.5 .* p0 .* exp(-YOUT(:,3)./hs).* YOUT(:,1).^2;  
subplot(3,2,5)
hold on
plot(TOUT,Pdyn.*S.*Cla.*(YOUT(:,5)-YOUT(:,2)))
plot(TOUT,Pdyn.*S.*Cdo)
legend('Portance','Trainee')
xlabel('Temps')
ylabel('Force(N)')


subplot(3,2,6)
plot(TOUT,YOUT(:,3))
xlabel('Temps')
ylabel('dTlim')

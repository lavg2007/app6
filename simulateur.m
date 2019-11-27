%% simulation 1

disp('Simulation 1...')  

tspan = [0 100];
[t1 y1] = ode45('funcODE1', tspan, cond_ini);
%%

figure
subplot(3,1,1)
plot(t1,rad2deg(y1(:,2)))
xlabel('Temps (s)')
ylabel('\gamma (deg)')
title('Angles de la capsule')
grid on

subplot(3,1,2)
hold on
plot(t1,rad2deg(y1(:,5)))
plot(t1,rad2deg(y1(:,5) - y1(:,2)))
xlabel('Temps (s)')
ylabel('Angle (deg)')
legend('\theta','\alpha')
grid on

subplot(3,1,3)
plot(t1,rad2deg(y1(:,6)))
xlabel('Temps (s)')
ylabel('q (deg/s)')
grid on
saveas(gcf, [pwd '\Figures\sim1_angles.png'])

figure
plot(y1(:,3)/1000,y1(:,1))
xlabel('Altitude (km)')
ylabel('Vitesse (m/s)')
title('Vitesse de la capsule')
grid on
saveas(gcf, [pwd '\Figures\sim1_v_h.png'])

Pdyn_sim1 = 0.5*rho(y1(:,3)).*y1(:,1).^2;
Daero_sim1 = Pdyn_sim1*S*C_D0;
Laero_sim1 = Pdyn_sim1.*S.*C_La.*(y1(:,5)-y1(:,2));

figure
subplot(2,1,1)
hold on
plot(t1,Laero_sim1)
plot(t1,Daero_sim1)
xlabel('Temps (s)')
ylabel('Force (N)')
title('Forces agissant sur la capsule')
legend('Portance (L_{aero})', 'Trainée (D_{aero})')
grid on

subplot(2,1,2)
plot(t1,y1(:,7))
xlabel('Temps (s)')
ylabel('\Delta t_{lim} (s)')
title('\Delta t_{lim} à D_{aero} > 2000 N')
grid on
saveas(gcf, [pwd '\Figures\sim1_force_tlim.png'])


%% simulation 2

disp('Simulation 2...')

[t2 y2] = ode45('funcODE2', tspan, cond_ini);



%%
tspan = [0 110];
hspan = [10 120];


figure
subplot(2,1,1)
plot(t2,rad2deg(y2(:,2)))
line(tspan, [rad2deg(gamma_ref1), rad2deg(gamma_ref1)])
xlim(tspan)
ylim([-25 -10])
xlabel('Temps (s)')
ylabel('\gamma (deg)')
title('Angle de vol')
grid on

subplot(2,1,2)
plot(y2(:,3)/1000,rad2deg(y2(:,2)))
xlabel('Altitude (km)')
ylabel('\gamma (deg)')
xlim(hspan)
ylim([-25 -10])
grid on
saveas(gcf, [pwd '\Figures\sim2_gamma.png'])

figure
subplot(2,1,1)
plot(t2,y2(:,1))
xlabel('Temps (s)')
ylabel('Vitesse (m/s)')
xlim(tspan)
title('vitesse')
grid on

subplot(2,1,2)
plot(y2(:,3)/1000,y2(:,1))
xlabel('Altitude (km)')
ylabel('Vitesse (m/s)')
xlim(hspan)
grid on
saveas(gcf, [pwd '\Figures\sim2_v.png'])

figure
plot(t2,y2(:,3)/1000)
xlabel('Temps (s)')
ylabel('Altitude (km)')
title('Altitude de la capsule')
xlim(tspan)
grid on
saveas(gcf, [pwd '\Figures\sim2_h.png'])

figure
subplot(2,1,1)
hold on
plot(t2,rad2deg(y2(:,5)))
plot(t2,rad2deg(y2(:,5)-y2(:,2)))
xlim(tspan)
ylim([-90 90])
xlabel('Temps (s)')
ylabel('Angle (deg)')
legend('\theta','\alpha')
title('Angle d''attaque et de tanguage')
grid on

subplot(2,1,2)
plot(t2, rad2deg(y2(:,6)))
xlabel('Temps (s)')
ylabel('Vitesse angulaire (deg)')
ylim([-90 90])
xlim(tspan)
grid on
saveas(gcf, [pwd '\Figures\sim2_angles.png'])

Pdyn_sim2 = 0.5*rho(y2(:,3)).*y2(:,1).^2;
Daero_sim2 = Pdyn_sim2*S*C_D0;
Laero_sim2 = Pdyn_sim2.*S.*C_La.*(y2(:,5)-y2(:,2));

figure
subplot(2,1,1)
hold on
plot(t2,Laero_sim2)
plot(t2,Daero_sim2)
xlabel('Temps (s)')
ylabel('Force (N)')
legend('Portance (L_{aero})', 'Trainée (D_{aero})')
title('Forces agissant sur la capsule')
xlim(tspan)
grid on

subplot(2,1,2)
plot(t2,y2(:,7))
xlabel('Temps (s)')
title('\Delta t_{lim} à D_{aero} > 2000 N')
xlim(tspan)
grid on

saveas(gcf, [pwd '\Figures\sim2_forces.png'])


%% simulation 3


disp('Simulation 3...')
disp(newline)

tspan = [0 125];
[t3 y3] = ode45('funcODE3', tspan, cond_ini);

hspan = [10 120];


figure
subplot(2,1,1)
plot(t3,rad2deg(y3(:,2)))
line(tspan, [rad2deg(gamma_ref2), rad2deg(gamma_ref2)])
xlim(tspan)
ylim([-25 -10])
xlabel('Temps (s)')
ylabel('\gamma (deg)')
title('Angle de vol')
grid on

subplot(2,1,2)
plot(y3(:,3)/1000,rad2deg(y3(:,2)))
xlabel('Altitude (km)')
ylabel('\gamma (deg)')
xlim(hspan)
ylim([-25 -10])
grid on
saveas(gcf, [pwd '\Figures\sim3_gamma.png'])

figure
subplot(2,1,1)
plot(t3,y3(:,1))
xlabel('Temps (s)')
ylabel('Vitesse (m/s)')
xlim(tspan)
title('vitesse')
grid on

subplot(2,1,2)
plot(y3(:,3)/1000,y3(:,1))
xlabel('Altitude (km)')
ylabel('Vitesse (m/s)')
xlim(hspan)
grid on
saveas(gcf, [pwd '\Figures\sim3_v.png'])

figure
plot(t3,y3(:,3)/1000)
xlabel('Temps (s)')
ylabel('Altitude (km)')
title('Altitude de la capsule')
xlim(tspan)
grid on
saveas(gcf, [pwd '\Figures\sim3_h.png'])

figure
subplot(2,1,1)
hold on
plot(t3,rad2deg(y3(:,5)))
plot(t3,rad2deg(y3(:,5)-y3(:,2)))
xlim(tspan)
ylim([-90 90])
xlabel('Temps (s)')
ylabel('Angle (deg)')
legend('\theta','\alpha')
title('Angle d''attaque et de tanguage')
grid on

subplot(2,1,2)
plot(t3, rad2deg(y3(:,6)))
xlabel('Temps (s)')
ylabel('Vitesse angulaire (deg)')
ylim([-90 90])
xlim(tspan)
grid on
saveas(gcf, [pwd '\Figures\sim3_angles.png'])

Pdyn_sim3 = 0.5*rho(y3(:,3)).*y3(:,1).^2;
Daero_sim3 = Pdyn_sim3*S*C_D0;
Laero_sim3 = Pdyn_sim3.*S.*C_La.*(y3(:,5)-y3(:,2));

figure
subplot(2,1,1)
hold on
plot(t3,Laero_sim3)
plot(t3,Daero_sim3)
xlabel('Temps (s)')
ylabel('Force (N)')
legend('Portance (L_{aero})', 'Trainée (D_{aero})')
title('Forces agissant sur la capsule')
xlim(tspan)
grid on

subplot(2,1,2)
plot(t3,y3(:,7))
xlabel('Temps (s)')
title('\Delta t_{lim} à D_{aero} > 2000 N')
xlim(tspan)
grid on

saveas(gcf, [pwd '\Figures\sim3_forces.png'])


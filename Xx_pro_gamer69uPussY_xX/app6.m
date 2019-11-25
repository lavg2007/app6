%% problematique APP6
close all
clear all
clc
%% Selection des figures
% # 1 Accelration de la sonde Russe
% # 2 Vitesse de la sonde Russe
% # 3 Postioion de la sonde Russe
% # 4 Vitesse de la sonde Russe et fct de la position
% # 5 Comparaison des accelerations
% # 6 Verfication de RAA acceleration
% # 7 Verfication de RAA vitesse
% # 8 Gamma ref en fonction de l'altitude
% # 9 Comparaison des vitesses
% # 10 Trouver dTlim
    name10 ='Trouver dTlim'; 
% # 11 Angle de commande en fct de h
%       [1 2 3 4 5 6 7 8 9 10 11];
plots = [0 0 0 0 0 0 0 0 0 0 0];
%% Identification des donnees des Russes
addpath('.\functions')
load('Accelero_Data_from_Moscow')
cte;
russes.acc.val = acc_mes;
clear acc_mes;
% plot des valeurs
%%%%%%%%% FIGURE %%%%%%%%%
if plots(1)
    figure(1)
    plot(t,russes.acc.val,'-x')
    title('Acceleration de la sonde Russe')
    saveas(gcf,'.\figures\Acceleration_de_la_sonde_Russe.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% obtenir les valeurs de vitesse
% [simpson_result,simpson_error] = simpson(t,russes.acc.val);
% [trapeze_result,trapeze_error] = trapeze(t,russes.acc.val) ;
    % la methode des trapeze donne une plus petite erreur, elle sera donc
    % utilise

% vitesse en tout point
for i = 2:numel(russes.acc.val)
    [russes.vit.val(i),russes.vit.error] = trapeze(t(1:i),russes.acc.val(1:i)');
    russes.vit.val(i) = vinit -russes.vit.val(i);
end
russes.vit.val(1) = vinit;
%%%%%%%%% FIGURE %%%%%%%%%
if plots(2)
    figure(2)
    hold on
    plot(t,russes.vit.val)
    title('Vitesse de la sonde Russe')
    saveas(gcf,'.\figures\Vitesse_de_la_sonde_Russe.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% obtenir les valeurs de position
% [simpson_result,simpson_error] = simpson(t,vitesse);
% [trapeze_result,trapeze_error] = trapeze(t,vitesse) ;
    % la methode des simpson donne une plus petite erreur, mais la methodes
    % des trapezes represente mieux les donnees
    

% vitesse en tout point
for i = 2:numel(russes.vit.val)
    [russes.h.val(i),russes.h.error] = trapeze(t(1:i),russes.vit.val(1:i));
    russes.h.val(i) =  hinit - russes.h.val(i);
end
russes.h.val(1) = hinit;
%%%%%%%%% FIGURE %%%%%%%%%
if plots(3)
    figure(3)
    hold on
    plot(t,russes.h.val)
    title('Position de la sonde Russe')
    saveas(gcf,'.\figures\Position_de_la_sonde_Russe.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FIGURE %%%%%%%%%
if plots(4)
    figure(4)
    hold on
    plot(russes.h.val,russes.vit.val)
    title('Vitesse de la sonde Russe en fct de la position')
    saveas(gcf,'.\figures\Vitesse_de_la_sonde_Russe_en_fct_de_la_position.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %% identification des parametres
Y = log(russes.acc.val'./russes.vit.val.^2);
X = russes.h.val;
[mc.coeff,mc.error] = MC(X,Y,2);
% traduction des coeff
p0 = exp(mc.coeff(1)) / ((S*Cdo) / (2*m));
hs = -1/mc.coeff(2);
    %% verification
verif.iden.accel.val = 1./(2*m) .* p0.*exp(-russes.h.val/hs) .* russes.vit.val.^2 *S *Cdo;
verif.iden.accel.error.rmsabs = sqrt(1/numel(russes.acc.val)*sum((verif.iden.accel.val - russes.acc.val').^2));
verif.iden.accel.error.rmsrela = sqrt(1/numel(russes.acc.val)*sum(((verif.iden.accel.val - russes.acc.val')./russes.acc.val').^2));
%%%%%%%%% FIGURE %%%%%%%%%
if plots(5)
    figure(5)
    hold on
    plot(t,russes.acc.val)
    plot(t,verif.iden.accel.val)
    title('Verification de l''acc')
    legend('Russes','Iden')
    saveas(gcf,'.\figures\Verification_de_la_acc.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loi de guidage
% RAA
pfin = p0 .* exp(-10000./hs); % derniere valeur de densite
p = p0 .* exp(-russes.h.val./hs);
    %% verification de la precision de la RAA
verif.RAA.vitesse.val = vinit*exp(B*0.5*hs*(p-p(1))./sin(-pi/2));
verif.RAA.accel.val = 0.5.*p.*verif.RAA.vitesse.val.^2 .* S.*Cdo./m;
verif.RAA.accel.R2 = sum((verif.RAA.accel.val-mean(russes.acc.val)).^2)/sum((russes.acc.val-mean(russes.acc.val)).^2);
verif.RAA.accel.RMS = sqrt(1/numel(russes.acc.val)*sum((verif.RAA.accel.val - russes.acc.val').^2));
%%%%%%%%% FIGURE %%%%%%%%%
if plots(6)
    figure(6)
    hold on
    plot(t,russes.acc.val','-x')
    plot(t,verif.RAA.accel.val,'-x')
    title('Verification de RAA acceleration')
    legend('Russe','RAA')
    saveas(gcf,'.\figures\Verification_RAA_Acc.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%% FIGURE %%%%%%%%%
if plots(7)
    figure(7)
    hold on
    plot(t,russes.vit.val','-x')
    plot(t,verif.RAA.vitesse.val,'-x')
    title('Verification de RAA vitesse')
    legend('Russe','RAA')
    saveas(gcf,'.\figures\Verification_RAA_Vit.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% trouver yref
% Changement de vitesse necessaire pour avoir la vitesse finale
dvaero.nominale = vfin.nominale - sqrt(russes.vit.val.^2 +2*um*(1/rfin - 1./(r+russes.h.val)));
dvaero.ideale = vfin.ideale - sqrt(russes.vit.val.^2 +2*um*(1/rfin - 1./(r+russes.h.val)));

yref.nominale = asin(B*0.5*hs*(pfin - p)./log(1+(dvaero.nominale./russes.vit.val)));
yref.ideale = asin(B*0.5*hs*(pfin - p)./log(1+(dvaero.ideale./russes.vit.val)));
% j'uitlise juste les 41 premieres valeurs pcq apres ca on descends en
% dessous de la position hfin de 10 000m en altitude
%%%%%%%%% FIGURE %%%%%%%%%
if plots(8)
    figure(8)
    hold on
    plot(russes.h.val(1:41),rad2deg(yref.nominale(1:41)) ,'-x')
    plot(russes.h.val(1:41),rad2deg(yref.ideale(1:41)) ,'-x')
    title('Angle de commande en fct de h')
    legend('nominale','ideale')
    saveas(gcf,'.\figures\Angle_de_commande_en_fct_de_h.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 
verif.guidage.vitesse.nominale.val = vinit*exp(B*0.5*hs*(p(1:41)-p(1))./sin(yref.nominale(1)));
verif.guidage.vitesse.ideale.val = vinit*exp(B*0.5*hs*(p(1:41)-p(1))./sin(yref.ideale(1)));
%%%%%%%%% FIGURE %%%%%%%%%
if plots(9)
    figure(9)
    hold on
    plot(russes.h.val(1:41),verif.guidage.vitesse.nominale.val ,'-x')
    plot(russes.h.val(1:41),verif.guidage.vitesse.ideale.val ,'-x')
    plot(russes.h.val(1:41),russes.vit.val(1:41),'-o')
    title('Verification de vitesse avec yref')
    legend('nominale','ideale','Russe')
    saveas(gcf,'.\figures\vitesse_et_yref.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Calcul de Daero
    Daero.positions = [120000:-10:10000];
    Daero.pinit = p0 .* exp(-120000./hs);

    Daero.p = p0*exp(-Daero.positions/hs) ;                                                    %ok
    Daero.v = vinit*exp(B*0.5*hs*(Daero.p-Daero.pinit)./sin(yref.nominale(1)));          %ok
    Daero.dp_dh = -p0/hs*exp(-Daero.positions./hs);                                           %ok
    Daero.dv_dh = (-B*vinit*p0)/(2*sin(yref.nominale(1))) * exp((B*hs/(2*sin...
        (yref.nominale(1))))*(p0*exp(-Daero.positions./hs)-Daero.pinit)-Daero.positions./hs );             %ok
    Daero.der = 0.5.*S.*Cdo.* Daero.v.*(Daero.v.*Daero.dp_dh+2.*Daero.p.*Daero.dv_dh);   %ok
    Daero.val = 0.5*Daero.p.*Daero.v.^2 * S*Cdo;
%ok
%     figure(100)
%     plot(Daero.positions,Daero.val)
    %% Determination de dlim
% ici il va falloir faire NR afin de trouver les point ou Daero est = 2000N
% il faut donc construire la derivee de Daero, ce qui est assez complexe...
% dp_dh = -p0/hs*exp(-russes.h.val./hs);
% dv_dh = (-B*vinit*p0)/(2*sin(yref.nominale(1))) * exp((B*hs/(2*sin(yref.nominale(1))))*(p0*exp(-russes.h.val./hs)-p(1))-russes.h.val./hs )
% der_Daero = 0.5.*S.*Cdo.*dp_dh(1:41).*verif.guidage.vitesse.nominale.val.^2 + p(1:41) .* 2.*dv_dh(1:41)
i = 0;
k =0;
for y = [yref.nominale(1),yref.ideale(1)]    
for search = [5.6e4,1.585e4]
    i= i+1;
    dp = p0*exp(-search/hs);                                               
    v = vinit*exp(B*0.5*hs*(dp-Daero.pinit)./sin(y));         
    dp_dh = -p0/hs*exp(-search./hs);                                          
    dv_dh = (-B*vinit*p0)/(2*sin(y)) * exp((B*hs/(2*sin...
        (y)))*(p0*exp(-search./hs)-Daero.pinit)-search./hs );        
    der = 0.5.*S.*Cdo.* v.*(v.*dp_dh+2.*dp.*dv_dh); 
    val = 0.5.*dp.*v.^2 * S*Cdo-2000;                  
                                      
    pos = search;
    while abs(val) > 1e-8
        k = k+1;
        p = p0*exp(-pos/hs);                                                     %ok
        v = vinit*exp(B*0.5*hs*(p-Daero.pinit)./sin(y));          %ok
        dp_dh = -p0/hs*exp(-pos./hs);                                           %ok
        dv_dh = (-B*vinit*p0)/(2*sin(y)) * exp((B*hs/(2*sin...
            (y)))*(p0*exp(-pos./hs)-Daero.pinit)-pos./hs );             %ok
        der = 0.5.*S.*Cdo.* v.*(v.*dp_dh+2.*p.*dv_dh);   %ok
        val = 0.5*p.*v.^2 * S*Cdo-2000 ;                   %ok
        pos = pos - val/der;
    end
Daero.tlim.iter(i) = k;
Daero.tlim.h(i) = pos;
Daero.tlim.v(i) = v;
end
end
% delta position
Daero.dhtlim.nominale = abs(Daero.tlim.h(1)- Daero.tlim.h(2));
Daero.dhtlim.ideale = abs(Daero.tlim.h(3)- Daero.tlim.h(4));
% Delta vitesse
Daero.dvtlim.nominale = abs(Daero.tlim.v(1)- Daero.tlim.v(2));
Daero.dvtlim.ideale = abs(Daero.tlim.v(3)- Daero.tlim.v(4));
% Delta temsps
Daero.dtlim.nominale = Daero.dhtlim.nominale/Daero.dvtlim.nominale;
Daero.dtlim.ideale = Daero.dhtlim.ideale/Daero.dvtlim.ideale;
%%%%%%%%% FIGURE %%%%%%%%%
if plots(10)
    figure(10)
    hold on
    line([1e4 6e4],[2e3 2e3])
    plot(Daero.positions,Daero.val)
    plot(Daero.tlim.h,2000,'or')
    legend('Daero nominale slm')
    title(name10)
    saveas(gcf,'.\figures\Daero_nominale_slm.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %% COMMANDE DYNAMIQUE
% voir eqny
%% simulation comparaison avec les Russes

tspan = [0, 30];
z0_1 = [6100, -90,120000,0,-90,0,0]';
reltol1 = 1e-08;
options = odeset('abstol' ,1e-06, 'reltol', reltol1);
[TOUT,YOUT] = ode45('eqny_no',tspan,z0_1,options);
% if plots(11)
%     figure(11)
%     plot( ,'-x')
%     title('Angle de commande en fct de h')
%     legend('nominale','ideale')
%     saveas(gcf,'.\figures\Angle_de_commande_en_fct_de_h.png')
% end
    %% affichage
disp('**---** IDENTIFICATION **---**')
disp('---- Error intergration ----')
disp(['vitesse = ' num2str(russes.vit.error)])
disp(['position = ' num2str(russes.h.error)])

disp('---- Identification ----')
disp(['p0 = ' num2str(p0)])
disp(['hs = ' num2str(hs)])

disp('---- Error identification ----')
disp(['R2 = ' num2str(mc.error.R2)])
disp(['RMS = ' num2str(mc.error.rms)])
% on observe que l'erreur RMS est tres pres de l'erreur de mesure des
% Russes. On confirme donc que la methode utilisee est valide
disp('---- Verification des erreur d''acc ----')
disp(['RMS de l''accel abs = ' num2str(verif.iden.accel.error.rmsabs)])
disp(['RMS de l''accel rela = ' num2str(verif.iden.accel.error.rmsrela)])
disp('')
disp('**---** LOI DE GUIDAGE **---**')
disp('---- Verification des erreur d''acc par RAA ----')
disp(['R2 = ' num2str(verif.RAA.accel.R2 )])
disp(['RMS = ' num2str(verif.RAA.accel.RMS)])
disp('---- dTlim ----')
disp(['       Ideale      Nominale'])
disp(['vfin = ' num2str(vfin.ideale) '         ' num2str(vfin.nominale) ])
disp(['yref = ' num2str(rad2deg(yref.ideale(1))) '   ' num2str(rad2deg(yref.nominale(1)))])
disp(['hmin = ' num2str(min([Daero.tlim.h(3) Daero.tlim.h(4)])) '  ' num2str(min([Daero.tlim.h(1) Daero.tlim.h(2)]))])
disp(['vmin = ' num2str(min([Daero.tlim.v(3) Daero.tlim.v(4)])) '    ' num2str(min([Daero.tlim.v(1) Daero.tlim.v(2)]))])
disp(['h/# =  ' num2str(1.585e4) '\' num2str(Daero.tlim.iter(4)) '    ' num2str(1.585e4) '\' num2str(Daero.tlim.iter(2))])
disp(['hmax = ' num2str(max([Daero.tlim.h(3) Daero.tlim.h(4)])) '  ' num2str(max([Daero.tlim.h(1) Daero.tlim.h(2)]))])
disp(['vmax = ' num2str(max([Daero.tlim.v(3) Daero.tlim.v(4)])) '   ' num2str(max([Daero.tlim.v(1) Daero.tlim.v(2)]))])
disp(['h/# =  ' num2str(5.6e4) '\' num2str(Daero.tlim.iter(3)) '    ' num2str(5.6e4) '\' num2str(Daero.tlim.iter(1))])
disp(['dtlim = ' num2str(Daero.dtlim.ideale) '     ' num2str(Daero.dtlim.nominale)])
disp('')
disp('**---** COMMANDE DYNAMIQUE DE TRANSLATION **---**')


% close all
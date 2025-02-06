%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P2in = [29.70; 29.72; 29.62; 29.63; 29.61; 29.64; 29.63; 29.57;...
 29.56; 29.64; 29.64; 29.64];
P2mm = [754.7; 754.6; 752.8; 752.4; 753.6; 752.8; 752.5; 751.0; 751.2;...
 752.7; 752.9; 752.4];
T2F = [73.1; 73.5; 73.2; 73.9; 73.7; 73.8; 74.1; 74.0; 74.2; 74.1;...
 74.2; 74.4];
P1in = [30.00; 29.70; 29.80; 29.80; 29.66; 29.81; 29.83; 29.73; 29.71;...
 29.74; 29.76; 29.71; 29.72; 29.70; 29.72; 29.73; 29.73; 29.73; 29.74];
P3in = [29.71; 29.68; 29.65; 29.71; 29.69; 29.70; 29.69; 29.70; 29.69;...
 29.70; 29.69; 29.68; 29.75; 29.72; 29.65; 29.65; 29.63; 29.71; 29.67];
n1 = numel(P1in);
n2 = numel(P2in);
n3 = numel(P3in);
Ni = n1 + n2 + n3;
%%%%%%%%%%%%%%%%%%%%%%%%%% Mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MP2in = sum(P2in)/12; %29.6333
MP2mm = sum(P2mm)/12/25.4; %752.8000 mm %29.6378
MT2F = sum(T2F)/12; %73.8500
%%%%%%%%%%%%%%%%%%% Offsets/Interpolation Inches %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Gravity Correction Inches %%%%%
x1g = 29; x2g = 30;
y1g = 32; y2g = 34;
f11g = 0.035; f21g = 0.036; f12g = 0.030; f22g = 0.031;
xig=MP2in; yig=32.7157;
Fig = (yig-y2g)/(y1g-y2g)*((xig-x2g)/(x1g-x2g)*f11g + ...
 (xig-x1g)/(x2g-x1g)*f21g) + (yig-y1g)/(y2g-y1g)*((xig-x2g)/...
 (x1g-x2g)*f12g + (xig-x1g)/(x2g-x1g)*f22g);
American Institute of Aeronautics and Astronautics
9
% Fig = 0.0338
%%%%% Tempreture Correction Inches %%%%%
x1 = 29; x2 = 30;
y1 = 72; y2 = 74;
f11 = 0.114; f21 = 0.118; f12 = 0.119; f22 = 0.123;
xi = MP2in; yi = MT2F;
Fi = (yi-y2)/(y1-y2)*((xi-x2)/(x1-x2)*f11 + (xi-x1)/(x2-x1)*f21) + ...
 (yi-y1)/(y2-y1)*((xi-x2)/(x1-x2)*f12 + (xi-x1)/(x2-x1)*f22);
% Fi = 0.1212
%%%%%%%%%%%%%%%%%%%% Applying Corrections Inches %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inches %%%%%
cgin = (32.1740 - Fig)*0.3048; % Corrected Gravity = 9.7963 m/s^2
chin = (P2in - Fi)/ 39.37; % Corrected Height in m
rho = 13595; % Density of Mercury in kg/m^3
cPin = rho*cgin*chin; % Corrected Pressure Pa
%%%%%%%%%%%%%%%%%%% Offsets/Interpolation MM %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Gravity Correction MM %%%%%
X1g = 29; X2g = 30;
Y1g = 32; Y2g = 34;
F11g = 0.035; F21g = 0.036; F12g = 0.030; F22g = 0.031;
Xig=MP2mm; Yig=32.7157;
Figg = (Yig-Y2g)/(Y1g-Y2g)*((Xig-X2g)/(X1g-X2g)*F11g + ...
 (Xig-X1g)/(X2g-X1g)*F21g) + (Yig-Y1g)/(Y2g-Y1g)*((Xig-X2g)/...
 (X1g-X2g)*F12g + (Xig-X1g)/(X2g-X1g)*F22g);
% Fig = 0.0338
%%%%% Tempreture Correction MM %%%%%
X1 = 29; X2 = 30;
Y1 = 72; Y2 = 74;
F11 = 0.114; F21 = 0.118; F12 = 0.119; F22 = 0.123;
Xi = MP2mm; Yi = MT2F;
Fii = (Yi-Y2)/(Y1-Y2)*((Xi-X2)/(X1-X2)*F11 + (Xi-X1)/(X2-X1)*F21) + ...
 (Yi-Y1)/(Y2-Y1)*((Xi-X2)/(X1-X2)*F12 + (Xi-X1)/(X2-X1)*F22);
% Fi = 0.1212
%%%%%%%%%%%%%%%%%%%% Applying Corrections MM %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inches %%%%%
cgmm = (32.1740 - Figg)*0.3048; % Corrected Gravity = 9.7963 m/s^2
chmm = (P2in - Fii)/39.37; % Corrected Height in m
rho = 13595; % Density of Mercury in kg/m^3
cPmm = rho*cgmm*chmm; % Corrected Pressure in
%%%%%%%%%%%%%%%%%%%%% Corrected Data Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 2.201;
R = 287.058;
T = (T2F-32)*(5/9)+273.15;
MT = mean(T);
SMT = sqrt((1/(12-1))*sum((T-MT).^2));
UT = 2.201*SMT;
SDTM = SMT/sqrt(12);
USDTM = t*SDTM;
D = cPin./(R*((T2F-32)*(5/9)+273.15));
Table = table(P2in,P2mm,T2F,cPin, cPmm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Table 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTOTS = mean(T2F);
SXTOTS = sqrt((1/(12-1))*sum((T2F-MTOTS).^2));
UTOTS = t*SXTOTS;
MPINOTS = mean(P2in);
SXPNOTS = sqrt((1/(12-1))*sum((P2in-MPINOTS).^2));
UPINOTS = t*SXPNOTS;
MPPOTS = mean(cPin);
SXPPOTS = sqrt((1/(12-1))*sum((cPin-MPPOTS).^2));
UPPOTS = t*SXPPOTS;
MDOTS = mean(D);
SXDOTS = sqrt((1/(12-1))*sum((D-MDOTS).^2));
UDOTS = sqrt(((1/(R^2*MT^2))*UPPOTS^2)+((MPPOTS^2/(R^2*MT^4))*UT^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Table 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SDMTE = SXTOTS/sqrt(12);
USDMTE = t*SDMTE;
SDMPINE = SXPNOTS/sqrt(12);
USDMPINE = t*SDMPINE;
SDMPPE = SXPPOTS/sqrt(12);
USDMPPE = t*SDMPPE;
SDMDE = MDOTS/sqrt(12);
UDE = sqrt(((1/(R^2*MT^2))*USDMPPE^2)+((MPPOTS^2/(R^2*MT^4))*USDTM^2));
%%%%%%%%%%%%%%%%%%%%%%%%%% Convergence Test %%%%%%%%%%%%%%%%%%%%%%%%%
American Institute of Aeronautics and Astronautics
11
Pressure1 = (P1in)';
Pressure2 = (P2in)';
Pressure3 = (P3in)';
Pressure = [Pressure1 Pressure2 Pressure3];
Mean(1,1) = 29.70;
Standard(1,1)= 0;
Standard_Mean(1,1) = 0;
for i = 2:Ni
 Mean(1,i) = ((sum(Pressure(1,1:i)))/i);
 Standard(1,i) = sqrt((1/(i-1))*sum((Pressure(1,i)-Mean(1,i)).^2));
 Standard_Mean(1,i) = Standard(1,i)./sqrt(i);
end
figure (1)
subplot(3,1,1)
plot(1:1:Ni,Mean(1:Ni), 'color', "#0072BD");
title('Convergence of Mean')
subtitle(['Measured Pressure [inHg] --- Convergance Towards $\bar{x}$' ...
 ' = 29.70'],'Interpreter','Latex')
ylabel('$\bar{x}$ [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:Ni])
grid on
subplot(3,1,2)
plot(2:1:Ni,Standard(2:Ni), 'color', "#A2142F");
title('Convergence of Standard Deviation')
subtitle('Measured Pressure [inHg] --- Convergance Towards Sx = 0', ...
 'Interpreter','Latex')
ylabel('Sx [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:Ni])
grid on
subplot(3,1,3)
plot(2:1:Ni,Standard_Mean(2:Ni), 'color', "#77AC30");
title('Convergence of Standard Deviation of the Means')
subtitle(['Measured Pressure [inHg] --- Convergance ' ...
 'Towards S$\bar{x}$= 0'], 'Interpreter','Latex')
ylabel('S$\bar{x}$ [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:Ni])
grid on
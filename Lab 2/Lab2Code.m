clc
clearvars
%%
%Data
P_0_Raw_Data = readtable('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '0 in H20 Test', 'Range', 'B3:AJ602');
P_2_Raw_Data = readtable('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '2 in H20 Test', 'Range', 'B3:AJ602');
P_5_Raw_Data = readtable('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '5 in H20 Test', 'Range', 'B3:AJ602');

P_0 = table2array(mean(P_0_Raw_Data, 1));
P_2 = table2array(mean(P_2_Raw_Data, 1));
P_5 = table2array(mean(P_5_Raw_Data, 1));
P_Ambient = readmatrix('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '0 in H20 Test', 'Range', 'G1:G1'); % inHg

T_0 = readmatrix('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '0 in H20 Test', 'Range', 'P1:P1'); %F
T_2 = readmatrix('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '2 in H20 Test', 'Range', 'P1:P1'); %F
T_5 = readmatrix('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '5 in H20 Test', 'Range', 'P1:P1'); %F
T_Ambient = readmatrix('Lab 2 Data 2025SP_UPDATED.xlsx', 'Sheet', '0 in H20 Test', 'Range', 'K1:K1'); %F

R = 287.058;
%%
%Corrections
% Latitude Correction
SD_lat = 32.7157;
Lat_correction29 = interp1([32 34],[0.035 0.030], SD_lat);
Lat_correction30 = interp1([32 34],[0.036 0.031], SD_lat);
Lat_correction = interp1([29 30], [Lat_correction29 Lat_correction30], P_Ambient, 'linear', 'extrap');

% Temperature Correction

temps = 71:2:75;
correction29 = [.114, .119, .124];
correction30 = [.118, .123, .128];
correction31 = [.122, .127, .133];

T_correct29 = interp1(temps, correction29, T_Ambient);
T_correct30 = interp1(temps, correction30, T_Ambient);
T_correct31 = interp1(temps, correction31, T_Ambient);
Temp_correction = interp1([29 30 31], [T_correct29 T_correct30 T_correct30], P_Ambient);

% Total inHg correction
Corrected_h = P_Ambient - Lat_correction - Temp_correction;
Corrected_h = Corrected_h * 0.0254;

% Corrected Pressure Calculation
Standard_gravity = 9.80665; % m/s^2
rhoHG = 13595; % kg/m^3
Corrected_Pressure_Amb = rhoHG.*Standard_gravity*Corrected_h; % kg/m^3 * m/s^2 * m = Pa
Corrected_Pressure_Amb = Corrected_Pressure_Amb / 6894.76;
Temps = [T_0; T_2; T_5];

% q indicated
q_indicated0 = [P_0(1:6)'; P_0(8:35)'];
q_indicated2 = [P_2(1:6)'; P_2(8:35)'];
q_indicated5 = [P_5(1:6)'; P_5(8:35)'];
% q true
q_true0 = q_indicated0 - P_0(7);
q_true2 = q_indicated2 - P_2(7);
q_true5 = q_indicated5 - P_5(7);

% Tunnel flow uniformity
q_uniform0 = (( q_true0 - mean(q_true0, 'all') ) / mean(q_true0, 'all')) * 100;
q_uniform2 = (( q_true2 - mean(q_true2, 'all') ) / mean(q_true2, 'all')) * 100;
q_uniform5 = (( q_true5 - mean(q_true5, 'all') ) / mean(q_true5, 'all')) * 100;
% Corrected Barometer Pressure for Density
P_corr0 = Corrected_Pressure_Amb + P_0(7) - P_0(7);P_corr0 = P_corr0 * 6894.76;
P_corr2 = Corrected_Pressure_Amb + P_2(7) - P_0(7);P_corr2 = P_corr2 * 6894.76;
P_corr5 = Corrected_Pressure_Amb + P_5(7) - P_0(7);P_corr5 = P_corr5 * 6894.76;

%Density Calcs
rho0 = P_corr0 / (R * ((Temps(1) - 32) * 5/9) + 273.15);
rho2 = P_corr2 / (R * ((Temps(2) - 32) * 5/9) + 273.15);
rho5 = P_corr5 / (R * ((Temps(3) - 32) * 5/9) + 273.15);
v0 = (sqrt((2 .* abs(q_true0 * 6894.76)) ./ rho0));
v2 = (sqrt((2 .* abs(q_true2 * 6894.76)) ./ rho2));
v5 = (sqrt((2 .* abs(q_true5 * 6894.76)) ./ rho5));
% Tunnel flow uniformity
v_uniform0 = (( v0 - mean(v0, 'all') ) / mean(v0, 'all')) * 100;
v_uniform2 = (( v2 - mean(v2, 'all') ) / mean(v2, 'all')) * 100;
v_uniform5 = (( v5 - mean(v5, 'all') ) / mean(v5, 'all')) * 100;
%unit Re = rho v / meu
Re0 = (rho0 * v0) / (18.34 * 10^-6);
Re2 = (rho2 * v2) / (18.34 * 10^-6);
Re5 = (rho5 * v5) / (18.6 * 10^-6);
disp('Res')
disp([mean(Re0, 'all'), mean(Re2, 'all'), mean(Re5, 'all')])
%%
% Measured Pressure Uniformity
fig = figure('Visible', 'off');
plot(P_0, 'o');
hold on;
plot(P_2, 'o');
plot(P_5, 'o');
yline([mean(P_0, 'all'), mean(P_2, 'all'), mean(P_5, 'all')], 'r');
grid on; grid minor;
xlabel('Pressure Port');ylabel('Pressure [psi]');
title('Flow uniformity - Measured Pressure');
legend('0.0 inH_2O', '2.0 inH_2O', '5.0 inH_2O', 'Average', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow uniformity_Measured Pressure.pdf', 'ContentType', 'vector');
close(fig);
%%
% Dynamic Pressure Uniformity
fig = figure('Visible', 'off');
plot(q_true0, 'o');
hold on;
plot(q_true2, 'o');
plot(q_true5, 'o');
yline([0.0, 0.0721824, 0.252638], 'k');
yline([mean(q_true0, 'all'), mean(q_true2, 'all'), mean(q_true5, 'all')], 'r');
grid on; grid minor;
xlabel('Pressure Port');
ylabel('Pressure [psi]');
title('Flow Uniformity - Dynamic Pressure');
legend('0.0 inH_2O', '2.0 inH_2O', '5.0 inH_2O', 'Average', '','','q_{setting}', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow Uniformity_Dynamic Pressure.pdf', 'ContentType', 'vector');
close(fig);
%%
% Airspeed Uniformity
fig = figure('Visible', 'off');
plot(v0, 'o');
hold on;
plot(v2, 'o');
plot(v5, 'o');
yline([mean(v0, 'all'), mean(v2, 'all'), mean(v5, 'all')], 'k');
yline([0, 29, 55], 'r');
grid on; grid minor;
xlabel('Pressure Port');
ylabel('Airspeed [m/s]');
title('Flow Uniformity - Test Section Airspeed');
legend('0.0 inH_2O', '2.0 inH_2O', '5.0 inH_2O', 'Average','','', 'V_{setting}', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow Uniformity_Test Section Airspeed.pdf', 'ContentType', 'vector' );
close(fig);
%%
% Dynamic Pressure
fig = figure('Visible', 'off');
plot(q_uniform2, 'ko')
hold on
plot(q_uniform5, 'rd')
grid on; grid minor;
xlabel('Pressure Port');
ylabel('$\frac{q - \bar{q}}{\bar{q}}$','Interpreter', 'latex');
title('Flow Uniformity - Dynamic Pressure Deviation from Average');
legend('2.0 inH_2O', '5.0 inH_2O', Location='best');
ax = gca;
ax.YAxis.TickLabelFormat = '%g%%';
hold off
exportgraphics(fig, 'Flow Uniformity_Dynamic Pressure Deviation from Average.pdf', 'ContentType', 'vector');
close(fig);
%%
% Airspeed
fig = figure('Visible', 'off');
plot(v_uniform2, 'ko')
hold on
plot(v_uniform5, 'rd')
grid on; grid minor;
xlabel('Pressure Port');
ylabel('$\frac{v - \bar{v}}{\bar{v}}$','Interpreter', 'latex');
title('Flow Uniformity - Airspeed Deviation from Average');
legend('2.0 inH_2O', '5.0 inH_2O', Location='best');
ax = gca;
ax.YAxis.TickLabelFormat = '%g%%';
hold off
exportgraphics(fig, 'Flow Uniformity_Airspeed Deviation from Average.pdf', 'ContentType', 'vector');
close(fig);
%%
x = [-19.5000, -16.5000, -13.5000, -10.5000, -7.5000, -4.5000, 0, 1.5000, 4.5000, 7.5000, 10.5000, 13.5000, 16.5000, 19.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000]';
z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.5000, 10.5000, 7.5000, 4.5000, 1.5000, -1.5000, -4.5000, -7.5000, -10.5000, -13.5000, 13.5000, 10.5000, 7.5000, 4.5000, 1.5000, -1.5000, -4.5000, -7.5000, -10.5000, -13.5000]';

%%
% Contour of Dynamic Pressure Deviation for 2in Plot %
fig = figure('Visible', 'off');
dq = q_uniform2;
[X,Z] = meshgrid(linspace(min(x), max(x), 100), linspace(min(z), max(z), 100));
F = scatteredInterpolant(x, z, dq, 'linear', 'none');
Dq = F(X, Z);
contourf(X, Z, Dq);
hold on; grid on; grid minor;
scatter(x, z, 'ro');
xlabel('x position [in]');
ylabel('z position [in]');
h = colorbar; 

ticks = h.Ticks;
percentLabels = cellstr(num2str(ticks', '%.1f%%'));
h.TickLabels = percentLabels; % Set custom tick labels with percentages
set(get(h, 'label'), 'string', '$\frac{q - \bar{q}}{\bar{q}}$','Interpreter', 'latex');
title('Dynamic Pressure Deviation [q=2.0 inH_2O]');
hold off;
exportgraphics(fig, 'Dynamic Pressure Deviation q=2.0 inH_2O.pdf', 'ContentType', 'vector');
close(fig);
%%
% Contour of Dynamic Pressure Deviation for 5in Plot %
fig = figure('Visible', 'off');
dq = q_uniform5;
[X,Z] = meshgrid(linspace(min(x), max(x), 100), linspace(min(z), max(z), 100));
F = scatteredInterpolant(x, z, dq, 'linear', 'none');
Dq = F(X, Z);
contourf(X, Z, Dq);
hold on; grid on; grid minor;
scatter(x, z, 'ro');
xlabel('x position [in]');
ylabel('z position [in]');
h = colorbar; 

ticks = h.Ticks;
percentLabels = cellstr(num2str(ticks', '%.1f%%'));
h.TickLabels = percentLabels; % Set custom tick labels with percentages
set(get(h, 'label'), 'string', '$\frac{q - \bar{q}}{\bar{q}}$','Interpreter', 'latex');
title('Dynamic Pressure Deviation [q=5.0 inH_2O]');
hold off;
exportgraphics(fig, 'Dynamic Pressure Deviation q=5.0 inH_2O.pdf', 'ContentType', 'vector');
close(fig);
%%
% Contour of Airspeed Distribution for 7in Plot, IAS = 54.4 m/s %
fig = figure('Visible', 'off');
dq = v5;
[X,Z] = meshgrid(linspace(min(x), max(x)), linspace(min(z), max(z)));
F = scatteredInterpolant(x, z, dq, 'linear', 'none');
Dq = F(X, Z);
contourf(X, Z, Dq);
hold on; grid on; grid minor;
scatter(x, z, 'ro');
xlabel('x position [in]');
ylabel('z position [in]');
h = colorbar;
set(get(h, 'label'), 'string', 'Measured Airspeed [m/s]');
title('Airspeed Distribution - IAS = 54.5 m/s');
hold off
exportgraphics(fig, 'Airspeed Distribution_54.5.pdf', 'ContentType', 'vector');
close(fig);
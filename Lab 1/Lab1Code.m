clc; clear; close;
% Lab 1
% Alex Aarhus
% AE 303
% Lab Data Parsing
P_Runs = readtable('Lab 2 Data SP24 Tue.xlsx', 'sheet', 'Group 3');
P_0_Raw = P_Runs{3:37,2:5}; P_2_Raw = P_Runs{3:37,6:9}; P_7_Raw = P_Runs{3:37,10:13};
% Average of each row for each pressure reading
P_0 = mean(P_0_Raw, 2);P_2 = mean(P_2_Raw, 2);P_7 = mean(P_7_Raw, 2);
% Ambient conditions
Pamb = 30.01; % InHg
Tamb = 75.2; % F
% Latitude Correction
SD_lat = 32.7157;
Lat_correction29 = interp1([32 34],[0.035 0.030], SD_lat);
Lat_correction30 = interp1([32 34],[0.036 0.031], SD_lat);
Lat_correction = interp1([29 30], [Lat_correction29 Lat_correction30], Pamb, 'linear', 'extrap');
% Temperature Correction
temps = 72:2:76;
correction29 = [.114, .119, .124];
correction30 = [.118, .123, .128];
correction31 = [.122, .127, .133];
T_correct29 = interp1(temps, correction29, Tamb);
T_correct30 = interp1(temps, correction30, Tamb);
T_correct31 = interp1(temps, correction31, Tamb);
Temp_correction = interp1([29 30 31], [T_correct29 T_correct30 T_correct30], Pamb);
% Total inHg correction
Corrected_h = Pamb - Lat_correction - Temp_correction;
Corrected_h = convlength(Corrected_h,'in','m'); % m
% Corrected Pressure Calculation
Standard_gravity = 9.80665; % m/s^2
rhoHG = 13595; % kg/m^3
Corrected_Pressure_Amb = rhoHG.*Standard_gravity*Corrected_h; % kg/m^3 * m/s^2 * m= Pa
Corrected_Pressure_Amb = convpres(Corrected_Pressure_Amb, "Pa", "psi");
Temps = [76.5667; 76.1489; 86.5359];
% q indicated
q_indicated0 = [P_0(1:6); P_0(8:35)];
q_indicated2 = [P_2(1:6); P_2(8:35)];
q_indicated7 = [P_7(1:6); P_7(8:35)];
% q true
q_true0 = q_indicated0 - P_0(7);
q_true2 = q_indicated2 - P_2(7);
q_true7 = q_indicated7 - P_7(7);
% Tunnel flow uniformity
q_uniform0 = (( q_true0 - mean(q_true0, 'all') ) / mean(q_true0, 'all')) * 100;
q_uniform2 = (( q_true2 - mean(q_true2, 'all') ) / mean(q_true2, 'all')) * 100;
q_uniform7 = (( q_true7 - mean(q_true7, 'all') ) / mean(q_true7, 'all')) * 100;
% Corrected Barometer Pressure for Density
P_corr0 = Corrected_Pressure_Amb + P_0(7) - P_0(7);P_corr0 = convpres(P_corr0, "psi", "Pa");
P_corr2 = Corrected_Pressure_Amb + P_2(7) - P_0(7);P_corr2 = convpres(P_corr2, "psi", "Pa");
P_corr7 = Corrected_Pressure_Amb + P_7(7) - P_0(7);P_corr7 = convpres(P_corr7, "psi", "Pa");
%Density Calcs
R = 287.1;
rho0 = P_corr0 / (R * convtemp(Temps(1), 'F', 'K'));
rho2 = P_corr2 / (R * convtemp(Temps(2), 'F', 'K'));
rho7 = P_corr7 / (R * convtemp(Temps(3), 'F', 'K'));
v0 = (sqrt((2 .* abs(convpres(q_true0,'psi', 'Pa'))) ./ rho0));
v2 = (sqrt((2 .* abs(convpres(q_true2,'psi', 'Pa'))) ./ rho2));
v7 = (sqrt((2 .* abs(convpres(q_true7,'psi', 'Pa'))) ./ rho7));
% Tunnel flow uniformity
v_uniform0 = (( v0 - mean(v0, 'all') ) / mean(v0, 'all')) * 100;
v_uniform2 = (( v2 - mean(v2, 'all') ) / mean(v2, 'all')) * 100;
v_uniform7 = (( v7 - mean(v7, 'all') ) / mean(v7, 'all')) * 100;
%unit Re = rho v / meu
Re0 = (rho0 * v0) / (18.34 * 10^-6);
Re2 = (rho2 * v2) / (18.34 * 10^-6);
Re7 = (rho7 * v7) / (18.6 * 10^-6);
disp('Res')
disp([mean(Re0, 'all'), mean(Re2, 'all'), mean(Re7, 'all')])
%% Visualization
% Measured Pressure Uniformity
fig = figure('Visible', 'off');
plot(P_0, 'o');
hold on;
plot(P_2, 'o');
plot(P_7, 'o');
yline([mean(P_0, 'all'), mean(P_2, 'all'), mean(P_7, 'all')], 'r');
grid on; grid minor;
xlabel('Pressure Port');ylabel('Pressure [psi]');
title('Flow uniformity - Measured Pressure');
legend('0.0 inH_2O', '2.0 inH_2O', '7.0 inH_2O', 'Average', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow uniformity_Measured Pressure.pdf', 'ContentType', 'vector');
close(fig);
% Dynamic Pressure Uniformity
fig = figure('Visible', 'off');
plot(q_true0, 'o');
hold on;
plot(q_true2, 'o');
plot(q_true7, 'o');
yline([0.0, 0.0721824, 0.252638], 'k');
yline([mean(q_true0, 'all'), mean(q_true2, 'all'), mean(q_true7, 'all')], 'r');
grid on; grid minor;
xlabel('Pressure Port');
ylabel('Pressure [psi]');
title('Flow Uniformity - Dynamic Pressure');
legend('0.0 inH_2O', '2.0 inH_2O', '7.0 inH_2O', 'Average', '','','q_{setting}', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow Uniformity_Dynamic Pressure.pdf', 'ContentType', 'vector');
close(fig);
% Airspeed Uniformity
fig = figure('Visible', 'off');
plot(v0, 'o');
hold on;
plot(v2, 'o');
plot(v7, 'o');
yline([mean(v0, 'all'), mean(v2, 'all'), mean(v7, 'all')], 'k');
yline([0, 29, 55], 'r');
grid on; grid minor;
xlabel('Pressure Port');
ylabel('Airspeed [m/s]');
title('Flow Uniformity - Test Section Airspeed');
legend('0.0 inH_2O', '2.0 inH_2O', '7.0 inH_2O', 'Average','','', 'V_{setting}', 'Location', 'northeast');
hold off;
exportgraphics(fig, 'Flow Uniformity_Test Section Airspeed.pdf', 'ContentType', 'vector' );
close(fig);
% Dynamic Pressure
fig = figure('Visible', 'off');
plot(q_uniform2, 'ko')
hold on
plot(q_uniform7, 'rd')
grid on; grid minor;
xlabel('Pressure Port');
ylabel('$\frac{q - \bar{q}}{\bar{q}}$','Interpreter', 'latex');
title('Flow Uniformity - Dynamic Pressure Deviation from Average'); 
%%%%%%%%%%%%%%nH_2O', '7.0 inH_2O', Location='best');
ax = gca;
ax.YAxis.TickLabelFormat = '%g%%';
hold off
exportgraphics(fig, 'Flow Uniformity_Dynamic Pressure Deviation from Average.pdf', 'ContentType', 'vector');
close(fig);
% Airspeed
fig = figure('Visible', 'off');
plot(v_uniform2, 'ko')
hold on
plot(v_uniform7, 'rd')
grid on; grid minor;
xlabel('Pressure Port');
ylabel('$\frac{v - \bar{v}}{\bar{v}}$','Interpreter', 'latex');
title('Flow Uniformity - Airspeed Deviation from Average');
legend('2.0 inH_2O', '7.0 inH_2O', Location='best');
ax = gca;
ax.YAxis.TickLabelFormat = '%g%%';
hold off
exportgraphics(fig, 'Flow Uniformity_Airspeed Deviation from Average.pdf', 'ContentType', 'vector');
close(fig);
%%
x = [-19.5000, -16.5000, -13.5000, -10.5000, -7.5000, -4.5000, 0, 1.5000, 4.5000, 7.5000, 10.5000, 13.5000, 16.5000, 19.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, -7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000, 7.5000]';
z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.5000, 10.5000, 7.5000, 4.5000, 1.5000, -1.5000, -4.5000, -7.5000, -10.5000, -13.5000, 13.5000, 10.5000, 7.5000, 4.5000, 1.5000, -1.5000, -4.5000, -7.5000, -10.5000, -13.5000]';
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
% Contour of Dynamic Pressure Deviation for 7in Plot %
fig = figure('Visible', 'off');
dq = q_uniform7;
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
title('Dynamic Pressure Deviation [q=7.0 inH_2O]');
hold off;
exportgraphics(fig, 'Dynamic Pressure Deviation q=7.0 inH_2O.pdf', 'ContentType', 'vector');
close(fig);
% Contour of Airspeed Distribution for 7in Plot, IAS = 54.4 m/s %
fig = figure('Visible', 'off');
dq = v7;
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
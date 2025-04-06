clc
clearvars
%%
%Parsing
Sph_4_Raw = readtable('Lab 3 Data 2025SP', 'Sheet', '4 in Sphere','Range', 'E7:ADX66');
Sph_4987_Raw = readtable('Lab 3 Data 2025SP', 'Sheet', '4.987 in Sphere', 'Range', 'E7:ADX66');
Sph_6_Raw = readtable('Lab 3 Data 2025SP', 'Sheet', '6 in Sphere','Range', 'E7:ADX66');
%%
temps = readtable('Lab 3 Data 2025SP', 'Sheet', 'TEMPERATURE DATA', 'Range', 'B6:G20');
temps = renamevars(temps, ["Var1","Var2","Var3","Var4","Var5","Var6"],["q4","Temp4","q4987", "Temp4987","q6", "Temp6"]);
% Ambient conditions
Pamb = 29.94; % InHg
Tamb = 75.3; % F
%% -- Correct Barometer Pressure --
Corrected_Pressure_Amb = Correct_Pres(Tamb, Pamb);
%% -- Convert q (H2O) to psi --
q_S = [temps.q4, temps.q4987, temps.q6]; % inH2O
q_S = q_S .* 0.03613; % psi
%% -- Obtain delta P from readings --
deltas_4 = dP(S_4_Raw, 62) - dmP(S_4_Raw, 61); % psi
deltas_4987 = dP(S_4987_Raw, 62) - dP(S_4987_Raw, 61); % psi
deltas_6 = dP(S_6_Raw, 62) - dP(S_6_Raw, 61); % psi
%% -- Compute delta P / q for each setting of q in H2O --
P_q_4 = deltas_4 ./ q_S(:, 1) ;
P_q_4987 = deltas_4987 ./ q_S(:, 2) ;
P_q_6 = deltas_6 ./ q_S(:, 3) ;
%% -- Calculate Re -- Use temp inside the test section --
% Ref NACA Report 1157
meu_0 = 0.000017907216866; % N-s / m^2
T0 = 518.6; T0 = T0 - 459.67; % R -> F
meu(:,1) = meu_0 .* (temps.Temp4 ./ T0) .^ 0.76;
meu(:,2) = meu_0 .* (temps.Temp4987 ./ T0) .^ 0.76;
meu(:,3) = meu_0 .* (temps.Temp6 ./ T0) .^ 0.76;
% Corrected Barometer Pressure for Density
P_0s = dP2(Sph_4_Raw); P_0 = P_0s(1);
P_corr4 = Corrected_Pressure_Amb + P_0s - P_0;
P_corr4 = P_corr4 * 6895;
P_0s = dP2(Sph_4987_Raw); P_0 = P_0s(1);
P_corr4987 = Corrected_Pressure_Amb + P_0s - P_0;
P_corr4987 = P_corr4987 * 6895;
P_0s = dP2(Sph_6_Raw); P_0 = P_0s(1);
P_corr6 = Corrected_Pressure_Amb + P_0s - P_0;
P_corr6 = P_corr6 * 6895;
% Density Calcs
R = 287.1;
rho4 = P_corr4 ./ (R * convtemp(temps.Temp4, 'F', 'K'));
rho4987 = P_corr4987 ./ (R * convtemp(temps.Temp4987, 'F', 'K'));
rho6 = P_corr6 ./ (R * convtemp(temps.Temp6, 'F', 'K'));
v4 = (sqrt((2 .* abs(convpres(q_S(:,1),'psi', 'Pa'))) ./ rho4));
v4987 = (sqrt((2 .* abs(convpres(q_S(:,2),'psi', 'Pa'))) ./ rho4987));
v6 = (sqrt((2 .* abs(convpres(q_S(:,3),'psi', 'Pa'))) ./ rho6));
% Unit Re = rho v / meu
Ds = convlength([4 4.987 6], 'in', 'm');
Re4 = (rho4 .* v4 .* Ds(1)) ./ (meu(:,1));
Re4987 = (rho4987 .* v4987 .* Ds(2)) ./ (meu(:,2));
Re6 = (rho6 .* v6 .* Ds(3)) ./ (meu(:,3));
%% -- Calculate TF --
% Polynomial fitting for different sphere sizes
polyDegree = 5;
fitSphere4 = polyfit(P_q_4(2:end), Re4(2:end), polyDegree);
fitSphere4987 = polyfit(P_q_4987(2:end), Re4987(2:end), polyDegree);
fitSphere6 = polyfit(P_q_6(2:end), Re6(2:end), polyDegree);
% Calculate critical Reynolds numbers
criticalRePoint = 1.22;
criticalRe4 = polyval(fitSphere4, criticalRePoint);
criticalRe4987 = polyval(fitSphere4987, criticalRePoint);
criticalRe6 = polyval(fitSphere6, criticalRePoint);
% Calculate unit Reynolds numbers
unitRe4 = criticalRe4 / Ds(1);
unitRe4987 = criticalRe4987 / Ds(2);
unitRe6 = criticalRe6 / Ds(3);
% Calculate turbulence factors
turbFactor4 = 385000 / criticalRe4;
turbFactor4987 = 385000 / criticalRe4987;
turbFactor6 = 385000 / criticalRe6;
%% -- Visualization for Reynolds vs Pressure Gradient --
fig = figure('Visible', 'off');
scatter(Re4(2:end), P_q_4(2:end), 'kd', 'DisplayName', 'D = 4.000 in');
hold on;
scatter(criticalRe4, criticalRePoint, 'filled', 'kd', 'DisplayName', 'Critical Re D =4.000 in');
scatter(Re4987(2:end), P_q_4987(2:end), 'rhexagram', 'DisplayName', 'D = 4.987 in');
scatter(criticalRe4987, criticalRePoint, 'filled', 'rhexagram', 'DisplayName','Critical Re D = 4.987 in');
scatter(Re6(2:end), P_q_6(2:end), 'bo', 'DisplayName', 'D = 6.000 in');
scatter(criticalRe6, criticalRePoint, 'filled', 'bo', 'DisplayName', 'Critical Re D =6 in');
xlabel('Reynolds Number');
ylabel('$\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
legend('show');
grid on; grid minor
title('Reynolds vs $\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
exportgraphics(fig, 'Reynolds_v_dp_q.pdf');
close(fig);
%% -- Visualization for Reynolds vs Pressure Gradient --
fig = figure;
scatter(Re4(2:end), P_q_4(2:end), 'kd', 'DisplayName', 'D = 4.000 in');
hold on;
scatter(criticalRe4, criticalRePoint, 'filled', 'kd', 'DisplayName', 'Critical Re D =4.000 in');
xlabel('Reynolds Number');
ylabel('$\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
legend('show');
title('Reynolds vs $\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
grid minor; grid on
exportgraphics(fig, 'Reynolds_v_dp_q_4in.pdf');
%%
fig = figure;
scatter(Re4987(2:end), P_q_4987(2:end), 'rhexagram', 'DisplayName', 'D = 4.987 in');
hold on;
scatter(criticalRe4987, criticalRePoint, 'filled', 'rhexagram', 'DisplayName','Critical Re D = 4.987 in');
x = xlabel('Reynolds Number')
ylabel('$\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
legend('show');
title('Reynolds vs $\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
grid minor; grid on
exportgraphics(fig, 'Reynolds_v_dp_q_4987in.pdf');
%%
fig = figure;
scatter(Re6(2:end), P_q_6(2:end), 'bo', 'DisplayName', 'D = 6.000 in');
hold on;
scatter(criticalRe6, criticalRePoint, 'filled', 'bo', 'DisplayName', 'Critical Re D =6 in');
xlabel('Reynolds Number');
ylabel('$\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
legend('show');
title('Reynolds vs $\frac{{\Delta}p}{q}$', 'Interpreter', 'Latex');
grid minor; grid on
exportgraphics(fig, 'Reynolds_v_dp_q_6in.pdf');
%% -- Visualization for Turbulence Factor and Percent Turbulence --
fig = figure('Visible', 'off');
scatter(unitRe4, turbFactor4, 'filled', 'kd', 'DisplayName', 'Turbulence Factor D = 4in');
hold on;
scatter(unitRe4987, turbFactor4987, 'filled', 'rhexagram', 'DisplayName', 'TurbulenceFactor D = 4.987 in');
scatter(unitRe6, turbFactor6, 'filled', 'bo', 'DisplayName', 'Turbulence Factor D = 6in');
xlabel('Critical Unit Reynolds Number');
ylabel('Turbulence Factor');
title('Turbulence Factor vs. Critical Unit Reynolds Number');
legend('show');
ylim([1.6, 1.8])
yyaxis right;
xRange = linspace(min(turbFactor4), max(turbFactor6), 100);
yLimitsRightAxis = [0.7, 1.0];
ylim(yLimitsRightAxis);
ylabel('Percent Turbulence');
grid on; grid minor
exportgraphics(fig, 'TF_v_Re.pdf');
close(fig);
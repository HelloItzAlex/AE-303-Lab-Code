function f_interp = bilinearInterpolation(x1, x2, y1, y2, fQ11, fQ21, fQ12, fQ22, x, y)
 % Ensure x1 < x2 and y1 < y2
 if x1 >= x2 || y1 >= y2
 error('x1 < x2 and y1 < y2.');
 end
 % Compute the weights
 f_interp = (fQ11 * (x2 - x) * (y2 - y) + fQ21 * (x - x1) * (y2 - y) + fQ12 * (x2 - x) * (y - y1) + fQ22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1));
end
%% Bilinear Interpolation - Temperature Correction inHg
% Define the coordinates of the known points
T1 = 72; 
T2 = 74;
B1 = 29; 
B2 = 30;
% Defining the Known Values
Q11 = 0.119; Q12 = 0.123;
Q21 = 0.114; Q22 = 0.118;
% Defining the points to interpolate
T = 73.66769; B = 29.81385;
% Performing the Interpolation
interpolatedValue = bilinearInterpolation(T1, T2, B1, B2, Q11, Q12, Q21, Q22, T, B);
disp(['Interpolated Value:', num2str(interpolatedValue)])
%% Bilinear Interpolation - Gravity Correction inHg
% Define the coordinates of the known points
H1 = 32; H2 = 34;
B1 = 29; B2 = 30;
% Defining the Known Values
Q11 = 0.035; Q12 = 0.036;
Q21 = 0.03; Q22 = 0.031;
% Defining the points to interpolate
H = 32.71574; B = 29.81385;
% Performing the Interpolation
interpolatedValue = bilinearInterpolation(H1, H2, B1, B2, Q11, Q12, Q21, Q22, H, B);
disp(['Interpolated Value:', num2str(interpolatedValue)])
%% Bilinear Interpolation - Temperature Correction Pa
% Define the coordinates of the known points
T1 = 72; T2 = 74;
B1 = 98204.98; B2 = 101591.36;
% Defining the Known Values
Q11 = 402.98; Q12 = 416.52;
Q21 = 386.05; Q22 = 399.59;
% Defining the points to interpolate
T = 73.66769; B = 100960.99;
% Performing the Interpolation
interpolatedValue = bilinearInterpolation(T1, T2, B1, B2, Q11, Q12, Q21, Q22, T, B);
disp(['Interpolated Value:', num2str(interpolatedValue)])
%% Bilinear Interpolation - Gravity Correction Pa
% Define the coordinates of the known points
H1 = 32; H2 = 34;
B1 = 98204.98; B2 = 101591.36;
% Defining the Known Values
Q11 = 118.52; Q12 = 121.91;
Q21 = 101.59; Q22 = 104.98;
% Defining the points to interpolate
H = 32.71574; B = 100960.99;
% Performing the Interpolation
interpolatedValue = bilinearInterpolation(H1, H2, B1, B2, Q11, Q12, Q21, Q22, H, B);
disp(['Interpolated Value:', num2str(interpolatedValue)])
%% Data Processing Code
% Raw Data
rawinHg1 = [29.90; 29.90; 29.91; 29.97; 29.99; 29.97; 29.88; 29.90; 29.88; 29.97;29.98; 29.99; 29.99; 29.99; 29.99; 29.99; 29.84; 29.90; 29.97; 29.99; 29.98;29.99];
rawinHg2 = [29.88; 29.88; 29.85; 29.89; 29.89; 29.97; 29.88; 29.99; 29.86; 29.87;29.85; 29.81; 29.86; 29.99; 29.87; 29.85; 29.87; 29.99; 29.79; 29.84; 29.73];
rawinHg3 = [29.88; 29.99; 29.98; 29.04; 28.88; 28.88; 28.88; 29.78; 28.78; 29.78;29.84; 29.97; 29.88; 30.05; 29.90; 29.88; 29.88; 29.89; 29.88; 29.86; 28.84;29.88];

rawmmHg1 = [759.8; 759.3; 759.2;759.1;759.1;759.3;759.8;760.0;758.8;759.3;759.4;758.7;759.2;759.2;759.2;759.1;759.1;759.3;759.6;759.6;759.8;759.1];
rawmmHg2 = [756.8; 756.2; 758.0; 759.9; 758.9; 758.6; 759.9; 759.3; 756.9; 758.6;757.2; 756.3; 756.7; 759.2; 758.9; 759.4; 759.7; 759.9; 756.1; 757.8; 757.8];
rawmmHg3 = [758.8;759.2;758.6;758.2;758.8;759.0;758.6;759.8;758.6;759.6;758.0;758.6;758.6;760.0;759.8;758.6;758.8;759.1;759.0;758.6;758.9;758.7];
T1 = [74;74;73;73;73;74;73.5;73.5;74;74;73;73.5;77.3;73.5;74;73;74;73.6;73;73;74;73];
T2 = [73;74;74.1;73.6;74;74;74.8;73;73;74;73.8;73.5;73.8;74;73.5;73.6;74;74;73.7;73.8;74.1];
T3 = [72.6;74;73;74.1;74;74;74.5;74.6;74;73;73.8;74;73.5;74;74;73;73;73;73.1;73;73;73];
EnsembleinHg = [29.90; 29.90; 29.91; 29.97; 29.99; 29.97; 29.88; 29.90; 29.88;
29.97; 29.98; 29.99; 29.99; 29.99; 29.99; 29.99; 29.84; 29.90; 29.97; 29.99;
29.98; 29.99; 29.88; 29.88; 29.85; 29.89; 29.89; 29.97; 29.88; 29.99; 29.86;
29.87; 29.85; 29.81; 29.86; 29.99; 29.87; 29.85; 29.87; 29.99; 29.79; 29.84;
29.73; 29.88; 29.99; 29.98; 29.04; 28.88; 28.88; 28.88; 29.78; 28.78; 29.78;
29.84; 29.97; 29.88; 30.05; 29.90; 29.88; 29.88; 29.89; 29.88; 29.86; 28.84;
29.88];
EnsemblemmHg = [759.8; 759.3; 759.2;
759.1;759.1;759.3;759.8;760.0;758.8;759.3;759.4;758.7;759.2;759.2;759.2;759.1;759.
1;759.3;759.6;759.6;759.8;759.1;756.8; 756.2; 758.0; 759.9; 758.9; 758.6; 759.9;
759.3; 756.9; 758.6; 757.2; 756.3; 756.7; 759.2; 758.9; 759.4; 759.7; 759.9;
756.1; 757.8; 757.8;
758.8;759.2;758.6;758.2;758.8;759.0;758.6;759.8;758.6;759.6;758.0;758.6;758.6;760.
0;759.8;758.6;758.8;759.1;759.0;758.6;758.9;758.7];
EnsembleTF = [74;74;73;73;73;74;73.5;73.5;74;74;73;73.5;77.3;73.5;74;73;74;73.6;73;73;74;73;73;
74;74.1;73.6;74;74;74.8;73;73;74;73.8;73.5;73.8;74;73.5;73.6;74;74;73.7;73.8;74.1;
72.6;74;73;74.1;74;74;74.5;74.6;74;73;73.8;74;73.5;74;74;73;73;73;73.1;73;73;73];
% Constants and Known Values
N1 = numel(rawinHg1);
N2 = numel(rawinHg2);
N3 = numel(rawinHg3);
EnsembleN = N1+N2+N3;
R = 287.058;
rhoHg = 13595;
g = 9.80665;
t_vp = 2.08;
t_vpEnsemble = 2.0;
TempCorrectioninHg = 0.11827;
GravityCorrectioninHg = 0.031289;
TempCorrectionPa = 400.4918;
GravityCorrectionPa = 105.9547;

% Converting mmHg to Pa
rawPa1 = rawmmHg1*133.32239;
rawPa2 = rawmmHg2*133.32239;
rawPa3 = rawmmHg3*133.32239;
% Finding the Mean of Raw Data
rawinHgAvg1 = sum(rawinHg1)/N1;
rawinHgAvg2 = sum(rawinHg2)/N2;
rawinHgAvg3 = sum(rawinHg3)/N3;
rawPaAvg1 = sum(rawPa1)/N1;
rawPaAvg2 = sum(rawPa2)/N2;
rawPaAvg3 = sum(rawPa3)/N3;
TAvg1 = sum(T1)/N1;
TAvg2 = sum(T2)/N2;
TAvg3 = sum(T3)/N3;
% Finding Standard Deviation of Raw Data
rawS_xinHg1 = sqrt((sum((rawinHg1 - rawinHgAvg1).^2)/(N1-1)));
rawS_xinHg2 = sqrt((sum((rawinHg2 - rawinHgAvg2).^2)/(N2-1)));
rawS_xinHg3 = sqrt((sum((rawinHg3 - rawinHgAvg3).^2)/(N3-1)));
rawS_xPa1 = sqrt((sum((rawPa1 - rawPaAvg1).^2)/(N1-1)));
rawS_xPa2 = sqrt((sum((rawPa2 - rawPaAvg2).^2)/(N2-1)));
rawS_xPa3 = sqrt((sum((rawPa3 - rawPaAvg3).^2)/(N3-1)));
S_xT1 = sqrt((sum((T1 - TAvg1).^2)/(N1-1)));
S_xT2 = sqrt((sum((T2 - TAvg2).^2)/(N2-1)));
S_xT3 = sqrt((sum((T3 - TAvg3).^2)/(N3-1)));
% Correcting Barometer (inHg) and (Pa) Data
CorrectedinHg1 = rawinHg1 - GravityCorrectioninHg - TempCorrectioninHg;
CorrectedinHg2 = rawinHg2 - GravityCorrectioninHg - TempCorrectioninHg;
CorrectedinHg3 = rawinHg3 - GravityCorrectioninHg - TempCorrectioninHg;
CorrectedinHgAvg1 = sum(CorrectedinHg1)/N1;
CorrectedinHgAvg2 = sum(CorrectedinHg2)/N2;
CorrectedinHgAvg3 = sum(CorrectedinHg3)/N3;
CorrectedPa1 = rawPa1 - GravityCorrectionPa - TempCorrectionPa;
CorrectedPa2 = rawPa2 - GravityCorrectionPa - TempCorrectionPa;
CorrectedPa3 = rawPa3 - GravityCorrectionPa - TempCorrectionPa;
CorrectedPaAvg1 = sum(CorrectedPa1)/N1;
CorrectedPaAvg2 = sum(CorrectedPa2)/N2;
CorrectedPaAvg3 = sum(CorrectedPa3)/N3;
% Finding Standard Deviation of Corrected Data
S_xinHg1 = sqrt((sum((CorrectedinHg1 - CorrectedinHgAvg1).^2)/(N1-1)));
S_xinHg2 = sqrt((sum((CorrectedinHg2 - CorrectedinHgAvg2).^2)/(N2-1)));
S_xinHg3 = sqrt((sum((CorrectedinHg3 - CorrectedinHgAvg3).^2)/(N3-1)));
S_xPa1 = sqrt((sum((CorrectedPa1 - CorrectedPaAvg1).^2)/(N1-1)));

S_xPa2 = sqrt((sum((CorrectedPa2 - CorrectedPaAvg2).^2)/(N2-1)));
S_xPa3 = sqrt((sum((CorrectedPa3 - CorrectedPaAvg3).^2)/(N3-1)));
% Find the sample standard deviation of the means
S_xbarinHg1 = S_xinHg1/sqrt(N1);
S_xbarinHg2 = S_xinHg2/sqrt(N2);
S_xbarinHg3 = S_xinHg3/sqrt(N3);
S_xbarPa1 = S_xPa1/sqrt(N1);
S_xbarPa2 = S_xPa2/sqrt(N2);
S_xbarPa3 = S_xPa3/sqrt(N3);
% Calculating the Ensemble Average
rawEnsembleAvginHg = (sum(EnsembleinHg)/(EnsembleN));
rawEnsembleStdinHg = sqrt((sum((EnsembleinHg - rawEnsembleAvginHg).^2))/(EnsembleN));
CorrectedEnsembleinHg = EnsembleinHg - GravityCorrectioninHg - TempCorrectioninHg;
CorrectedEnsembleMeaninHg = (sum(CorrectedEnsembleinHg)/(EnsembleN));
EnsembleStdinHg = sqrt((sum((CorrectedEnsembleinHg - CorrectedEnsembleMeaninHg).^2))/(EnsembleN));
RawEnsemblePa = EnsemblemmHg*133.32239;
CorrectedEnsemblePa = RawEnsemblePa - GravityCorrectionPa - TempCorrectionPa;
CorrectedEnsembleMeanPa = (sum(CorrectedEnsemblePa)/(EnsembleN));
EnsembleStdPa = sqrt((sum((CorrectedEnsemblePa - CorrectedEnsembleMeanPa).^2))/(EnsembleN));
% Standard Deviation of Corrected Ensemble Average
StdEnsembleAvginHg = CorrectedEnsembleMeaninHg/sqrt(EnsembleN);
StdEnsembleAvgPa = CorrectedEnsembleMeanPa/sqrt(EnsembleN);
EnsembleT = ((EnsembleTF-32)/1.8)+273.15;
EnsembleAvgT = sum(EnsembleT)/EnsembleN;
StdEnsembleT = sqrt((sum((EnsembleT - EnsembleAvgT).^2)/(EnsembleN-1)));
% Finding Density
rhokgm_inHg1 = CorrectedinHg1./(R.*T1)*3386.38;
rhokgm_inHg2 = CorrectedinHg2./(R.*T2)*3386.38;
rhokgm_inHg3 = CorrectedinHg3./(R.*T3)*3386.38;
rhokgm_Pa1 = CorrectedPa1./(R.*T1);
rhokgm_Pa2 = CorrectedPa2./(R.*T2);
rhokgm_Pa3 = CorrectedPa3./(R.*T3);
rhoensemble_kgm_inHg = CorrectedEnsembleinHg./(R.*EnsembleT)*3386.38;
rhoensemble_kgm_Pa = CorrectedEnsemblePa./(R.*EnsembleT);
% Finding Density Averages
rhoAvgkgm_inHg1 = sum(rhokgm_inHg1)/N1;
rhoAvgkgm_inHg2 = sum(rhokgm_inHg2)/N2;
rhoAvgkgm_inHg3 = sum(rhokgm_inHg3)/N3;

rhoAvgkgm_Pa1 = sum(rhokgm_Pa1)/N1;
rhoAvgkgm_Pa2 = sum(rhokgm_Pa2)/N2;
rhoAvgkgm_Pa3 = sum(rhokgm_Pa3)/N3;
rhoensembleAvg_kgm_inHg = sum(rhoensemble_kgm_inHg)/EnsembleN;
rhoensembleAvg_kgm_Pa = sum(rhoensemble_kgm_Pa)/EnsembleN;
% Table 1 Uncertainty for One-Time Sample
Uncertainty_T1 = t_vp*S_xT1;
Uncertainty_T2 = t_vp*S_xT2;
Uncertainty_T3 = t_vp*S_xT3;
Uncertainty_rawinHg1 = t_vp*rawS_xinHg1;
Uncertainty_rawinHg2 = t_vp*rawS_xinHg2;
Uncertainty_rawinHg3 = t_vp*rawS_xinHg3;
Uncertainty_CorrectedPa1 = t_vp*S_xPa1;
Uncertainty_CorrectedPa2 = t_vp*S_xPa2;
Uncertainty_CorrectedPa3 = t_vp*S_xPa3;
% Table 2 Uncertainty of Ensemble Data
Uncertainty_rawEnsembleinHg = t_vpEnsemble*rawEnsembleStdinHg;
% Error Propagation of One-Time Values
ErrorProp_OneTimepinHg1 = t_vp*S_xinHg1;
ErrorProp_OneTimepinHg2 = t_vp*S_xinHg2;
ErrorProp_OneTimepinHg3 = t_vp*S_xinHg3;
ErrorProp_OneTimeT1 = t_vp*S_xT1;
ErrorProp_OneTimeT2 = t_vp*S_xT2;
ErrorProp_OneTimeT3 = t_vp*S_xT3;
ErrorProp_OneTimeRhoinHg1 = sqrt(((rhoAvgkgm_inHg1/CorrectedinHgAvg1)^2*(ErrorProp_OneTimepinHg1^2))+(rhoAvgkgm_inHg1/TAvg1)^2*(ErrorProp_OneTimeT1^2));
ErrorProp_OneTimeRhoinHg2 = sqrt(((rhoAvgkgm_inHg2/CorrectedinHgAvg2)^2*(ErrorProp_OneTimepinHg2^2))+(rhoAvgkgm_inHg2/TAvg2)^2*(ErrorProp_OneTimeT2^2));
ErrorProp_OneTimeRhoinHg3 = sqrt(((rhoAvgkgm_inHg3/CorrectedinHgAvg3)^2*(ErrorProp_OneTimepinHg3^2))+(rhoAvgkgm_inHg3/TAvg3)^2*(ErrorProp_OneTimeT3^2));
ErrorProp_OneTimepPa1 = t_vp*S_xPa1;
ErrorProp_OneTimepPa2 = t_vp*S_xPa2;
ErrorProp_OneTimepPa3 = t_vp*S_xPa3;
ErrorProp_OneTimeRhoPa1 = sqrt(((rhoAvgkgm_Pa1/CorrectedPaAvg1)^2*(ErrorProp_OneTimepPa1^2))+(rhoAvgkgm_Pa1/TAvg1)^2*(ErrorProp_OneTimeT1^2));
ErrorProp_OneTimeRhoPa2 = sqrt(((rhoAvgkgm_Pa2/CorrectedPaAvg2)^2*(ErrorProp_OneTimepPa2^2))+(rhoAvgkgm_Pa2/TAvg2)^2*(ErrorProp_OneTimeT2^2));
ErrorProp_OneTimeRhoPa3 = sqrt(((rhoAvgkgm_Pa3/CorrectedPaAvg3)^2*(ErrorProp_OneTimepPa3^2))+(rhoAvgkgm_Pa3/TAvg3)^2*(ErrorProp_OneTimeT3^2));
% Error Propagation / Uncertainty of the Ensemble
ErrorProp_EnsembleAvgpinHg = t_vpEnsemble*StdEnsembleAvginHg;
ErrorProp_EnsembleAvgT = t_vpEnsemble*StdEnsembleT;
ErrorProp_EnsembleAvgRhoinHg = sqrt(((rhoensembleAvg_kgm_inHg/CorrectedEnsembleMeaninHg)^2*(ErrorProp_EnsembleAvgpinHg^2))+(rhoensembleAvg_kgm_inHg/EnsembleAvgT)^2*(ErrorProp_EnsembleAvgT^2));
ErrorProp_EnsembleAvgpPa = t_vpEnsemble*StdEnsembleAvgPa;
ErrorProp_EnsembleAvgRhoPa = sqrt(((rhoensembleAvg_kgm_Pa/CorrectedEnsembleMeanPa)^2*(ErrorProp_EnsembleAvgpPa^2))+(rhoensembleAvg_kgm_Pa/EnsembleAvgT)^2*(ErrorProp_EnsembleAvgT^2));
% Resolution of Measurement
resolutionT = 0.1;
resolutioninHg = 0.01;
resolutionPa = 10;
partialDevrho_p = 1/(R*EnsembleAvgT);
partialdevrho_T = -CorrectedEnsembleMeanPa/(R*(EnsembleAvgT^2));
resolutionrho = sqrt(((partialDevrho_p^2)*resolutionPa)+((partialdevrho_T^2)*resolutionT));
% Convergence Test + Plot
n_points = (1:length(CorrectedEnsembleinHg))';
runningMeanCorrectedinHg = arrayfun(@(n) mean(CorrectedEnsembleinHg(1:n)),n_points);
RunningStdDevinHg = arrayfun(@(n) std(CorrectedEnsembleinHg(1:n),1), n_points);
RunningStd_of_MeaninHg = RunningStdDevinHg./sqrt(n_points);
MeanConvergence = runningMeanCorrectedinHg(end)
StdDevConvergence = RunningStdDevinHg(end)
StdofMeanConvergence = RunningStd_of_MeaninHg(end)
figure (1)
subplot(3,1,1)
hold on;
plot(n_points, runningMeanCorrectedinHg, 'LineWidth', 1.5, 'MarkerSize', 6,'color', "#0072BD");
yline(runningMeanCorrectedinHg(end), 'r--', 'LineWidth', 1.5);
title('Convergence of Mean')
subtitle(['Measured Pressure [inHg] --- Convergance Towards $\bar{x}$' ' =29.949'],'Interpreter','Latex')
ylabel('$\bar{x}$ [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:EnsembleN])
grid on
subplot(3,1,2)
plot(n_points, RunningStdDevinHg, 'LineWidth', 1.5, 'MarkerSize', 6, 'color',"#A2142F");
yline(RunningStdDevinHg(end), 'r--', 'LineWidth', 1.5);
title('Convergence of Standard Deviation')
subtitle('Measured Pressure [inHg] --- Convergance Towards Sx = 0.3047','Interpreter','Latex')

ylabel('Sx [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:EnsembleN])
grid on
subplot(3,1,3)
plot(n_points, RunningStd_of_MeaninHg, 'LineWidth', 1.5, 'MarkerSize', 6, 'color', "#77AC30");
yline(RunningStd_of_MeaninHg(end), 'r--', 'LineWidth', 1.5);
title('Convergence of Standard Deviation of the Means')
subtitle(['Measured Pressure [inHg] --- Convergance ''Towards S$\bar{x}$=0.0378'], 'Interpreter','Latex')
ylabel('S$\bar{x}$ [inHg]', 'Interpreter','Latex')
xlabel('Number of Samples', 'Interpreter','Latex')
xticks([0:5:EnsembleN])
grid on
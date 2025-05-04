clc; 
clearvars;
close all; 
clear all;
%%
sheetNames = {'q=0 AoA=0', 'q=5 AoA=-5', 'q=5 AoA=0', 'q=5 AoA=5', 'q=5 AoA=10', 'q=5 AoA=15', 'q=5 AoA=20 (1)', 'q=5 AoA=20 (2)'}; 
sheetNamesCleaned = matlab.lang.makeValidName(sheetNames); 

for i = 1:length(sheetNames)
    

    fullData.(sheetNamesCleaned{i}).Samples = mean(table2array(readtable('Lab 4 Data SP2025 Thursday.xlsx', 'Sheet', sheetNames{i}, 'Range', 'B11:EEW62'))* 6894.76, 2);
    fullData.(sheetNamesCleaned{i}).Static_Pressure = mean(table2array(readtable('Lab 4 Data SP2025 Thursday.xlsx', 'Sheet', sheetNames{i}, 'Range', 'B9:EEW9'))* 6894.76, 2);
    fullData.(sheetNamesCleaned{i}).Total_Pressure = mean(table2array(readtable('Lab 4 Data SP2025 Thursday.xlsx', 'Sheet', sheetNames{i}, 'Range', 'B10:EEW10'))* 6894.76, 2);
    fullData.(sheetNamesCleaned{i}).Temperature = (table2array(readtable('Lab 4 Data SP2025 Thursday.xlsx', 'Sheet', sheetNames{i}, 'Range', 'E7:E7')) + 459.67) * (5/9); 
end

% Ambient conditions
pAmb = 29.93; % InHg
tAmb = 77.8; % F
%% Load Xfoil Data
fileNames = {'01-NACA43012A_a-05.cp', '02-NACA43012A_a+00.cp','03-NACA43012A_a+05.cp', '04-NACA43012A_a+10.cp','05-NACA43012A_a+15.cp', '06-NACA43012A_a+20.cp', '07-NACA43012A_a+20.cp'};
xFoilFieldNames = {'X_a_N5', 'X_a_0', 'X_a_5', 'X_a_10', 'X_a_15', 'X_a_20_1', 'X_a_20_2'};

for i = 1:length(fileNames)
    [dataStruct] = read_txt_File(fileNames{i});
    Xfoil.(xFoilFieldNames{i}) = dataStruct;
end

xCoordsUpper = [0; 0.013; 0.025; 0.048; 0.073; 0.097; 0.15; 0.2; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
xCoordsLower = [0.015; 0.029; 0.055; 0.08; 0.105; 0.157; 0.207; 0.257; 0.306; 0.407; 0.507; 0.608; 0.708; 0.812; 0.912; 1.0]; 
yCoordsUpper = [0; 0.0394; 0.0514; 0.0678; 0.0794; 0.0868; 0.0933; 0.0927; 0.0895; 0.0857; 0.077; 0.0656; 0.052; 0.037; 0.0244; 0.0125];
yCoordsLower = -1.*[0.0089; 0.01156; 0.016; 0.019; 0.0216; 0.0262; 0.0299; 0.0330; 0.0353; 0.0393; 0.0402; 0.0389; 0.0353; 0.0259; 0.0132; 0];

% Manometer 
% Ports 1-16: Airfoil lower surface                                        | 
% Ports 17-32: Airfoil upper surface
% Ports 33-50: Wake rake
%% Water Manometer
load("Manometer.mat")
W = W.*248.84; % inH2O -> Pa
%SCV
% Ports 1-16: Airfoil lower surface
% Ports 17-32: Airfoil upper surface
% Ports 33-52: Wake rake

sheetNames = sheetNames(2:end);
for i = 1:length(sheetNames)

    pS = fullData.(sheetNamesCleaned{i}).Samples;
    q = fullData.(sheetNamesCleaned{i}).Total_Pressure - fullData.(sheetNamesCleaned{i}).Static_Pressure;
    pInf = pS(33) - q;
    cpi.(sheetNamesCleaned{i}) = (pS(1:32) - pInf) ./  q;
   
    Pinf_M = 5*248.84;
    cpi_M.(sheetNamesCleaned{i}) = ((W(26, i) - W(1:32, i)) - Pinf_M) ./  q;
end


titles = {'Config 1, AoA -5^{\circ}, Rake pos 0.00in, Rake Angle 17^{\circ}';'Config 2, AoA 0^{\circ}, Rake pos 0.25in, Rake Angle 6^{\circ}';'Config 3, AoA 5^{\circ}, Rake pos 1.0in, Rake Angle 6^{\circ}';'Config 4, AoA 10^{\circ}, Rake pos 2.25in, Rake Angle 10^{\circ}';
'Config 5, AoA 15^{\circ}, Rake pos 3.00in, Rake Angle 16^{\circ}';'Config 6, AoA 20^{\circ}, Rake pos 1.00in, Rake Angle 90^{\circ}';'Config 7, AoA 20^{\circ}, Rake pos 3.00in, Rake Angle 90^{\circ}'};

figure; 
lowerLabels = arrayfun(@num2str, 1:16, 'UniformOutput', false);upperLabels = arrayfun(@num2str, 17:32, 'UniformOutput', false);
% Main plotting loop
for i = 1:length(xFoilFieldNames)
    f = figure;
    plot(Xfoil.(xFoilFieldNames{i}).x, Xfoil.(xFoilFieldNames{i}).Cp, 'k-');
    hold on;
    plot(xCoordsLower, cpi.(sheetNamesCleaned{i})(1:16), '-o');
    plot(xCoordsUpper(1:16), cpi.(sheetNamesCleaned{i})(17:32), '-d');
    annotatePlot(xCoordsLower, cpi.(sheetNamesCleaned{i})(1:16),lowerLabels, 'bottom', 'right');
    annotatePlot(xCoordsUpper(1:16), cpi.(sheetNamesCleaned{i})(17:32), upperLabels, 'bottom', 'right');
    title(titles{i});
    legend('Theoretical', 'SCV Lower', 'SCV Upper'); 
    set(gca, 'YDir', 'reverse');
    grid on, grid minor;
    pdfFileName = sprintf('Cp_plot_%d.png', i); 
    print(f, pdfFileName, '-dpng'); 
end

c = 0.3048; 
AoA = [-5:5:20, 20]; 
%% Finding Cn Ca Cm
dYdXUpper = ( yCoordsUpper(16)-yCoordsUpper(1) ) / (xCoordsUpper(16)-xCoordsUpper(1)); 
dYdXLower = (yCoordsLower(16)-yCoordsLower(1) ) / (xCoordsLower(16)-xCoordsLower(1)); 
CnSCV = zeros(1, 7);CaSCV = zeros(1, 7);CmSCV = zeros(1, 7);Cl_SCV = zeros(1, 7);Cd_SCV = zeros(1, 7);

for  i = 1:length(sheetNames)
    CpLower =  cpi.(sheetNamesCleaned{i})(1:16); Cp_Upper = cpi.(sheetNamesCleaned{i})(17:32);
    CnSCV(i) = trapz( xCoordsLower , CpLower ) - trapz( xCoordsUpper , Cp_Upper );
    CaSCV(i) = trapz( xCoordsUpper , Cp_Upper .* dYdXUpper ) - trapz( xCoordsLower , CpLower .* dYdXLower );
    CmSCV(i) = trapz( xCoordsUpper , Cp_Upper .* xCoordsUpper ) - trapz( xCoordsLower , CpLower .* xCoordsLower )...
               +trapz( xCoordsUpper , Cp_Upper .* dYdXUpper .* yCoordsUpper ) - trapz( xCoordsLower , CpLower.* dYdXLower .* yCoordsLower );
    Cl_SCV(i) = CnSCV(i) * cosd( AoA(i) ) - CaSCV(i) * sind( AoA(i) );  
    Cd_SCV(i) = CnSCV(i) * sind( AoA(i) ) + CaSCV(i) * cosd( AoA(i) ); 
end 
xAc = 0.238; y_ac = 0.070;
Cm_ac_SCV = CmSCV + ( Cl_SCV .* xAc .* cosd(AoA) ) - ( Cd_SCV .* y_ac .* cosd( AoA ) ) + ( Cl_SCV .* y_ac .* sind( AoA ) ) + ( Cd_SCV .* xAc .* sind(AoA) );

%% Wake rake 
centerDistances = convlength([-6.75; -4.75; -3.75; -3.25; -2.75; -2.25; -1.75; -1.25;-0.75; -0.25; .25; .75; 1.25; 1.75; 2.25; 2.75; 3.25; 3.75; 4.75;6.75], 'in', 'm');
wakeAngles = [17 6 6 10 16 90 90];
CdWakeSCV = zeros(7, 1); 
for i = 1:length(sheetNames)
    q1 = fullData.(sheetNamesCleaned{i}).Total_Pressure - fullData.(sheetNamesCleaned{i}).Static_Pressure; 
    q2 = fullData.q_0AoA_0.Samples(33:52) ;
    
    wakeData = sqrt( (q2 ./ q1) ) - (q2 ./ q1);
    CdWakeSCV(i) = (2 / 0.03048) * trapz( centerDistances .* sind(wakeAngles(i)), wakeData);
end

% Calculation of wake velocity
R = 287.1; 

for i = 1:length(sheetNames)
      Temp = fullData.(sheetNamesCleaned{i}).Temperature; % K 
      P = pAmb + fullData.(sheetNamesCleaned{i}).Static_Pressure - fullData.q_0AoA_0.Static_Pressure; % Pa
      rho  = pAmb / (R * Temp);
      v_SCV(i, :) = sqrt( (2 .*  abs(fullData.(sheetNamesCleaned{i}).Samples(33:52))) ./ rho ); % m/s
end

%% Visualization 
load('NACA_43012A.mat')
f = figure;
    plot(NACA_43012A(:, 1), NACA_43012A(:, 2))
    axis equal; grid on; grid minor
    pdfFileName = sprintf('NACA_43012A.png'); 
    print(f, pdfFileName, '-dpng'); 
    close(f); 

% Calibration plot
f = figure; 
    plot (fullData.q_0AoA_0.Samples, 'o')
    ylim([-100 100])
    xlim([1 52])
    grid on; grid minor
    AoA_Plot = AoA + [zeros(1, 6) 5];
    pdfFileName = sprintf('Calibration_Plot.png'); 
    print(f, pdfFileName, '-dpng'); 

    %% cl vs alpha
load('AeroCoefficients.mat')
f = figure; 
    plot(AoA_Plot, Cl_SCV, '--o', 'LineWidth', 2, 'DisplayName', 'SCV'); hold on 
    plot(AoA_Plot, NACA_LiftCoefficients, '--s', 'LineWidth', 2, 'DisplayName', 'NACA'); hold on 
    plot(AoA_Plot, Theory_LiftCoefficients, '--^', 'LineWidth', 2, 'DisplayName', 'Theoretical'); hold on 
    xlabel('Angle of Attack (degrees)');ylabel('Coefficient Values');
    title('c_{l} vs. Angle of Attack');
    legend('show');grid on; grid minor; 
    xticks(-5:5:25);xticklabels({'-5', '0', '5', '10', '15', '20 (1)', '20 (2)'});
    pdfFileName = sprintf('Cl_v_Alpha.png'); 
    print(f, pdfFileName, '-dpng');

%% cl vs cd
f = figure; 
    plot(Cl_SCV, Cd_SCV, '--o', 'LineWidth', 2, 'DisplayName', 'SCV'); hold on 
    plot(NACA_LiftCoefficients, NACA_DragCoefficients, '--s', 'LineWidth', 2, 'DisplayName', 'NACA'); hold on 
    plot(Theory_LiftCoefficients, Theory_DragCoefficients, '--^', 'LineWidth', 2, 'DisplayName', 'Theoretical'); hold on 
    plot(Cl_SCV, CdWakeSCV, '--^', 'LineWidth', 2, 'DisplayName', 'Wake Rake'); hold on 
    xlabel('c_{l}');ylabel('c_{d}');
    title('c_{l} vs. c_{d}');
    legend('show');grid on; grid minor; 
    pdfFileName = sprintf('Cl_v_Cd.png'); 
    print(f, pdfFileName, '-dpng'); 

%% cl vs cm_ac
f = figure; 
    plot(Cl_SCV, Cm_ac_SCV, '--o', 'LineWidth', 2, 'DisplayName', 'SCV'); hold on 
    plot(NACA_LiftCoefficients, NACA_MomentCoefficients, '--s', 'LineWidth', 2, 'DisplayName', 'NACA'); hold on 
    plot(Theory_LiftCoefficients, Theory_DragCoefficients, '--^', 'LineWidth', 2, 'DisplayName', 'Theoretical'); hold on 
    xlabel('c_{l}');ylabel('c_{m_ac}'); 
    title('c_{l} vs. c_{m_ac}');
    legend('show');grid on; grid minor; 
    pdfFileName = sprintf('Cl_v_Cm_ac.png'); 
    print(f, pdfFileName, '-dpng');

%% Wake Velocity Profile
figure; 
titles = {'Config 1, AoA -5^{\circ}, Rake pos 0.00in, Rake Angle 17^{\circ}';'Config 2, AoA 0^{\circ}, Rake pos 0.25in, Rake Angle 6^{\circ}';'Config 3, AoA 5^{\circ}, Rake pos 1.0in, Rake Angle 6^{\circ}';'Config 4, AoA 10^{\circ}, Rake pos 2.25in, Rake Angle 10^{\circ}';
'Config 5, AoA 15^{\circ}, Rake pos 3.00in, Rake Angle 16^{\circ}';'Config 6, AoA 20^{\circ}, Rake pos 1.00in, Rake Angle 90^{\circ}';'Config 7, AoA 20^{\circ}, Rake pos 3.00in, Rake Angle 90^{\circ}'};

for i = 1:length(sheetNames)
    f = figure; 
    plot(centerDistances, v_SCV(i, :), '--.', 'LineWidth', 2, 'DisplayName', 'SCV'); hold on 
    ylabel('Local Wind Speed (m/s)');xlabel('Distance from Center Line');
    title(titles{i});
    legend('show');grid on; grid minor; 

    pdfFileName = sprintf('Wake_Velocity_%d.png', i); 
    print(f, pdfFileName, '-dpng'); 
end

close all
---------------------------------------------------------------------------
% Functions

function [dataStruct] = read_txt_File(filePath)
    % Open the file for reading
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Failed to open the file.');
    end
  
    header = {};
    data = [];
    
    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        if ~isempty(line) && ~isletter(line(1))
            numData = sscanf(line, '%f %f %f');
            data = [data; numData'];
        else
            header{end+1} = line;
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Convert the data array into a struct with fields x, y, Cp
    dataStruct = struct('x', data(:,1), 'y', data(:,2), 'Cp', data(:,3));
end

function annotatePlot(X, Y, labels, va, ha)
    for j = 1:length(Y)
        text(X(j), Y(j), labels{j}, 'VerticalAlignment', va, 'HorizontalAlignment', ha);
    end
end

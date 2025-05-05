clc; 
clearvars;
close all; 
clear all;
%%

dataRaw = readtable('Lab 5 Thursday SP25 data Sheet.xlsx', 'Sheet','Ensemble Average FORCE data');
fieldNamesdata = {'a', 'b', 'T', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};
fieldNamesConfigs = {'Config1_Mon_Ton_Woff', 'Config2_Mon_Ton_Won', 'Config3_Mon_Toff_Won', 'Config4_Moff_Toff_Won','Config5_Moff_Toff_Woff'};
columns = 5:13; 
rows = [3 14 20 31 37 48 54 56 62 64];
for j = 1:length(fieldNamesConfigs)
    rowStart = rows(2*j - 1);
    rowEnd   = rows(2*j);
    for i = 1:length(fieldNamesdata)
        data.(fieldNamesConfigs{j}).(fieldNamesdata{i}) = table2array(dataRaw(rowStart:rowEnd, columns(i)));
        if i == 3
            data.(fieldNamesConfigs{j}).(fieldNamesdata{i}) = convtemp(data.(fieldNamesConfigs{j}).(fieldNamesdata{i}), 'F','K');
        elseif any(i == [4 5 6])
            data.(fieldNamesConfigs{j}).(fieldNamesdata{i}) = convforce(data.(fieldNamesConfigs{j}).(fieldNamesdata{i}), 'lbf','N');
        elseif any(i == [7 8 9])
            data.(fieldNamesConfigs{j}).(fieldNamesdata{i}) = (data.(fieldNamesConfigs{j}).(fieldNamesdata{i})).*0.112985;
        end
    end
    if any(i == [2 3])
        data.(fieldNamesConfigs{j}).Tamb = table2array(dataRaw(rowStart-3, 15));
        data.(fieldNamesConfigs{j}).Tamb = convtemp(data.(fieldNamesConfigs{j}).Tamb, 'F', 'K');
        data.(fieldNamesConfigs{j}).Pamb = table2array(dataRaw(rowStart-2, 15));
    else 
        data.(fieldNamesConfigs{j}).Tamb = table2array(dataRaw(rowStart-2, 15));
        data.(fieldNamesConfigs{j}).Tamb = convtemp(data.(fieldNamesConfigs{j}).Tamb, 'F', 'K');
        data.(fieldNamesConfigs{j}).Pamb = table2array(dataRaw(rowStart-1, 15));
    end
end
for i = 1:length(fieldNamesdata)
    data.(fieldNamesConfigs{4}).(fieldNamesdata{i}) = [ones(9,1)*data.(fieldNamesConfigs{4}).(fieldNamesdata{i})(1); data.(fieldNamesConfigs{4}).(fieldNamesdata{i})];
    data.(fieldNamesConfigs{5}).(fieldNamesdata{i}) = [ones(9,1)*data.(fieldNamesConfigs{5}).(fieldNamesdata{i})(1); data.(fieldNamesConfigs{5}).(fieldNamesdata{i})];
end

%% -- Global Variables --
% Wind Tunnel Test Setting
q = 7 * 248.84; % pa
% Full Model Airplane: DC-6B
S = 93.81 * 0.0006451606243503233; % reference area, m^2
cBar = convlength(3.466, 'in', 'm'); % reference length, in
b = convlength(27.066, 'in', 'm'); % reference wingspan, in
%% -- Calculate cl, cm, cn, cd --
% F = [F(model on, wind on) - F(model on, wind off)]- [F(model off, wind on) - F(model off, wind off)]
fxTOn = (data.Config2_Mon_Ton_Won.Fx - data.Config1_Mon_Ton_Woff.Fx) - (data.Config4_Moff_Toff_Won.Fx - data.Config5_Moff_Toff_Woff.Fx);
fzTOn = (data.Config2_Mon_Ton_Won.Fz - data.Config1_Mon_Ton_Woff.Fz) - (data.Config4_Moff_Toff_Won.Fz - data.Config5_Moff_Toff_Woff.Fz);
myTOn = (data.Config2_Mon_Ton_Won.My - data.Config1_Mon_Ton_Woff.My) - (data.Config4_Moff_Toff_Won.My - data.Config5_Moff_Toff_Woff.My);
mzTOn = (data.Config2_Mon_Ton_Won.Mz - data.Config1_Mon_Ton_Woff.Mz) - (data.Config4_Moff_Toff_Won.Mz - data.Config5_Moff_Toff_Woff.Mz);
fxTOff = (data.Config3_Mon_Toff_Won.Fx - data.Config1_Mon_Ton_Woff.Fx) - (data.Config4_Moff_Toff_Won.Fx - data.Config5_Moff_Toff_Woff.Fx);
fzTOff = (data.Config3_Mon_Toff_Won.Fz - data.Config1_Mon_Ton_Woff.Fz) - (data.Config4_Moff_Toff_Won.Fz - data.Config5_Moff_Toff_Woff.Fz);
myTOff = (data.Config3_Mon_Toff_Won.My - data.Config1_Mon_Ton_Woff.My) - (data.Config4_Moff_Toff_Won.My - data.Config5_Moff_Toff_Woff.My);
mzTOff = (data.Config3_Mon_Toff_Won.Mz - data.Config1_Mon_Ton_Woff.Mz) - (data.Config4_Moff_Toff_Won.Mz - data.Config5_Moff_Toff_Woff.Mz);
clTOn = fzTOn/(q*S);
clTOff = fzTOff/(q*S);
cdTOn = fxTOn/(q*S);
cdTOff = fxTOff/(q*S);
cmTOn = myTOn/(q*S*cBar);
cmTOff = myTOff/(q*S*cBar);
cnTOn = mzTOn/(q*S*b);
cnTOff = mzTOff/(q*S*b);
%% Calculate Oswald efficiency factor e, cl/cd max, dcl/dα, clmax and stalling angle of atack, dcm/dα, dcn/dbeta
AoA = data.Config1_Mon_Ton_Woff.a;
beta = data.Config1_Mon_Ton_Woff.b;
aspectRatio = b / cBar;
oswaldEfficiency = 1.78 * (1 - 0.045 * (aspectRatio ^ 0.68)) - 0.64;
% cl/cd max
cl_cdMaxOn = max(clTOn ./ cdTOn);
cl_cdMaxOff = max(clTOff ./ cdTOff);
% dcl/dalpha
dcl_dalphaOn = polyfit(AoA, clTOn, 1);
dcl_dalphaOff = polyfit(AoA, clTOff, 1);
% cl_max
clMaxOn = max(clTOn);
clMaxOff = max(clTOff);
% stall_angle (angle at cl_max)
stallAngleOn = AoA(clTOn == clMaxOn);
stallAngleOff = AoA(clTOff == clMaxOff);
% dcm/dalpha
dcm_dalphaOn = polyfit(AoA, cmTOn, 1);
dcm_dalphaOff = polyfit(AoA, cmTOff, 1);
% dcn/dbeta
dcn_dbetaOn = polyfit(beta, cnTOn, 1);
dcn_dbetaOff = polyfit(beta, cnTOff, 1);
% Visualization
f = figure;
range = 1:10;
figure (1)
plot(AoA(range), clTOn(range), '-o'); hold on
plot(AoA(range), clTOff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('\alpha'); ylabel('C_L')
title('C_L vs \alpha')
grid on
% Print_PNG(f, 'cl_vs_AoA.png')
f = figure;
plot(AoA(range), cmTOn(range), '-o'); hold on
plot(AoA(range), cmTOff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('\alpha'); ylabel('C_M')
title('C_M vs \alpha')
grid on
% Print_PNG(f, 'cm_vs_AoA.png')
f = figure;
plot([beta(11:12); beta(10); beta(13:end)], [cnTOn(11:12); cnTOn(10); cnTOn(13:end)], '-o'); hold on
plot([beta(11:12); beta(10); beta(13:end)], [cnTOff(11:12); cnTOff(10); cnTOff(13:end)], '-d')
legend('Tail on', 'Tail off')
xlabel('\beta'); ylabel('C_N')
title('C_N vs \beta')
grid on
% Print_PNG(f, 'cn_vs_beta.png')
f = figure;
plot(cdTOn(range), clTOn(range), '-o'); hold on
plot(cdTOff(range), clTOff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('C_D'); ylabel('C_L')
title('C_D vs. C_L')
grid on
% Print_PNG(f, 'cd_vs_cl.png')
% close all
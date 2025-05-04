clc; 
clearvars;
close all; 
clear all;
%%

dataRaw = readtable('Lab 5 Thursday SP25 Data Sheet.xlsx', 'Sheet','Ensemble Average FORCE DATA');
fieldNamesData = {'a', 'b', 'T', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};
fieldNamesConfigs = {'Config1_Mon_Ton_Woff', 'Config2_Mon_Ton_Won', 'Config3_Mon_Toff_Won', 'Config4_Moff_Toff_Won','Config5_Moff_Toff_Woff'};
columns = 5:13; 
rows = [3 16 21 34 39 52 57 61 65 69];
for j = 1:length(fieldNamesConfigs)
    for i = 1:length(fieldNamesData)
        data.(fieldNamesConfigs{j}).(fieldNamesData{i}) = table2array(dataRaw(rows(j):rows(j+1), columns(i)));
        % if i == 3
        %     data.(fieldNamesConfigs{j}).(fieldNamesData{i}) = convtemp(data.(fieldNamesConfigs{j}).(fieldNamesData{i}), 'F','K');
        % elseif any(i == [4 5 6])
        %     data.(fieldNamesConfigs{j}).(fieldNamesData{i}) = convforce(data.(fieldNamesConfigs{j}).(fieldNamesData{i}), 'lbf','N');
        % elseif any(i == [7 8 9])
        %     data.(fieldNamesConfigs{j}).(fieldNamesData{i}) = (data.(fieldNamesConfigs{j}).(fieldNamesData{i})).*0.112985;
        % end
    end
    data.(fieldNamesConfigs{j}).Tamb = table2array(dataRaw(rows(j)-2, 15));
    data.(fieldNamesConfigs{j}).Pamb = table2array(dataRaw(rows(j)-1, 15));
    data.(fieldNamesConfigs{j}).Pamb = Corr_Press(data.(fieldNamesConfigs{j}).Tamb, data.(fieldNamesConfigs{j}).Pamb);
    data.(fieldNamesConfigs{j}).Tamb = convtemp(data.(fieldNamesConfigs{j}).Tamb, 'F', 'K');
    rows = rows(2:end);
end
for i = 1:length(fieldNamesData)
    data.(fieldNamesConfigs{4}).(fieldNamesData{i}) = [ones(9,1)*data.(fieldNamesConfigs{4}).(fieldNamesData{i})(1); data.(fieldNamesConfigs{4}).(fieldNamesData{i})];
    data.(fieldNamesConfigs{5}).(fieldNamesData{i}) = [ones(9,1)*data.(fieldNamesConfigs{5}).(fieldNamesData{i})(1); data.(fieldNamesConfigs{5}).(fieldNamesData{i})];
end

%% -- Global Variables --
% Wind Tunnel Test Setting
q = 5 * 248.84; % pa
% Full Model Airplane: DC-6B
S = 93.81 * 0.0006451606243503233; % reference area, m^2
c_bar = convlength(3.466, 'in', 'm'); % reference length, in
b = convlength(27.066, 'in', 'm'); % reference wingspan, in
%% -- Calculate cl, cm, cn, cd --
% F = [F(model on, wind on) - F(model on, wind off)]- [F(model off, wind on) - F(model off, wind off)]
Fx_Ton = (Data.Config2_Mon_Ton_Won.Fx - Data.Config1_Mon_Ton_Woff.Fx) - (Data.Config4_Moff_Toff_Won.Fx - Data.Config5_Moff_Toff_Woff.Fx);
Fz_Ton = (Data.Config2_Mon_Ton_Won.Fz - Data.Config1_Mon_Ton_Woff.Fz) - (Data.Config4_Moff_Toff_Won.Fz - Data.Config5_Moff_Toff_Woff.Fz);
My_Ton = (Data.Config2_Mon_Ton_Won.My - Data.Config1_Mon_Ton_Woff.My) - (Data.Config4_Moff_Toff_Won.My - Data.Config5_Moff_Toff_Woff.My);
Mz_Ton = (Data.Config2_Mon_Ton_Won.Mz - Data.Config1_Mon_Ton_Woff.Mz) - (Data.Config4_Moff_Toff_Won.Mz - Data.Config5_Moff_Toff_Woff.Mz);
Fx_Toff = (Data.Config3_Mon_Toff_Won.Fx - Data.Config1_Mon_Ton_Woff.Fx) - (Data.Config4_Moff_Toff_Won.Fx - Data.Config5_Moff_Toff_Woff.Fx);
Fz_Toff = (Data.Config3_Mon_Toff_Won.Fz - Data.Config1_Mon_Ton_Woff.Fz) - (Data.Config4_Moff_Toff_Won.Fz - Data.Config5_Moff_Toff_Woff.Fz);
My_Toff = (Data.Config3_Mon_Toff_Won.My - Data.Config1_Mon_Ton_Woff.My) - (Data.Config4_Moff_Toff_Won.My - Data.Config5_Moff_Toff_Woff.My);
Mz_Toff = (Data.Config3_Mon_Toff_Won.Mz - Data.Config1_Mon_Ton_Woff.Mz) - (Data.Config4_Moff_Toff_Won.Mz - Data.Config5_Moff_Toff_Woff.Mz);
Cl_Ton = Fz_Ton/(q*S);
Cl_Toff = Fz_Toff/(q*S);
Cd_Ton = Fx_Ton/(q*S);
Cd_Toff = Fx_Toff/(q*S);
Cm_Ton = My_Ton/(q*S*c_bar);
Cm_Toff = My_Toff/(q*S*c_bar);
Cn_Ton = Mz_Ton/(q*S*b);
Cn_Toff = Mz_Toff/(q*S*b);
%% Calculate Oswald efficiency factor e, cl/cd max, dcl/dα, clmax and stalling angle of atack, dcm/dα, dcn/dbeta
AoA = data.Config1_Mon_Ton_Woff.a;
Beta = data.Config1_Mon_Ton_Woff.b;
aspect_ratio = b / c_bar;
oswald_efficiency = 1.78 * (1 - 0.045 * (aspect_ratio ^ 0.68)) - 0.64;
% cl/cd max
cl_cd_max_on = max(Cl_Ton ./ Cd_Ton);
cl_cd_max_off = max(Cl_Toff ./ Cd_Toff);
% dcl/dalpha
dcl_dalpha_on = polyfit(AoA, Cl_Ton, 1);
dcl_dalpha_off = polyfit(AoA, Cl_Toff, 1);
% cl_max
cl_max_on = max(Cl_Ton);
cl_max_off = max(Cl_Toff);
% stall_angle (angle at cl_max)
stall_angle_on = AoA(Cl_Ton == cl_max_on);
stall_angle_off = AoA(Cl_Toff == cl_max_off);
% dcm/dalpha
dcm_dalpha_on = polyfit(AoA, Cm_Ton, 1);
dcm_dalpha_off = polyfit(AoA, Cm_Toff, 1);
% dcn/dbeta
dcn_dbeta_on = polyfit(Beta, Cn_Ton, 1);
dcn_dbeta_off = polyfit(Beta, Cn_Toff, 1);
% Visualization
f = figure;
range = 1:10;
figure (1)
plot(AoA(range), Cl_Ton(range), '-o'); hold on
plot(AoA(range), Cl_Toff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('\alpha'); ylabel('C_L')
title('C_L vs \alpha')
grid on
Print_PNG(f, 'cl_vs_AoA.png')
f = figure;
plot(AoA(range), Cm_Ton(range), '-o'); hold on
plot(AoA(range), Cm_Toff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('\alpha'); ylabel('C_M')
title('C_M vs \alpha')
grid on
Print_PNG(f, 'cm_vs_AoA.png')
f = figure;
plot([Beta(11:12); Beta(10); Beta(13:end)], [Cn_Ton(11:12); Cn_Ton(10); Cn_Ton(13:end)], '-o'); hold on
plot([Beta(11:12); Beta(10); Beta(13:end)], [Cn_Toff(11:12); Cn_Toff(10); Cn_Toff(13:end)], '-d')
legend('Tail on', 'Tail off')
xlabel('\beta'); ylabel('C_N')
title('C_N vs \beta')
grid on
Print_PNG(f, 'cn_vs_beta.png')
f = figure;
plot(Cd_Ton(range), Cl_Ton(range), '-o'); hold on
plot(Cd_Toff(range), Cl_Toff(range), '-d')
legend('Tail on', 'Tail off')
xlabel('C_D'); ylabel('C_L')
title('C_D vs. C_L')
grid on
Print_PNG(f, 'cd_vs_cl.png')
close all
%% --- Functions --- %%
% Function for correcting barometer pressure
function Corrected_Pressure_Amb = Corr_Press(Tamb, Pamb)
% Latitude Correction
SD_lat = 32.7157;
Lat_correction29 = interp1([32 34],[0.035 0.030], SD_lat);
Lat_correction30 = interp1([32 34],[0.036 0.031], SD_lat);
Lat_correction = interp1([29 30], [Lat_correction29 Lat_correction30], Pamb, 'linear', 'extrap');
% Temperature Correction
temps = 72:2:78;
correction29 = [.114, .119, .124 .129];
correction30 = [.118, .123, .128 .134];
correction31 = [.122, .127, .133 .138];
T_correct29 = interp1(temps, correction29, Tamb);
T_correct30 = interp1(temps, correction30, Tamb);
T_correct31 = interp1(temps, correction31, Tamb);
Temp_correction = interp1([29 30 31], [T_correct29 T_correct30 T_correct31], Pamb);
% Total inHg correction
Corrected_h = Pamb - Lat_correction - Temp_correction;
Corrected_h = convlength(Corrected_h,'in','m'); % m
% Corrected Pressure Calculation
Standard_gravity = 9.80665; % m/s^2
rhoHG = 13595; % kg/m^3
Corrected_Pressure_Amb = rhoHG.*Standard_gravity*Corrected_h; % kg/m^3 * m/s^2 * m = Pa
end
function Print_PNG(f, title)
pdfFileName = sprintf(title);
print(f, pdfFileName, '-dpng');
end
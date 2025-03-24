function [Corrected_Pressure_Amb] = Correct_Pres(Tamb,Pamb)
    % Latitude Correction
    SD_lat = 32.7157;
    Lat_correction29 = interp1([32 34],[0.035 0.030], SD_lat);
    Lat_correction30 = interp1([32 34],[0.036 0.031], SD_lat);
    Lat_correction = interp1([29 30], [Lat_correction29 Lat_correction30], Pamb, 'linear','extrap');
    % Temperature Correction
    temps = 72:2:76;
    correction29 = [.114, .119, .124];
    correction30 = [.118, .123, .128];
    correction31 = [.122, .127, .133];
    T_correct29 = interp1(temps, correction29, Tamb);
    T_correct30 = interp1(temps, correction30, Tamb);
    T_correct31 = interp1(temps, correction31, Tamb);
    Temp_correction = interp1([29 30 31], [T_correct29 T_correct30 T_correct31], Pamb);
    % Total inHg correction
    Corrected_h = Pamb - Lat_correction - Temp_correction;
    Corrected_h = Corrected_h/39.37; % m
    % Corrected Pressure Calculation
    Standard_gravity = 9.80665; % m/s^2
    rhoHG = 13595; % kg/m^3
    Corrected_Pressure_Amb = rhoHG.*Standard_gravity*Corrected_h; % kg/m^3 * m/s^2 * m =Pa
    Corrected_Pressure_Amb = Corrected_Pressure_Amb/6895;
end

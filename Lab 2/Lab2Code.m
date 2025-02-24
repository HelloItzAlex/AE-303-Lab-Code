%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = N(1:35);
q01 = x0S1(1:35); q02 = x0S2(1:35); q03 = x0S3(1:35); q04 = x0S4(1:35);
q21 = x2S1(1:35); q2_7 = x2S4(1:35); q23 = x2S4(1:35); q24 = x2S4(1:35);
q71 = x7S1(1:35); q72 = x7S2(1:35); q73 = x7S3(1:35); q74 = x7S4(1:35);
pa = 30.02; % inHg
pa_Pa = pa*3386.39; % Pa
q22 = 2*248.84; %Pa
q77 = 7*248.84; %pa
T0 = 298.1223889; %K
T2 = 298.4571667; %K
T7 = 303.8085556; %K
R = 287.058;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q0 = (q01 + q02 + q03 + q04)/4;
q2 = (q21 + q2_7 + q23 + q24)/4;
q7 = (q71 + q72 + q73 + q74)/4;
q0d = q0 - q0(7);
q2d = q2 - q2(7);
q7d = q7 - q7(7);
mq0 = (sum(q0(1:6)) + sum(q0(8:35)))/34;
mq2 = (sum(q2(1:6)) + sum(q2(8:35)))/34;
mq7 = (sum(q7(1:6)) + sum(q7(8:35)))/34;
mq0d = (sum(q0d(1:6)) + sum(q0d(8:35)))/34;
mq2d = (sum(q2d(1:6)) + sum(q2d(8:35)))/34;
mq7d = (sum(q7d(1:6)) + sum(q7d(8:35)))/34;
%%%%%%%%%%%%%%%%%%%%% Measured Pressure Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(n,q0, 'o', 'color', "#A2142F");
yline(mq0,'LineWidth',1)
hold on
plot(n,q2, 'o', 'color', "#0072BD");
yline(mq2,'LineWidth',1)
hold on
plot(n,q7, 'o', 'color', "#77AC30");
yline(mq7,'LineWidth',1)
title('Flow Uniformity - Measured Pressure')
ylabel('Measured Pressure [psi]', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
xlim([0 36])
ylim([-0.05 0.25])
xticks(1:2:37)
yticks(-0.05:0.025:0.25)
legend('0.0 $H_{2}$O', '', '2.0 in$H_{2}$O', '', '7.0 in$H_{2}$O',...
 'Average', 'Location','northeastoutside', 'Interpreter','Latex');
grid on
%%%%%%%%%%%%%%%%%%%%% Dynamic Pressure Greph %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(n,q0d, 'o', 'color', "#A2142F");
yline(mq0d,'LineWidth',1);
hold on
plot(n,q2d, 'o', 'color', "#0072BD");
yline(mq2d, 'LineWidth',1)
yline(0.0721824, 'LineWidth',1, 'color', 'r')
hold on
plot(n,q7d, 'o', 'color', "#77AC30");
yline(mq7d,'LineWidth',1)
yline(0.252638, 'LineWidth',1, 'color', 'r')
title('Flow Uniformity - Dynamic Pressure')
ylabel('Dynamic Pressure [psi]', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
xlim([0 36])
ylim([-0.05 0.3])
xticks(1:2:37)
yticks(-0.05:0.025:0.3)
legend('0.0 $H_{2}$O', '', '2.0 in$H_{2}$O', '', '', '7.0 in$H_{2}$O',...
 'Average', '$q_{setting}$', 'Location','northeastoutside',...
 'Interpreter','Latex');
grid on
%%%%%%%%%%%%%%%%%%%%% Test Section Airspeed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p0 = pa_Pa + (q0(7)*6894.76 - q0(7)*6894.76);
p2 = pa_Pa + (q2(7)*6894.76 - q0(7)*6894.76);
p7 = pa_Pa + (q7(7)*6894.76 - q0(7)*6894.76);
d0 = p0/(R*T0);
d2 = p2/(R*T2);
d7 = p7/(R*T7);
q0_7 = q0d; q2_7 = q2d; q7_7 = q7d;
q0_7(7)=NaN; q2_7(7)=NaN; q7_7(7)=NaN;
Pa_q0 = q0_7*6894.76;
Pa_q2 = q2_7*6894.76;
Pa_q7 = q7_7*6894.76;
v0 = ((2*Pa_q0)/d0).^(1/2);
v2 = ((2*Pa_q2)/d2).^(1/2);
v7 = ((2*Pa_q7)/d7).^(1/2);
mv2 = (sum(v2(1:6)) + sum(v2(8:35)))/34;
mv7 = (sum(v7(1:6)) + sum(v7(8:35)))/34;
v22 = ((2*q22)/d2).^(1/2);
v77 = ((2*q77)/d7).^(1/2);
%%%%%%%%%%%%%%%%%%% Test Section Airspeed Graph %%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(n,v2, 'o', 'color', "#0072BD");
yline(mv2,'LineWidth',1)
yline(v22, 'LineWidth',1, 'color', 'r');
hold on
plot(n,v7, 'o', 'color', "#77AC30");
yline(mv7,'LineWidth',1)
yline(v77, 'LineWidth',1, 'color', 'r')
title('Flow Uniformity - Test Section Airspeed')
ylabel('Test Section Airspeed [m/$s^{2}$]', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
xlim([0 36])
ylim([20 60])
xticks(1:2:37)
yticks(20:2.5:60)
legend('2.0 in$H_{2}$O', '', '', '7.0 in$H_{2}$O','Average',...
 '$v_{setting}$', 'Location','northeastoutside',...
 'Interpreter','Latex');
grid on
%%%%%%%%%%%%%%%%%%% Pressure/Velocity Deviation %%%%%%%%%%%%%%%%%%%%%%%%%%
dq2 = (q2_7-mq2d)/mq2d*100; dq7 = (q7_7-mq7d)/mq7d*100;
dv2 = (v2-mv2)/mv2*100; dv7 = (v7-mv7)/mv7*100;
%%%%%%%%%%%%%%%%%% Pressure Deviation Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Combined %%%%%
figure (4)
plot(n,dq2, 'o', 'color', "#0072BD");
hold on
plot(n,dq7, 'diamond', 'color', "#77AC30");
title('Flow Uniformity - Dynamic Pressure Deviation From Average')
ylabel('$\frac{q-\bar{q}}{\bar{q}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-4 3])
xticks(1:2:37)
yticks(-4:0.5:3)
legend('2.0 in$H_{2}$O', '7.0 in$H_{2}$O', 'Location','southeast',...
 'Interpreter','Latex');
grid on
%%%%% 2 Only %%%%%%
figure (5)
plot(n,dq2, 'o', 'color', "#0072BD");
title('Flow Uniformity - Dynamic Pressure Deviation From Average')
subtitle('2.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('$\frac{q-\bar{q}}{\bar{q}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-4 3])
xticks(1:2:37)
yticks(-4:0.5:3)
grid on
%%%%% 7 Only %%%%%
figure (6)
plot(n,dq7, 'diamond', 'color', "#77AC30");
title('Flow Uniformity - Dynamic Pressure Deviation From Average')
subtitle('7.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('$\frac{q-\bar{q}}{\bar{q}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-4 3])
xticks(1:2:37)
yticks(-4:0.5:3)
grid on
%%%%%%%%%%%%%%%%%%% Velocity Deviation Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (7)
plot(n,dv2, 'o', 'color', "#0072BD");
hold on
plot(n,dv7, 'diamond', 'color', "#77AC30");
title('Flow Uniformity - Airspeed Deviation From Average')
ylabel('$\frac{v-\bar{v}}{\bar{v}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-2 1.5])
xticks(1:2:37)
yticks(-2:0.25:1.5)
legend('2.0 in$H_{2}$O', '7.0 in$H_{2}$O', 'Location','southeast',...
 'Interpreter','Latex');
grid on
%%%%% 2 Only %%%%%%
figure (8)
plot(n,dv2, 'o', 'color', "#0072BD");
title('Flow Uniformity - Airspeed Deviation From Average')
subtitle('2.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('$\frac{v-\bar{v}}{\bar{v}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-2 1.5])
xticks(1:2:37)
yticks(-2:0.25:1.5)
grid on
%%%%% 7 Only %%%%%%
figure (9)
plot(n,dv7, 'diamond', 'color', "#77AC30");
title('Flow Uniformity - Airspeed Deviation From Average')
subtitle('7.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('$\frac{v-\bar{v}}{\bar{v}}$', 'Interpreter','Latex')
xlabel('Pressure Port Number', 'Interpreter','Latex')
ytickformat('percentage')
xlim([0 36])
ylim([-2 1.5])
xticks(1:2:37)
yticks(-2:0.25:1.5)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%% Raynald's Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ku2 = 0.000018387;
ku7 = 0.0000186397;
Re2 = d2*mv2/ku2;
Re7 = d7*mv7/ku7;
%%%%%%%%%%%%%%%%%%%%%%%%% Cotour Graph Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [-19.5, -16.5, -13.5, -10.5, -7.5, -4.5, -1.5, 0, 1.5, 4.5,...
 7.5, 10.5, 13.5, 16.5, 19.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5,...
 -7.5, -7.5, -7.5, -7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5,...
 7.5, 7.5,]';
z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.5, 10.5, 7.5,...
 4.5, 1.5, -1.5, -4.5, -7.5, -10.5, -13.5, 13.5, 10.5, 7.5, 4.5,...
 1.5, -1.5, -4.5, -7.5, -10.5, -13.5,]';
dq2(7) = dq2(9); dq7(7) = dq7(9);
dq_q2 = dq2; dq_q7 = dq7;
v2(7) = v2(9); v7(7) = v7(9);
dv_v2 = v2; dv_v7 = v7;
x1 = x; z1 = z;
x1(7) = NaN; z1(7) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure 2.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
[X,Z] = meshgrid(linspace(min(x),max(x)), linspace(min(z),max(z)));
F = scatteredInterpolant(x,z,dq_q2,'linear','none');
Dq = F(X,Z);
contourf(X,Z,Dq)
hold on
plot(x1,z1, 'o', 'color', "r");
c = colorbar;
c.Ruler.TickLabelFormat='%g%%';
title('Dynamic Pressure Deviation from Average')
subtitle('q = 2.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('z Position [in]', 'Interpreter','Latex')
xlabel('x Position [in]', 'Interpreter','Latex')
ylabel(c,'$\frac{v-\bar{v}}{\bar{v}}$','Position',[4.25 -.75],...
 'FontSize',18,'Rotation',270,...
 'Interpreter','Latex');
xlim([-20 20])
ylim([-15 15])
xticks(-19.5:3:19.5)
yticks(-13.5:3:13.5)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure 7.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
[X,Z] = meshgrid(linspace(min(x),max(x)), linspace(min(z),max(z)));
F = scatteredInterpolant(x,z,dq_q7,'linear','none');
Dq = F(X,Z);
contourf(X,Z,Dq)
hold on
plot(x1,z1, 'o', 'color', "r");
c = colorbar;
c.Ruler.TickLabelFormat='%g%%';
title('Dynamic Pressure Deviation from Average')
subtitle('q = 7.0 in$H_{2}$O', 'Interpreter','Latex')
ylabel('z Position [in]', 'Interpreter','Latex')
xlabel('x Position [in]', 'Interpreter','Latex')
ylabel(c,'$\frac{v-\bar{v}}{\bar{v}}$','Position',[4.25 -.3],...
 'FontSize',18,'Rotation',270,...
 'Interpreter','Latex');
xlim([-20 20])
ylim([-15 15])
xticks(-19.5:3:19.5)
yticks(-13.5:3:13.5)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%% Airspeed 2.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12)
[X,Z] = meshgrid(linspace(min(x),max(x)), linspace(min(z),max(z)));
F = scatteredInterpolant(x,z,dv_v2,'linear','none');
Dq = F(X,Z);
contourf(X,Z,Dq)
hold on
plot(x1,z1, 'o', 'color', "r");
c = colorbar;
title('Airspeed Distrobution')
subtitle('IAS = 28.97 m/s', 'Interpreter','Latex')
ylabel('z Position [in]', 'Interpreter','Latex')
xlabel('x Position [in]', 'Interpreter','Latex')
ylabel(c,'Airspeed [m/s]','Position',[3.5 29.28],...
 'FontSize',11,'Rotation',270,...
 'Interpreter','Latex');
xlim([-20 20])
ylim([-15 15])
xticks(-19.5:3:19.5)
yticks(-13.5:3:13.5)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%% Airspeed 2.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)
[X,Z] = meshgrid(linspace(min(x),max(x)), linspace(min(z),max(z)));
F = scatteredInterpolant(x,z,dv_v7,'linear','none');
Dq = F(X,Z);
contourf(X,Z,Dq)
hold on
plot(x1,z1, 'o', 'color', "r");
c = colorbar;
title('Airspeed Distrobution')
subtitle('IAS = 54.74 m/s', 'Interpreter','Latex')
ylabel('z Position [in]', 'Interpreter','Latex')
xlabel('x Position [in]', 'Interpreter','Latex')
ylabel(c,'Airspeed [m/s]','Position',[3.5 55.4],...
 'FontSize',11,'Rotation',270,...
 'Interpreter','Latex');
xlim([-20 20])
ylim([-15 15])
xticks(-19.5:3:19.5)
yticks(-13.5:3:13.5)
grid on
%% Load parameters
parameters

%% Internal Model     % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     dX =  J(X, Y) - K.*X - phi.*U;
%     dY = -J(X, Y);
%% Find equilibrium solution
Xequil = -phi*U/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);
% [0.1696664 6.181387558794222]
%% Figure 1A
kick = [.001 0]; % increase cytosolic Ca by 1e-3 uM
Tkick = 2;
Tfinal = 20;
[T1,Y1] = ode45(@InternalModel,[0 Tkick],[Xequil Yequil]);
[T2,Y2] = ode45(@InternalModel,[T1(end) Tfinal],Y1(end,:) + kick);

T = [T1; T2];
Y = [Y1; Y2];
%%
clf('reset')
subplot(1,2,1)
plot(T, Y(:,1), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)' );

subplot(1,2,2)
plot(T, Y(:,2), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Stored Ca (\muM)');





%% Single Cell Model  % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%       dX =  J(X, Y) - K.*X - phi.*ICa(V);
%       dY = -J(X, Y);
%       dV = -1/C * (ICa(V) + IK(V) + Ileak(V));
%% Find equilibrium solution
syms Vsym
sol = solve(-1/C * (ICa(Vsym) + IK(-phi*ICa(Vsym)/K, Vsym) + Ileak(Vsym)));
Vequil = double(sol);
Xequil = -phi*ICa(Vequil)/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);
% [ -59.000020714943, 0.16999759457746, 6.17997307521213216]
%% Figure 1B
kick = [0 .1 0]; % increase cytosolic Ca by 0.1 uM
Tkick = 2;
Tfinal = 10;
options = odeset('RelTol',1e-5,'AbsTol',[1e-4 1e-4 1e-4]);
[T1,Y1] = ode15s(@SingleCellModel,[0 Tkick],[Vequil Xequil Yequil], options);
[T2,Y2] = ode15s(@SingleCellModel,[T1(end) Tfinal],Y1(end,:)+kick, options);

T = [T1; T2];
Y = [Y1; Y2];

clf('reset')
subplot(1,2,2)
plot(T, Y(:,1), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');

subplot(1,2,1)
plot(T, Y(:,2), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');






%% Single Cell Model with shunt % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%       dX =  J(X, Y) - K.*X - phi.*ICa(V);
%       dY = -J(X, Y);
%       dV = -1/C * (ICa(V) + IK(V) + Ileak(V));
% Increasing the shunt (leak) conductance causes oscillations to reappear
gleak = 2701; % The default value
gleak = gleak + 2e4; % gleak increases from 2701 to 22701 uS/cm^2
Ileak  = @(V) gleak * (V - Vleak);
%% Find equilibrium solution
syms Vsym
sol = solve(-1/C * (ICa(Vsym) + IK(-phi*ICa(Vsym)/K, Vsym) + Ileak(Vsym)));
Vequil = double(sol);
Xequil = -phi*ICa(Vequil)/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);
SingleCellModel(0, [Vequil, Xequil, Yequil])
% [ -55.681795622534693,   0.232706368015332,   5.723197194336485]
%% Figure 1C                                        % Already starts to oscillate without any kick . . .
Ileak  = @(V) gleak * (V - Vleak);

kick = [0 0 0]; % increase cytosolic Ca by 1e-3 uM
Tkick = 2;
Tfinal = 3;
options = odeset('RelTol',1e-3,'AbsTol',[1e-10 1e-10 1e-10]);
[T1,Y1] = ode45(@SingleCellModel,[0 Tkick],[Vequil Xequil Yequil], options);
[T2,Y2] = ode45(@SingleCellModel,[T1(end) Tfinal],Y1(end,:)+kick, options);

T = [T1; T2];
Y = [Y1; Y2];

clf('reset')
subplot(2,2,1)
plot(T, Y(:,1))
title( 'Membrane Potential [mV]' );

subplot(2,2,2)
plot(T, Y(:,2))
title( 'Cytosolic Ca [uM]' );

subplot(2,2,3)
plot(T, Y(:,3))
title( 'Stored Ca [uM]' );






%% Shunt Attempt 2 % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%       dX =  J(X, Y) - K.*X - phi.*ICa(V);
%       dY = -J(X, Y);
%       dV = -1/C * (ICa(V) + IK(V) + Ileak(V) + Ishunt(V));
% Add a high conductance with a reversal potential equal to the resting
% potential
%% Find equilibrium solution
global Vshunt gshunt Ishunt
Vshunt = -59;
% High shunt conductance
gshunt = 2e4;
Ishunt = @(V) gshunt * (V - Vshunt);

syms Vsym
sol = solve(-1/C * (ICa(Vsym) + IK(-phi*ICa(Vsym)/K, Vsym) + Ileak(Vsym) + Ishunt(Vsym)));
Vequil = double(sol);
Xequil = -phi*ICa(Vequil)/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);
% [ -59.000020714943368,   0.169997594577456,   6.179973075212130]
% Figure 1C                                      
kick = [0 1e-3 0]; % increase cytosolic Ca by 1e-3 uM
Tkick = 2;
Tfinal = 20;
options = odeset('RelTol',1e-3,'AbsTol',[1e-6 1e-4 1e-4]);
[T1,Y1] = ode15s(@SingleCellModel_Shunt,[0 Tkick],[Vequil Xequil Yequil], options);
[T2,Y2] = ode15s(@SingleCellModel_Shunt,[T1(end) Tfinal],Y1(end,:)+kick, options);

T = [T1; T2];
Y = [Y1; Y2];

clf('reset')
subplot(1,2,2)
plot(T, Y(:,1), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');

subplot(1,2,1)
plot(T, Y(:,2), 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');






%% Two Cell Model % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     dV1 = -1/C * (ICa(V1) + IK(X1,V1) + Ileak(V1) + Icoupled(V1, V2));
%     dX1 =  J(X1, Y1) - K.*X1 - phi.*ICa(V1);
%     dY1 = -J(X1, Y1);
%     
%     dV2 = -1/C * (ICa(V2) + IK(X2,V2) + Ileak(V2) + Icoupled(V2, V1));
%     dX2 =  J(X2, Y2) - K.*X2 - phi.*ICa(V2);
%     dY2 = -J(X2, Y2);
%% Equilibrium Solution
% Will be the same as in the single cell case with the values the same for
% both cells
% [ -59.000020714943, 0.16999759457746, 6.17997307521213216]
% global gcoupled Icoupled
% gcoupled = 10000;
% Icoupled  = @(V1, V2) gcoupled * (V1 - V2);
syms Vsym
sol = solve(-1/C * (ICa(Vsym) + IK(-phi*ICa(Vsym)/K, Vsym) + Ileak(Vsym)));
Vequil = double(sol);
Xequil = -phi*ICa(Vequil)/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);

V1equil = Vequil;
X1equil = Xequil;
Y1equil = Yequil;

V2equil = Vequil;
X2equil = Xequil;
Y2equil = Yequil;
%% 
kick = [0 .1 0 0 0 0];
Tkick = 1;
Tfinal = 15;
options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[T1,Y1] = ode15s(@TwoCoupledCellsModel,[0 Tkick],[V1equil X1equil Y1equil V2equil X2equil Y2equil], options);
[T2,Y2] = ode15s(@TwoCoupledCellsModel,[T1(end) Tfinal],Y1(end,:)+kick, options);

T = [T1; T2];
Y = [Y1; Y2];
%% Figure 3
clf('reset')
subplot(1,2,1)
StackedPlot(T, Y(:,[2 5]), {[10 1.6], 0.2, '\muM'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');

subplot(1,2,2)
StackedPlot(T, Y(:,[1 4]), {[10 -51], 2, 'mV'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');


parameters;
%% Find equilibrium solution
syms Vsym
sol = solve(-1/C * (ICa(Vsym) + IK(-phi*ICa(Vsym)/K, Vsym) + Ileak(Vsym)));
Vequil = double(sol);
Xequil = -phi*ICa(Vequil)/K;
Yequil = V2(Xequil)/(V3(Xequil) + K5);






%% Two Cell Model
Ncells = 2;

kick = [0 .1 0 0 0 0];

len = 3*Ncells;
initialValues = zeros( len, 1 );

initialValues(1:3:len) = Vequil;
initialValues(2:3:len) = Xequil;
initialValues(3:3:len) = Yequil;

% Run the dynamics of the system
Tkick = 1;
% Length of time to allow for equilibration
Tequilibrate = 30;
% Length of time to record the actual data
Trecord = 10;
options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[T1,Y1] = ode15s(@NetworkModel,[0 Tkick], initialValues, options);
endValues = Y1(end,:);
endValues = endValues + kick;
% Skip some period of time to allow the system to reach equilbrium
[T2,Y2] = ode15s(@NetworkModel,[0 Tequilibrate], endValues', options);
endValues = Y2(end,:);
[T,Y] = ode15s(@NetworkModel,[0 Trecord], endValues', options);
%% Generate Figure 3
clf('reset')
subplot(1,2,1)
StackedPlot(T, Y(:,[2 5]), {[7 1.5], 0.2, '\muM'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');

subplot(1,2,2)
StackedPlot(T, Y(:,[1 4]), {[7 -51], 1, 'mV'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');







%% 6 Cell Model
Ncells = 6;
% For kick_4A let Tequilibrate >= 26
kick_4A = [0.0, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
kick_4B = [0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0];

len = 3*Ncells;
initialValues = zeros( len, 1 );

initialValues(1:3:len) = Vequil;
initialValues(2:3:len) = Xequil;
initialValues(3:3:len) = Yequil;
%% Run the dynamics of the system
%kick = zeros(1, len);
%kick(2:3:len) = .1.*(sign(round(randn(1,Ncells))));
%kick(2:3:len) = [0 -.1 .1 -.1 .1 -.1];
kick = kick_4B;

Tkick = 1;
% Length of time to allow for equilibration
Tequilibrate = 30;
% Length of time to record the actual data
Trecord = 10;
options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[T1,Y1] = ode15s(@NetworkModel,[0 Tkick], initialValues, options);
endValues = Y1(end,:);
endValues = endValues + kick;
% Skip some period of time to allow the system to reach equilbrium
[T2,Y2] = ode15s(@NetworkModel,[0 Tequilibrate], endValues', options);
endValues = Y2(end,:);
[T,Y] = ode15s(@NetworkModel,[0 Trecord], endValues', options);

%% Display Stacked Plots for Figures A and B
Vmembrane = Y(:, 1:3:end );
cytoCa = Y(:, 2:3:end );

% Reordering for figure 4A
%cytoCa = Y(:, [8 14 17 2 5 11] );
%Vmembrane = Y(:, [8 14 17 2 5 11]-1 );

clf('reset')
subplot(1,2,1)
StackedPlot(T, Vmembrane, {[7 -37], 2, 'mV'}) % -45 for A  -37 for B
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');

subplot(1,2,2)
StackedPlot(T, cytoCa, {[7 4.8], .5, '\mu M'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');





%% Multiple Kicks for Figure 4C
kick = zeros(1, len);
kick(2) = .1;  % kick to the cytosolic ca concentration of cell 1\

kickStart = zeros(1, len);
kickStart(2:3:len) = .1.*(sign(round(randn(1,Ncells))));
kickStart = [0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0];

% Time between kicks
Tkick = 10;
numKicks = 6;
% Time to allow the system to equillibrate
Tequilibrate = 20;

options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[T,Y] = ode15s(@NetworkModel,[0 Tequilibrate], initialValues+kickStart', options);
endValues = Y(end,:);
% Start actually taking the data
T = [];
Y = [];
for i = 1: numKicks
    [Ttmp,Ytmp] = ode15s(@NetworkModel,[Tkick*(i-1) Tkick*i], endValues', options);
    % Store the values
    T = [T; Ttmp];
    Y = [Y; Ytmp];
    % Get the end values and apply the kick
    endValues = Ytmp(end,:) + kick;
end

%% Display the figure
Vmembrane = Y(:, 1 );

clf('reset')
plot(T, Vmembrane, 'LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');
%% Stacked Figure
Vmembrane = Y(:, 1:3:end );
cytoCa = Y(:, 2:3:end );
storCa = Y(:, 3:3:end );

clf('reset')
subplot(1,2,1)
StackedPlot(T, Vmembrane, {[7 -37], 2, 'mV'}) % -45 for A  -37 for B
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Membrane Potential (mV)');

subplot(1,2,2)
StackedPlot(T, cytoCa, {[7 4.8], .5, '\mu M'})
set(gca, 'FontSize', 12)
xlabel( 'time (s)' );
ylabel( 'Cytosolic Ca (\muM)');










%% Strongly Coupled Network  % % % %  % % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%
%
% Membrane potential is synchronized for all cells.
tic
Ncells = 200;

len = 1 + 2*Ncells;
initialValues = zeros( len, 1 );

initialValues(1)       = Vequil;
initialValues(2:2:len) = Xequil;
initialValues(3:2:len) = Yequil;

kick = zeros(1, len);
kick(2:2:len) = .1.*(sign(round(2*randn(1,Ncells))));

Tkick = 1;
% Length of time to allow for equilibration
Tequilibrate = 30;
% Length of time to record the actual data
Trecord = 20;

options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[T1,Y1] = ode15s(@NetworkModel_strong, [0 Tkick], initialValues, options);
endValues = Y1(end,:);
endValues = endValues + kick;
% Skip some period of time to allow the system to reach equilbrium
[T2,Y2] = ode15s(@NetworkModel_strong, [0 Tequilibrate], endValues', options);
endValues = Y2(end,:);
[T,Y] = ode15s(@NetworkModel_strong, [0 Trecord], endValues', options);
toc

%% Display the results
numDisplay = 25;

Vmembrane = Y(:, 1 );
cytoCa = Y(:, 2:2:end );

cytoCa = cytoCa(:, 1:25);

clf('reset')
subplot(2,1,1)
plot(T, Vmembrane, 'LineWidth', 2)
ylabel( 'Membrane Potential (mV)' );

subplot(2,1,2)
StackedPlot(T, cytoCa)
ylabel( 'Cytosolic Ca (uM)' );
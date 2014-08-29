% global namespace pollution!!!!!!!!!!!!!
global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
global V2 V3 J ICa IK Ileak h_inf m_inf sigma
% All equation parameters
%   Taken from a reference assuming a constnat inositol triphosphate
%   concentration of 0.14 microMolar

% Parameters for the Calcium Dynamics Cell Model
K  = 10; %Hz
K2 = 0.2; % microMolar
K4 = 0.69; % microMolar
K5 = 1; % Hz

Vm2 = 50; % microMolar Hz
Vm3 = 600; % Hz

phi = 9.221e-3; % conversion factor


U = -184; % nA/cm^2 % the constant calcium current

% Parameters with membrane potential interaction included
%   Taken form in vitro measurements in the inferior olive assuming
%   instantaneous activation.

% Voltages in mV and conductances in microsiemens/cm^2
% concentrations in microMolar

% internal variables
Vm = -61;
Tm = 4.2;

Vh = -85.5;
Th = 8.6;

% calcium current
VCa = 120;
gCa = 100;

% potassium current
VK = -85;
gK = 2000;

% leak current
Vleak = -55;
gleak = 2701;

C = 1; % microFarads/cm^2

beta = 2.5;
Xcrit = 0.4334;


% Functions in the dynamics models

% Calcium Dynamics
V2 = @(X) Vm2*(X.*X)./(K2*K2 + X.*X);
V3 = @(X) Vm3*(K4.*X).^3./(X + K4).^6;

J = @(X,Y) -V2(X) + ( V3(X) + K5).*Y;

% Membrane potential dependent Calcium current
m_inf = @(V) 1./(1 + exp(-(V - Vm)./Tm));
h_inf = @(V) 1./(1 + exp((V - Vh)./Th));

sigma = @(X) .5*(1 + tanh( beta.*(X - Xcrit) ) );

ICa =    @(V)   gCa * m_inf(V).^3 .* h_inf(V) .* (V - VCa);
IK  =    @(X,V) gK * sigma(X) .* (V - VK);
Ileak  = @(V)   gleak * (V - Vleak);


% 303 : Hodgkin and Huxley

% Here we will modify the passive membrane code to integrate the equations for the
% m, n and h variables.
%
% For later use, we'll retain the equation describing a passive membrane
% driven by a current source.
%
% parallel RC circuit
%
% I_R = V / R
% I_C = C dV/dt
% --> C dV/dt = I - V/R
%
%  Choose some reasonable constants
%
C = 1e-6; % 1 uF (per cm^2)
R = 1e4; % 10 kOhms (per cm^2)
%
%  Establish the time base (in sec)
dt = 0.0001;
tFinal = 0.1;
t = 0:dt:tFinal;        % this creates a list of times, from 0 to tFinal in steps of dt.
%
% Invent a current input and solve for the voltage output
%
% By commenting out one or the other, could choose either the step current 
% we looked at in class:
%
%I = (t>0.01 & t<0.02) * 1e-6;           % this is a current in nanoamps
% 
% or a sinusoidal current: 
%
I = 6*ones(size(t));
%
V = zeros(size(I));
tau = R*C % 10 msec = 0.01 sec
%
% Integrate equations using Euler method
%
% In the Euler method, we write:
%
%   dV/dt = [V(t+dt) - V(t)]/dt = -V(t)/RC  + I(t)/C
%
% So 
%
%   dV = V(t+dt) - V(t) =  dt* [-V(t)/RC + I(t)/C]
% 
% and then
%
%   V(t+dt) = V(t) + dV = V(t) + dt*[-V(t)/R + I(t)]/C.
%
% So the integration simply steps through time to whatever point you want
% to reach.

V(1) = 0;   % Initialize voltage to zero.
%
for i=1:(length(t)-1)
    dV = dt * (I(i) - V(i)/R) / C;
    V(i+1) = V(i) + (dV/dt) * dt;
end

% 
% Plot results
%
subplot(2,1,1);
plot(t, 1e6 * I);
%ylim([-0.05 1.05]);
xlabel('time (msec)')
ylabel('current (\mu A)')
subplot(2,1,2);
plot(10^3*t, 1e3 * V);
xlabel('time (msec)')
ylabel('voltage (mV)')
% 
V(1) = 0;   % Initialize voltage to zero.
%
% How to drive with a conductance input?
%
% Now we need to replace I(t) by g(t)*(V(t)-Vg), where now you control g(t). 
%   Give it a try!


%===========================================================================
%
% Voltage clamp experiments on m, n and h.
%
% Now we are solving a similar first order equation, but where the time
% constant depends on V. 
%
% Note that since we are starting with the case that we are doing a voltage
% clamp experiment,  V is not one of the dynamical variables and has no
% corresponding equation; we just input it, as we did the current in the 
% example above.
%
% Let's write out the equation for m:
%
% We can rearrange the HH form 
%
%   dm/dt = alpha_m (1 - m) - beta_m m
% 
% to give
%  
%   dm/dt = -m *(alpha_m + beta_m) + alpha_m
%
% Now the Euler method gives us (making explicit the V dependence of alpha and beta):
% 
% dm = m(d+dt) - m(t) = dt*[-m(t)*(alpha(V(t)) + beta(V(t))) + alpha(V(t))]
%
%
% Note that the rate constants are all in units of milliseconds!
% so let's establish a time base in msec
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now, let's proceed to put in the V dynamics.
%
%
dt = 0.05;
tFinal = 500;
t = 0:dt:tFinal;
% 
clear V;
%
% Integrate equations using Euler method
%
%
m = zeros(size(t));   % Initialize everything to zero.
n = zeros(size(t));
h = zeros(size(t));
V = zeros(size(t));
n(1) = 0.6;
h(1) = 0.3;
%
% Define constants
%
gkmax = 36; %*0.01 * 1000 ; % cm--> 10^7 g --> 10^-2 dt --> 10^3
gnamax = 120;
gl = 0.3;
%Vk = -77;
%Vl = -54.387;
%Vna = 50;
Vk = 12;
Vna = -115;
Vl = -10.613;
C = 1.0;
%
% Set up your experiment with your desired current input:
%
Iext = zeros(size(t));
Iext(500:end) = -5;   % don't forget to make the current negative...
%
% Now integrate to obtain the inactivation variables...
%
for i=1:(length(t)-1)
    alpham = 0.1*(V(i)+25)/(exp((V(i)+25)/10)-1);
    betam = 4*exp(V(i)/18);
    alphan = 0.01*(V(i)+10)/(exp((V(i)+10)/10)-1);
    betan = 0.125*exp(V(i)/80);
    alphah = 0.07*exp(V(i)/20);
    betah = 1./(exp((V(i)+30)/10)+1);
    dm = dt * (-m(i)*(alpham + betam) + alpham);
    m(i+1) = m(i) + dm;
    dn = dt * (-n(i)*(alphan + betan) + alphan);
    n(i+1) = n(i) + dn;
    dh = dt * (-h(i)*(alphah + betah) + alphah);
    h(i+1) = h(i) + dh;
    Ik = gkmax*n(i)^4*(V(i) - Vk);
    Ina = gnamax*m(i)^3*h(i)*(V(i)-Vna);
    Il = gl*(V(i) - Vl);
    V(i+1) = V(i) + dt*(Iext(i) - Ik - Ina - Il)/C;
end

figure;
subplot(2,1,1);
plot(t,Iext)
xlabel('time (msec)')
ylabel('input current (nA)')
subplot(2,1,2);
plot(t,-V)
xlabel('time (msec)')
ylabel('voltage (mV)')

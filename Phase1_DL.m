clear;clc;clear

%% Time %%
% The time step is set, time vector is created and
% the length of the vector is taken for the later computation.

ts  =  input('What would you want the time step to be in millisecond?\n'); % Time-step (ms)

tt  =  100;  % Total simulation time (ms)
t_vec  =  0 : ts : tt;  % Time vector
len  =  length(t_vec);  % The length of the time vector


%% Injected Current %%
% The injected current step function is generated using for-loop for the
% duration given.

Iinj_comp  =  input('What is the magnitude (uA/cm^2) and duration (ms) of the injected current?\ne.g.)[5 0.5]\n');
% The injected current input is given as a vector

mag_Iinj  =  Iinj_comp(1);  % Magnitude of the injected current (uA/cm^2)
dur_Iinj  =  Iinj_comp(2);  % Duration of the injected current (ms)

Iinj  =  zeros(1,len);   % Creates the same size of a vector as t_range

for ll=1:floor(dur_Iinj/ts)  
    
    % Feeds the values of the step into I_inj vector for the duration
    % that's given
    
    Iinj(ll)  =  mag_Iinj;
    
end



%% Constants %%
% Constants are defined at this stage. 


gK_bar  =  36;   % The maximum conductance value of K (mS/cm^2)
gNa_bar  =  120; % The maximum conductance value of Na (mS/cm^2)
gL_bar  =  0.3;  % The maximum conductance value of L (mS/cm^2)

Erest  =  -70;   % Rest membrane voltage (mV)

EK  =  -12;      % Nernst potential of K (mV)
ENa  =  115;     % Nernst potential of Na (mV)
EL  =  10.6;     % Nernst potential of L (mV)

Cm  =  1.0;      % Membrane capacitance (uF/cm^2)


%% Initialization and Gating Variables %%
% The membrane voltage is initialized as well as m, n and h. Also, the
% gating variables are computed. The gating variables are initialized to
% find initial conditions of m, n and h.


Vm(1)  =  0;  % Displacement of membrane voltage initialization (mV)

am  =  0.1 * ((25 - Vm(1)) / (exp((25 - Vm(1)) / 10) - 1));  % Gating variable
bm  =  4 * exp(-Vm(1) / 18);                                 % Gating variable
an  =  0.01 * ((10 - Vm(1)) / (exp((10 - Vm(1)) / 10) - 1)); % Gating variable
bn  =  0.125 * exp(-Vm(1) / 80);                             % Gating variable
ah  =  0.07 * exp(-Vm(1) / 20);                              % Gating variable
bh  =  1 / (exp((30 - Vm(1)) / 10) + 1);                     % Gating variable

m(1)  =  am / (am + bm);  % Initialization of m
n(1)  =  an / (an + bn);  % Initialization of n
h(1)  =  ah / (ah + bh);  % Initialization of h


%% Solving Equations %%
% The differential equations of m, n, h and membrane voltage (Vm) are
% solved using Euler's method. Every iteration, the gating variables,
% conductances and currents are computed and fed into the differential
% equation.


for kk  =  1 : len - 1
    
    am  =  0.1 * ((25 - Vm(kk)) / (exp((25 - Vm(kk)) / 10) - 1));   % Gating variable
    bm  =  4 * exp(-Vm(kk) / 18);                                   % Gating variable
    an  =  0.01 * ((10 - Vm(kk)) / (exp((10 - Vm(kk)) / 10) - 1));  % Gating variable
    bn  =  0.125 * exp(-Vm(kk) / 80);                               % Gating variable
    ah  =  0.07 * exp(-Vm(kk) / 20);                                % Gating variable
    bh  =  1 / (exp((30 - Vm(kk)) / 10) + 1);                       % Gating variable
    
    gNa(kk)  =  m(kk)^3 * gNa_bar * h(kk);  % Conductance of Na (mS/cm^2)
    gK(kk)  =  n(kk)^4 * gK_bar;            % Conductance of K  (mS/cm^2)
    gL(kk)  =  gL_bar;                      % Conductance of L  (mS/cm^2)
    
    INa(kk)   =  gNa(kk) * (Vm(kk) - ENa);               % Na channel current (uA/cm^2)
    IK(kk)    =  gK(kk) * (Vm(kk) - EK);                 % K channel current (uA/cm^2)
    IL(kk)    =  gL(kk) * (Vm(kk) - EL);                 % L channel current (uA/cm^2)
    Iion(kk)  =  Iinj(kk) - IK(kk) - INa(kk) - IL(kk);   % Total current (uA/cm^2)
    
    % Differential is computed to be used for Euler's method
    dm  =  am * (1 - m(kk)) - bm * m(kk);  
    dn  =  an * (1 - n(kk)) - bn * n(kk);
    dh  =  ah * (1 - h(kk)) - bh * h(kk);
    
    m(kk+1)  =  m(kk) + ts * dm;  % Differential equation of m
    n(kk+1)  =  n(kk) + ts * dn;  % Differential equation of n
    h(kk+1)  =  h(kk) + ts * dh;  % Differential equation of h
    
    % Differential is computed to be used for Euler's method
    dv  =  Iion(kk) / Cm;
    
    Vm(kk+1)  =  Vm(kk) + ts * dv;  % Displacement value of membrane voltage (mV)
    
end
    

%% Plot %%
% Vm_plot is actual membrane voltage. Figure(1) graphs membrane voltage,
% figure(2) conductances of each channel and figure(3) injected current for
% the purpose of checking.


Vm_plot  =  Vm - 70;  % Actual membrane voltage (mV)


figure(1)    
plot(t_vec, Vm_plot)    
title('Membrane Voltage vs. Time Plot')
xlabel('Time (ms)')
ylabel('Vm (mV)')
axis([0, 100, -100, 80])
legend('Vm')

% Conductances and currents are vectors with one length smaller because of
% the iterations in the for-loop. Therefore, g_time_vec was taken from
% t_vec with one length smaller.

g_time_vec  =  t_vec(1 : len - 1);

figure(2)
plot(g_time_vec, gNa, 'b', g_time_vec, gK, '-r', g_time_vec, gL, 'g')
title('Conductances of each channels vs. Time Plot')
xlabel('Time (ms)')
ylabel('g (ms/cm^2)')
legend('gNa','gK','gL')

figure(3)
plot(t_vec, Iinj)
title('Injected Current vs. Time plot')
xlabel('Time (ms)')
ylabel('Current (uA)')
legend('Injected current')







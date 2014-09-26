%% Time %%

ts  =  0.05; % Time-step : 0.05 (ms)

tt  =  100;  % Total simulation time : 100 (ms)
t_vec  =  0 : ts : tt;  % Time vector
len  =  length(t_vec);  % The length of the time vector

%% Injected Current %%

I_inj  =  0 * 10^-3;  % Injected current set to 0 (A/cm^2)

%% Constants %%

gK_bar  =  36 * 10^-3;  % The maximum conductance value of K (S/cm^2)
gNa_bar  =  120 * 10^-3;  % The maximum conductance value of Na (S/cm^2)
gL_bar  =  0.3 * 10^-3;  % The maximum conductance value of L (S/cm^2)

EK  =  -12 * 10^-3;  % Nernst potential of K (V)
ENa  =  115 * 10^-3;  % Nernst potential of Na (V)
EL  =  10.6 * 10^-3;  % Nernst potential of L (V)

Vrest  =  -70 * 10^-3;  % Rest membrane voltage (V)
Cm  =  1.0;  % Membrane capacitance (F/cm^2)

%% Gating Variables %%

Vm(1)  =  Vrest;
am  =  0.1 * ((25 - Vm(1)) / (exp((25 - Vm(1)) / 10) - 1));
bm  =  4 * exp(-Vm(1) / 18);
an  =  0.01 * ((10 - Vm(1)) / (exp((10 - Vm(1)) / 10) - 1));
bn  =  0.125 * exp(-Vm(1) / 80);
ah  =  0.07 * exp(-Vm(1) / 20);
bh  =  1 / (exp((30 - Vm(1)) / 10) + 1);

m(1)  =  am / (am + bm);
n(1)  =  an / (an + bn);
h(1)  =  ah / (ah + bh);


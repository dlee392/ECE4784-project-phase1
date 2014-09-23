%% Membrane Voltage Vm = input %%
V_m=input('What is the membrane voltage?\n');

%% Simulation Parameter %%
%% 1. Total simulation time %%

tt=100*10^-3; % seconds

%% Contasnts Defintion %%
%% 1. Conductance (mS/cm^2) %%
%% g_k, g_Na : f(t), g_L : contant %%

g_k=36;
g_Na=120;
g_l=0.3;

%% 2. Nernst Potential (mV) %%
%% Constants %%

E_k=-12;
E_Na=115;
E_l=10.6;

%% 3. Membrane Rest Potential (mV) %%

V_rest=-70;


%% Equations %%
%% 1. Gating Variables %%

alpha_m=0.1*((25-V_m)/(exp((25-V_m)/10)-1));
beta_m=4*exp(-V_m/18);
alpha_n=.01*((10-V_m)/exp((10-V_m)/10)-1);
beta_n=.125*exp(-V_m/80);
alpha_h=.07*exp(-V_m/20);
beta_h=1/(exp((30-V_m)/10)+1);

%% 2. Derivatives and ODE45 %%

% DVm==Iion/Cm

tr=[0 tt]; % seconds
initial
[t,y]

%% Voltage measurement as displacements from the resting potential %%
% E_r=abs(V_rest);
% V_Na=E_Na-E_r;
% V_k=E_k-E_r;
% V_l=E_l-E_r;
% V=E-E_r; %% What is E? %%


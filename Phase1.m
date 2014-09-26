clear;
clc;

%% Inputs %%
% This section asks for inputs for name, time step and injected current.  
% Name is for a self-padding on the back effect after the program finishes 
% running. Time step is in millisecond and input of the injected current   
% will be divided into magnitude and duration for vectorization purpose    
% later.                                                                   


% name      =  input('What is your name?\n','s');
ts_input  =  input('What is the unit time-step in millisecond?\n'); % time-step
ts        =  ts_input;  % time_step in millisecond

I_inj_comp =  input('What is the magnitude (mA/cm^2) and duration (ms) of the injected current?\ne.g.)[mag dur] and 0 < dur < 100 \n'); % components of I_inj
mag_I_inj  =  I_inj_comp(1) * 10^-3;  % Magnitude of the injected current (A/cm^2)
dur_I_inj  =  I_inj_comp(2);          % Duration of the injected current (ms)

% ts=.598;
% mag_I_inj=0;
% dur_I_inj=0;



%% Simulation time %%
% The simultion time is vectorized with the step size of time step input
% and the length of the vector is taken for the computation later on.


tt       =  100;              % Total simulation time (ms)
t_range  =  0:ts:tt;          % Time range vector
len      =  length(t_range);  % Length of the t_range vector



%% Vectorization of the injected current %%

I_inj  =  zeros(1,len);   % Creates the same size of a vector as t_range

for ll=1:ceil(dur_I_inj/ts)  % Feeds the value into I_inj vector for the duration of the time
    
    I_inj(ll)  =  mag_I_inj;
    ll=ll+1;
end

        
%% Constants Defintion %%

g_K_bar   =  36*10^-3;     % Maximum conductance of K (S/cm^2) 
g_Na_bar  =  120*10^-3;    % Maximum conductance of Na (S/cm^2) 
g_L_bar   =  0.3*10^-3;    % Maximum conductance of L (S/cm^2) 

E_K   =  -12*10^-3;        % Nernst Potential of K (V) 
E_Na  =  115*10^-3;        % Nernst Potential of Na (V) 
E_L   =  10.6*10^-3;       % Nernst Potential of L (V) 

V_rest  =  -70*10^-3;     % Membrane Rest Potential (V) 

C_m  =  1;  % Membrane capacitance (uF/cm^2) 


%% Initial conditions %%

V_m(1)=V_rest; % Vm = Vrest at t=0

m_alpha(1)=  0.1*((25-V_m(1))/(exp((25-V_m(1))/10)-1));   % Gating variable of m
m_beta(1) =  4*exp(-V_m(1)/18);                           % Gating variable of m
n_alpha(1)=  0.01*((10-V_m(1))/(exp((10-V_m(1))/10)-1));  % Gating variable of n
n_beta(1) =  0.125*exp(-V_m(1)/80);                       % Gating variable of n
h_alpha(1)=  0.07*exp(-V_m(1)/20);                        % Gating variable of h
h_beta(1) =  1/(exp((30-V_m(1))/10)+1);                   % Gating variable of h

m(1) =  m_alpha(1)/(m_alpha(1)+m_beta(1)); % m
n(1) =  n_alpha(1)/(n_alpha(1)+n_beta(1)); % n
h(1) =  h_alpha(1)/(h_alpha(1)+h_beta(1)); % h


%% Solving quations using for loop and Euler's method %%

for kk=1:len-1
    
    m_alpha(kk) =  0.1*((25-V_m(kk))/(exp((25-V_m(kk))/10)-1));  % Gating variable of m
    m_beta(kk)  =  4*exp(-V_m(kk)/18);                           % Gating variable of m
    n_alpha(kk) =  0.01*((10-V_m(kk))/(exp((10-V_m(kk))/10)-1)); % Gating variable of n
    n_beta(kk)  =  0.125*exp(-V_m(kk)/80);                       % Gating variable of n
    h_alpha(kk) =  0.07*exp(-V_m(kk)/20);                        % Gating variable of h
    h_beta(kk)  =  1/(exp((30-V_m(kk))/10)+1);                   % Gating variable of h
    
    I_Na(kk)    =  m(kk)^3*g_Na_bar*h(kk)*(V_m(kk)-E_Na);    % Solves for I_Na (A/cm^2)
    I_K(kk)     =  n(kk)^4*g_K_bar*(V_m(kk)-E_K);            % Solves for I_K (A/cm^2)
    I_L(kk)     =  g_L_bar*(V_m(kk)-E_L);                    % Solves for I_L (A/cm^2)
    I_ion(kk)   =  -I_K(kk)-I_Na(kk)-I_L(kk)+I_inj(kk);       % Solves for I_ion (A/cm^2)
    
    m(kk+1)     =  m(kk)+ts*(m_alpha(kk)*(1-m(kk))-m_beta(kk)*m(kk));  % Solves for m 
    n(kk+1)     =  n(kk)+ts*(n_alpha(kk)*(1-n(kk))-n_beta(kk)*n(kk));  % Solves for n 
    h(kk+1)     =  h(kk)+ts*(h_alpha(kk)*(1-h(kk))-h_beta(kk)*h(kk));  % Solves for h 
    
    V_m(kk+1)   =  V_m(kk)+ts*(I_ion(kk)/C_m);         % Solves for membrane voltage
    
    kk=kk+1;
end


%% Plot %%
figure(1)
plot(t_range, I_inj)  % plot of the injected current against time to make sure it's right

figure(2)
plot(t_range, V_m*10^3)    % plot of the membrane voltage against time
axis([0 100 -80 60])
xlabel('Time (ms)')
ylabel('Vm (mV)')

%% Self-padding on the back %%
% fprintf('Good job, %s \n',name)
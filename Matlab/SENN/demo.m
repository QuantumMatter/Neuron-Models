t = 0:1e-3:5;                   % msec
t = t * 1e-3;                   % msec -> sec
pa = 10e-9;                     % A
pd = 1e-3;                      % sec
pw = 0.1e-3;                    % sec
stim = pa * (stepfun(t, pd) - stepfun(t, pd+pw));

V = zeros(1, length(t));        % V    Deviation from resting potential
iNa = zeros(1, length(t));      % A / cm^2
iK = zeros(1, length(t));       % A / cm^2
h = zeros(1, length(t));        % prob, [0,1]
m = zeros(1, length(t));        % prob, [0,1]
n = zeros(1, length(t));        % prob, [0,1]

T = 310.15;     % K

PNa = 51.5e-6;  % m/s      Sodium permeability constant
PK = 2.04e-6;   % m/s      Potassium permeability constant

Na_out = 142.0e-3; % M
Na_in = 10e-3;     % M
K_out = 4.2e-3;    % M
K_in = 141.0e-3;   % M

c_m = 0.02;     % F/m^2     Membrane Capacitance per unit area
dk = 10e-6;     % m        diameter of node k (1.4-2.81)
lk = 10e-6;     % m        length of node k
C_m = c_m * pi * dk * lk;

rhok = 0.7;     % Ohm * m   Axoplasm resistivity
L = 300e-6;     % m        internode distance
Ga = (pi * dk ^ 2) / (4 * rhok * L); % Siemens

gL = 728;       % S/m^2     Leak conductance per unit area
GL = gL * pi * dk * lk; % Siemens

VL = 0;         % V        Leak current equilibrium potential

[alpha_h, beta_h] = dynamics.h_dynamics(V(1), T);
[alpha_m, beta_m] = dynamics.m_dynamics(V(1), T);
[alpha_n, beta_n] = dynamics.n_dynamics(V(1), T);

h(1) = alpha_h / (alpha_h + beta_h);
m(1) = alpha_m / (alpha_m + beta_m);
n(1) = alpha_n / (alpha_n + beta_n);

for i = 2:length(t)

    dt = t(i) - t(i-1);

    [Nai, Ki, hi, mi, ni] = SENNChannelStep(...
        V(i-1), dt, ...
        h(i-1), m(i-1), n(i-1), ...
        PNa, PK, ...
        Na_in, Na_out, ...
        K_in, K_out, ...
        T ...
    );

    % uA / m^2 -> uA / cm^2
    iNa(i) = Nai;
    iK(i) = Ki;
    h(i) = hi;
    m(i) = mi;
    n(i) = ni;

    Iion = (pi * dk * lk) * (iNa(i) + iK(i));
    Ileak = GL * (V(i-1) - VL);
    Istim = stim(i);

    % A / F -> V/sec
    dvdt = (-1/C_m) * (Iion + Ileak - Istim);
    V(i) = V(i-1) + dvdt * dt;

end

figure()
plot(t*1e3, V*1e3); hold on
title('Membrane Potential vs Time')
ylabel('Membrane Potential (mV)')
xlabel('Time (msec)')
yyaxis right
plot(t*1e3, stim)

figure()
plot(t, m); hold on
plot(t, n);
plot(t, h); hold off
legend('m', 'n', 'h')
ylabel('Probability')
xlabel('Time (msec)')
title('Gating Variables vs Time')

figure()
plot(t, iNa * 1e-3); hold on
plot(t, iK * 1e-3); hold off
legend('iNa', 'iK')
ylabel('Current (mA/cm^2)')
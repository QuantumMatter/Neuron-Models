t = 0:1e-3:10;                % msec
% stim = zeros(1, length(t));     % uA
pa = 1e-6;
pd = 1;
pw = 1;
stim = -pa * (stepfun(t, pd) - stepfun(t, pd+pw));

V = zeros(1, length(t));        % mV    Deviation from resting potential
iNa = zeros(1, length(t));      % uA
iK = zeros(1, length(t));       % uA
h = zeros(1, length(t));        % prob, [0,1]
m = zeros(1, length(t));        % prob, [0,1]
n = zeros(1, length(t));        % prob, [0,1]

T = 310.15;     % K

PNa = 51.5;     % um/s      Sodium permeability constant
PK = 2.0;       % um/s      Potassium permeability constant

Na_out = 142.0; % mol/m^3
Na_in = 10;     % mol/m^3
K_out = 4.2;    % mol/m^3
K_in = 141.0;   % mol/m^3

c_m = 0.02;     % F/m^2     Membrane Capacitance per unit area
dk = 2;         % um        diameter of node k (1.4-2.81)
lk = 1;         % um        length of node k
C_m = c_m * pi * (dk * 1e-6) * (lk * 1e-6);
C_m = C_m * 1e6; % uF

rhok = 0.7;     % Ohm * m   Axoplasm resistivity
L = 300;        % um        internode distance
Ga = (pi * (dk * 1e-6) ^ 2) / (4 * rhok * (L * 1e-6)); % Siemens

gL = 728;       % S/m^2     Leak conductance per unit area
GL = gL * pi * (dk * 1e-6) * (lk * 1e-6); % Siemens

VL = 0;         % mV        Leak current equilibrium potential

[alpha_h, beta_h] = dynamics.h_dynamics(V(1), T);
[alpha_m, beta_m] = dynamics.m_dynamics(V(1), T);
[alpha_n, beta_n] = dynamics.n_dynamics(V(1), T);

% V(1) = 0.001;
h(1) = alpha_h / (alpha_h + beta_h);
m(1) = alpha_m / (alpha_m + beta_m);
n(1) = alpha_n / (alpha_n + beta_n);

for i = 2:length(t)

    dt = t(i) - t(i-1);
    dt = dt * 1e-3;

    [Nai, Ki, hi, mi, ni] = SENNChannelStep(...
        V(i-1), dt, ...
        h(i-1), m(i-1), n(i-1), ...
        PNa, PK, ...
        Na_in, Na_out, ...
        K_in, K_out, ...
        T ...
    );

    iNa(i) = Nai;
    iK(i) = Ki;
    h(i) = hi;
    m(i) = mi;
    n(i) = ni;

    
    Iion = pi * (dk * 1e-6) * (lk * 1e-6) * (iNa(i) + iK(i));

    dvdt = (1 / C_m) * (stim(i) - GL * (V(i-1) - VL) - Iion);
    V(i) = V(i-1) + dvdt * dt;

    % disp(dvdt)
    % disp([GL*(V(i-1)-VL) iNa(i) iK(i)])

end

figure()
plot(t, V)

figure()
plot(t, stim)

figure()
plot(t, h); hold on
plot(t, m);
plot(t, n); hold off
legend('H', 'M', 'N')
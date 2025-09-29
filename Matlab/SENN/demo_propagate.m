dt = 1e-6;
t = 0:dt:2e-3; % sec

pa = -2e-3;                     % A
pd = 0.5e-3;                    % sec
pw = 0.04e-3;                   % sec
stim = pa * (stepfun(t, pd) - stepfun(t, pd+pw));

N = 51;

V = zeros(N, length(t));        % V    Deviation from resting potential
iNa = zeros(N, length(t));      % A / cm^2
iK = zeros(N, length(t));       % A / cm^2
h = zeros(N, length(t));        % prob, [0,1]
m = zeros(N, length(t));        % prob, [0,1]
n = zeros(N, length(t));        % prob, [0,1]
Iinj = zeros(N, length(t));     % A

T = 310.15;     % K

PNa = 51.5e-6;  % m/s      Sodium permeability constant
PK = 2.04e-6;   % m/s      Potassium permeability constant

Na_out = 142.0; % M
Na_in = 10;     % M
K_out = 4.2;    % M
K_in = 141.0;   % M

c_m = 0.02;     % F/m^2    Membrane Capacitance per unit area
dk = 10e-6;     % m        diameter of node k (1.4-2.81)
lk = 10e-6;     % m        length of node k
C_m = c_m * pi * dk * lk;

rhoe = 3.5;     % Ohm * m  Extracellular resistivity
rhok = 0.7;     % Ohm * m  Axoplasm resistivity
L = 300e-6;     % m        internode distance
Ga = (pi * dk ^ 2) / (4 * rhok * L); % Siemens

gL = 728;       % S/m^2    Leak conductance per unit area
GL = gL * pi * dk * lk; % Siemens

z = 1e-3;          % m        electrode distance from axon
x = L * (-((N-1)/2):1:(N-1)/2);
r = sqrt(x.*x + z*z);     % distance of each node from the electrode

VL = 0;         % V        Leak current equilibrium potential

[alpha_h, beta_h] = dynamics.h_dynamics(0, T);
[alpha_m, beta_m] = dynamics.m_dynamics(0, T);
[alpha_n, beta_n] = dynamics.n_dynamics(0, T);

h(:,1) = alpha_h / (alpha_h + beta_h);
m(:,1) = alpha_m / (alpha_m + beta_m);
n(:,1) = alpha_n / (alpha_n + beta_n);

for i = 2:length(t)

    [Nai, Ki, hi, mi, ni] = SENNChannelStep(...
        V(:,i-1), -0.0846, dt, ...
        h(:,i-1), m(:,i-1), n(:,i-1), ...
        PNa, PK, ...
        Na_in, Na_out, ...
        K_in, K_out, ...
        T ...
    );

    Nai = (pi * dk * lk) * Nai; % convert to total current
    Ki = (pi * dk * lk) * Ki;   % convert to total current
    iL = GL * (V(:,i-1) - VL);   % leak current; A

    iNa(:,i) = Nai;
    iK(:,i) = Ki;
    h(:,i) = hi;
    m(:,i) = mi;
    n(:,i) = ni;

    % Activating function
    V_e = (rhoe * stim(i)) ./ (4 * pi * r); % extracellular potential; V
    f = second_derivative(V_e); % activating function; V/m^2
    f = f';

    obj = second_derivative(V(:,i-1)); % axial current; A

    % Iinj(:,i) = stim(i);
    Iinj(:,i) = Ga * (f + obj);
    Iion = iNa(:,i) + iK(:,i);

    dvdt = (1/C_m) .* (Iinj(:,i) - Iion - iL);
    V(:,i) = V(:,i-1) + dvdt * dt;
    
end

function [f] = second_derivative(x)
    N = length(x);
    f = zeros(size(x));

    f(1) = x(2) - x(1);
    f(2:N-1) = x(1:N-2) - 2*x(2 : N-1) + x(3:N);
    f(N) = x(N-1) - x(N);
end

% Plot a few adjacent nodes
figure()
for i = 1:5:26
    plot(t*1e3, V(i,:)*1e3, 'DisplayName', sprintf("Node %i", i-1)); hold on
end
legend show
xlabel('Time (msec)')
ylabel('Voltage (mV)')
grid on
yyaxis right
plot(t*1e3, Iinj(26,:)*1e6, 'k', 'DisplayName', 'Stimulus Current at Node 25'); hold off
ylabel('Stimulus Current (uA)')
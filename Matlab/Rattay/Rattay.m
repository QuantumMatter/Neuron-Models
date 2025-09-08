%% Inputs and Outputs
% Inputs
% t - time vector;      msec
% i_elec - stimulus;     uA, same dimensions as t
% z - electrode distance from axon; mm
% N - total number of nodes
% dx - distance between nodes; um
% L - node length; um
% D - fiber diameter; um
% rho_i - intracellular resistivity; kOhm / cm
% rho_e - extracellular resistivity; kOhm / cm
% gNa_bar - max sodium conductance;      mS/cm^2
% gK_bar - max potassium conductance;    mS/cm^2
% gL_bar - max leak conductance;         mS/cm^2
% k_conc_out - external potassium concentration;    mM
% na_conc_out - external sodium concentration;      mM
% V - membrane potential;           mV
% n,m,h - gate probabilities;       %
% INa,IK,IL,Iinj - ion currents;    uA/cm^2

function [V, f, Iinj] = Rattay(...
    t, i_elec, z, ...
    N, L, D, rho_e, rho_i, ...
    gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out)
    

    %% Validate inputs

    
    %% Constants
    c_m = 1;            % uF/cm^2; membrane capacitance per unit area

    na_conc_in = 12;    % mM; intra-cellular sodium concentration
    k_conc_in = 140;    % mM; intra-cellular potassium concentration
    cl_conc_in = 4;     % mM; intra-cellular chloride concentration
    cl_conc_out = 110;  % mM; extra-cellular chloride concentration

    P_Na = 0.03;    % Relative resting permeability of sodium
    P_K = 1;        % Relative resting permeability of potassium
    P_Cl = 0.1;     % Relative resting permeability of chloride

    R = 8.314;      % mJ/K*mmol;    gas constant
    T = 310;        % K;            temperature
    F = 96.4853;    % C/mmol;       Faraday's 

    %% Initialize Variables - Ouputs
    V = zeros(N, length(t));
    f = zeros(N, length(t));
    

    %% Initialize Variables - Local
    n = zeros(N, length(t));
    m = zeros(N, length(t));
    h = zeros(N, length(t));
    INa = zeros(N, length(t));
    IK = zeros(N, length(t));
    IL = zeros(N, length(t));
    Iinj = zeros(N, length(t));

    d = 0.7 * D;    % axon diameter;        um
    dx = 100 * D;   % inter-node spacing;   um
    Ga = (pi * power(d * 1e-4, 2)) / (4 * rho_i * (dx * 1e-4)); % axial conductivity, mS
    C_m = c_m * pi * (d * 1e-4) * (L * 1e-4);                     % nodal membrane capacitance; uF

    x = dx * (-(N-1)/2 : 1 : (N-1)/2);      % x coordinates of nodes; um
    % r = sqrt(x.*x + power(z * 1e3, 2));     % distance of each node from the electrode

    %% Calculate Initial Values

    % Resting potential
    V_rest = (R*T/F) * log( ...
        (P_Na*na_conc_out + P_K*k_conc_out + P_Cl*cl_conc_in ) ...
        / ...
        (P_Na*na_conc_in  + P_K*k_conc_in  + P_Cl*cl_conc_out) ...
    );
    V(:,1) = V_rest;

    % Nernst potentials
    E_Na = (R*T/F) * log(na_conc_out / na_conc_in);
    E_K = (R*T/F) * log(k_conc_out / k_conc_in);

    % Channel conductance probabilities
    n(:,1) = calc_n(V_rest);
    m(:,1) = calc_m(V_rest);
    h(:,1) = calc_h(V_rest);

    % Channel conductances
    gNa_0 = calc_gNa(gNa_bar, m(1,1), h(1,1));
    gK_0 = calc_gK(gK_bar, n(1,1));

    % Leak potential
    E_L = V_rest + ( ...
        (gK_0*(V_rest - E_K) + gNa_0*(V_rest - E_Na)) ...
        / ...
        gL_bar ...
    );
    % fprintf('E_K: %f, E_Na: %f, E_L: %f\r\n', E_K, E_Na, E_L);
    
    % Ion currents
    INa(:,1) = gNa_0 * (V_rest - E_Na);
    IK(:,1) = gK_0 * (V_rest - E_K);
    IL(:,1) = gL_bar * (V_rest - E_L);

    %% Simulate the ODEs with Euler's method
    for i = 2:length(t)

        dt = t(i) - t(i-1);

        [mi, hi, ni, ina_step, ik_step, il_step] = HH_step(dt, V(:,i-1), m(:,i-1), h(:,i-1), n(:,i-1), gNa_bar, gK_bar, gL_bar, E_Na, E_K, E_L);

        m(:,i) = mi;
        n(:,i) = ni;
        h(:,i) = hi;

        INa(:,i) = ina_step;
        IK(:,i) = ik_step;
        IL(:,i) = il_step;

        % uA / cm^2
        i_ionic = INa(:,i) + IK(:,i) + IL(:,i);
        I_ionic = pi * (d * 1e-4) * (L * 1e-4) * i_ionic;

        % V_e = (rho_e * i_elec(i)) ./ (4 * pi * (r * 1e-4));
        % f(1,i) = V_e(2) - V_e(1);
        % f(2:N-1 , i) = V_e(1:N-2) - 2*V_e(2 : N-1) + V_e(3:N);
        % f(N,i) = V_e(N-1) - V_e(N);

        % mV / msec
        [fi] = activation_func(x, z, i_elec(i), rho_e, rho_i, d, L);
        f(:,i) = fi';

        % obj = zeros(N);
        % obj(1) = V(2, i-1) - V(1, i-1);
        % obj(2:N-1) = V(1:N-2, i-1) - 2*V(2:N-1, i-1) + V(3:N, i-1);
        % obj(N) = V(N-1,i-1) - V(N,i-1);
        [obj] = second_derivative(V(:,i-1));

        % Iinj(:,i) = f(:,i) .* C_m;
        % 
        % dvdt = f(:,i) + (Ga * obj - I_ionic) / C_m;
        % 
        % V(:,i) = V(:,i-1) + dvdt .* dt;

        % cap: i = C * dvdt -> dvdt = i / C
        % mV/msec + uA / uF
        % dvdt = f(:,i); % + I_ionic / C_m;
        dvdt = Ga * (f(:,i) + obj) / C_m;
        % disp(dvdt)
        V(:,i) = V(:,i-1) + dvdt * dt;

    end
end

function [f] = second_derivative(x)
    N = length(x);
    f = zeros(size(x));

    f(1) = x(2) - x(1);
    f(2:N-1) = x(1:N-2) - 2*x(2 : N-1) + x(3:N);
    f(N) = x(N-1) - x(N);
end

% Activation Function - Descibes the effect of the electrode on the nodes
% x - The x-coordinates of each of the nodes;           um
% z - The distance from the axon to the electrode;      mm
% i_elec - The electrode stimulus;                      uA
% rho_e - The resistivity of the extra-cellular space;  kOhm * cm
% rho_i - The resistivity of the intra-cellular space;  kOhm * cm
% d - The diameter of the axon;                         um
% L - Length of each node;                              um
% f - The value of the activation function at each x;   mV / msec
function [f] = activation_func(x, z, i_elec, rho_e, rho_i, d, L)
    
    % Assume that all x is equally spaced; um
    dx = x(2) - x(1);
    
    % axial conductivity, mS
    % cm^2 / ((kOhm * cm) * cm) = 1/kOhm
    Ga = (pi * power(d * 1e-4, 2)) / (4 * rho_i * (dx * 1e-4)); 

    % nodal membrane capacitance; uF
    % uF/cm^2 * cm  * cm
    C_m = 1  * pi * (d * 1e-4) * (L * 1e-4);

    % distance of each node from the electrode; um
    r = sqrt(x.*x + power(z * 1e3, 2));

    % mV = (kOhm * cm * uA) / (cm)
    V_e = (rho_e * i_elec) ./ (4 * pi * (r * 1e-4));

    % mV
    [f] = second_derivative(V_e);

    % mS / uF = 1 / (uF * kOhm) = 1 / msec
    % f = Ga / C_m * f;
end

function [mi, hi, ni, I_Na, I_K, I_L] = HH_step(dt, V, m, h, n, gNa_bar, gK_bar, gL_bar, E_Na, E_K, E_L)

    % I_Na = zeros(size(v));
    % I_K = zeros(size(V));
    % I_L = zeros(size(V));

    [alpha_m, beta_m] = calc_m_rates(V);
    [alpha_h, beta_h] = calc_h_rates(V);
    [alpha_n, beta_n] = calc_n_rates(V);

    dmdt = alpha_m .* (1 - m) - beta_m .* m;
    dhdt = alpha_h .* (1 - h) - beta_h .* h;
    dndt = alpha_n .* (1 - n) - beta_n .* n;

    mi = m + dt * dmdt;
    hi = h + dt * dhdt;
    ni = n + dt * dndt;

    gNa = calc_gNa(gNa_bar, m, h);
    gK = calc_gK(gK_bar, n);

    I_Na = gNa .* (V - E_Na);
    I_K = gK .* (V - E_K);
    I_L = gL_bar .* (V - E_L);

end

% Calculate the forward and reverse rate constants for the 'm' gate
% (Na+ activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_m, beta_m] = calc_m_rates(V)
    alpha_m = (0.1 * (V + 40)) ./ (1 - exp(-(V+40)/10));
    beta_m = 4 * exp(-(V+65)/18);
end

% Calculate the forward and reverse rate constants for the 'h' gate
% (Na+ de-activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_h, beta_h] = calc_h_rates(V)
    alpha_h = 0.07 * exp(-(V+65)/20);
    beta_h = 1 ./ (1 + exp(-(V+35)/10));
end

% Calculate the forward and reverse rate constants for the 'n' gate
% (K+ activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_n, beta_n] = calc_n_rates(V)
    alpha_n = (0.01 * (V + 55)) ./ (1 - exp(-(V+55)/10));
    beta_n = 0.125 * exp(-(V + 65) / 80);
end

function [m] = calc_m(V)
    [alpha_m, beta_m] = calc_m_rates(V);
    m = alpha_m ./ (alpha_m + beta_m);
end

function [h] = calc_h(V)
    [alpha_h, beta_h] = calc_h_rates(V);
    h = alpha_h ./ (alpha_h + beta_h);
end

function [n] = calc_n(V)
    [alpha_n, beta_n] = calc_n_rates(V);
    n = alpha_n ./ (alpha_n + beta_n);
end

function [gNa] = calc_gNa(gNa_bar, m, h)
    gNa = gNa_bar * power(m, 3) .* h;
end

function [gK] = calc_gK(gK_bar, n)
    gK = gK_bar * power(n, 4);
end
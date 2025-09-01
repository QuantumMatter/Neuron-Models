% t - time vector;          usec
% pd - pulse delay;         usec
% pa - pulse amplitude;     uA/cm^2
% pw - pulse width;         usec

% gNabar - max sodium conductance;      mS/cm^2
% gKbar - max potassium conductance;    mS/cm^2
% gLbar - max leak conductance;         mS/cm^2

% k_conc_out - external potassium concentration;    mM
% na_conc_out - external sodium concentration;      mM

% V - membrane potential;           mV
% n,m,h - gate probabilities;       %
% INa,IK,IL,Iinj - ion currents;    uA/cm^2

function [V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa, pw, gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out)

    %% Constants
    c_m = 1;            % uF/cm^2; membrane capacitance per unit area
    
    na_conc_in = 12;    % mM; intra-cellular sodium concentration
    k_conc_in = 140;    % mM; intra-cellular potassium concentration
    cl_conc_in = 4;     % mM; intra-cellular chloride concentration
    % ??? Should this be an argument?
    cl_conc_out = 110;  % mM; extra-cellular chloride concentration
    
    P_Na = 0.03;    % Relative resting permeability of sodium
    P_K = 1;        % Relative resting permeability of potassium
    P_Cl = 0.1;     % Relative resting permeability of chloride
    
    R = 8.314;      % mJ/K*mmol;    gas constant
    T = 310;        % K;            temperature
    F = 96.4853;    % C/mmol;       Faraday's constant

    %% Initialize variables
    V = zeros(size(t));
    n = zeros(size(t));
    m = zeros(size(t));
    h = zeros(size(t));
    INa = zeros(size(t));
    IK = zeros(size(t));
    IL = zeros(size(t));
    Iinj = zeros(size(t));

    %% Start by finding the initial values

    % Resting potential
    V_rest = (R*T/F) * log( ...
        (P_Na*na_conc_out + P_K*k_conc_out + P_Cl*cl_conc_in ) ...
        / ...
        (P_Na*na_conc_in  + P_K*k_conc_in  + P_Cl*cl_conc_out) ...
    );
    V(1) = V_rest;

    % Nernst potentials
    E_Na = (R*T/F) * log(na_conc_out / na_conc_in);
    E_K = (R*T/F) * log(k_conc_out / k_conc_in);

    % Channel conductance probabilities
    n(1) = calc_n(V_rest);
    m(1) = calc_m(V_rest);
    h(1) = calc_h(V_rest);

    % Channel conductances
    gNa_0 = calc_gNa(gNa_bar, m(1), h(1));
    gK_0 = calc_gK(gK_bar, n(1));

    % Leak potential
    E_L = V_rest + ( ...
        (gK_0*(V_rest - E_K) + gNa_0*(V_rest - E_Na)) ...
        / ...
        gL_bar ...
    );
    % fprintf('E_K: %f, E_Na: %f, E_L: %f\r\n', E_K, E_Na, E_L);
    
    % Ion currents
    INa(1) = gNa_0 * (V_rest - E_Na);
    IK(1) = gK_0 * (V_rest - E_K);
    IL(1) = gL_bar * (V_rest - E_L);

    %% Simulate the ODEs with Euler's Method
    for i = 2:length(t)

        dt = t(i) - t(i-1); % Should be 1 usec

        [alpha_m, beta_m] = calc_m_rates(V(i-1));
        [alpha_h, beta_h] = calc_h_rates(V(i-1));
        [alpha_n, beta_n] = calc_n_rates(V(i-1));

        dmdt = alpha_m * (1 - m(i-1)) - beta_m*m(i-1);
        dhdt = alpha_h * (1 - h(i-1)) - beta_h*h(i-1);
        dndt = alpha_n * (1 - n(i-1)) - beta_n*n(i-1);

        m(i) = m(i-1) + dt * dmdt;
        h(i) = h(i-1) + dt * dhdt;
        n(i) = n(i-1) + dt * dndt;

        % Calculate the currents of each ion channel
        gNa = calc_gNa(gNa_bar, m(i-1), h(i-1));
        gK = calc_gK(gK_bar, n(i-1));

        INa(i) = gNa * (V(i-1) - E_Na);
        IK(i) = gK * (V(i-1) - E_K);
        IL(i) = gL_bar * (V(i-1) - E_L);

        % Calculate the difference between the current stimulus
        % and the currents of the channels
        % I_stim = I_cm + I_L + I_K + I_Na
        if (t(i) > pd) && (t(i) < (pd + pw))
            Iinj(i) = pa;
        end

        I_cm = Iinj(i) - INa(i) - IK(i) - IL(i);

        % Calculate the new membrane potential
        % I_cm = C_M * dvdt
        dvdt = I_cm / c_m; % uA / uA/cm^2
        V(i) = V(i-1) + dt * dvdt;

    end

end

% Calculate the forward and reverse rate constants for the 'm' gate
% (Na+ activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_m, beta_m] = calc_m_rates(V)
    alpha_m = (0.1 * (V + 40)) / (1 - exp(-(V+40)/10));
    beta_m = 4 * exp(-(V+65)/18);
end

% Calculate the forward and reverse rate constants for the 'h' gate
% (Na+ de-activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_h, beta_h] = calc_h_rates(V)
    alpha_h = 0.07 * exp(-(V+65)/20);
    beta_h = 1 / (1 + exp(-(V+35)/10));
end

% Calculate the forward and reverse rate constants for the 'n' gate
% (K+ activation gate)
% Units are ms^-1
% V - membrane potential in millivolts
function [alpha_n, beta_n] = calc_n_rates(V)
    alpha_n = (0.01 * (V + 55)) / (1 - exp(-(V+55)/10));
    beta_n = 0.125 * exp(-(V + 65) / 80);
end

function [m] = calc_m(V)
    [alpha_m, beta_m] = calc_m_rates(V);
    m = alpha_m / (alpha_m + beta_m);
end

function [h] = calc_h(V)
    [alpha_h, beta_h] = calc_h_rates(V);
    h = alpha_h / (alpha_h + beta_h);
end

function [n] = calc_n(V)
    [alpha_n, beta_n] = calc_n_rates(V);
    n = alpha_n / (alpha_n + beta_n);
end

function [gNa] = calc_gNa(gNa_bar, m, h)
    gNa = gNa_bar * power(m, 3) * h;
end

function [gK] = calc_gK(gK_bar, n)
    gK = gK_bar * power(n, 4);
end



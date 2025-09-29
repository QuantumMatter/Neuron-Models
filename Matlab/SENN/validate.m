% Compares my Matlab implementation to Evan's Java implementation
% Specifically, it compares the following
%   - Channel dynamics, forward and backward rates vs Vm
%       - (a_m, b_m, ...) vs (Vm, T)
%   - Current for a given state
%       - (i_Na, i_K, i_L) vs (Vm, m, h, n)
%   - Response to constant stimulus

clear java;
javaaddpath('.\+ref\');

%%

function validate_dynamics()

    V = -100:0.1:100; 
    V = V / 1000; % mV to V
    T = GlobalConstants.T;

    % Preallocate
    david = zeros(length(V), 6);
    evan = zeros(length(V), 6);
    
    for i = 1:length(V)
        v = V(i);

        % David's implementation
        [a_m, b_m] = dynamics.m_dynamics(v * 1000, T);
        [a_h, b_h] = dynamics.h_dynamics(v * 1000, T);
        [a_n, b_n] = dynamics.n_dynamics(v * 1000, T);

        david(i, :) = [a_m, b_m, a_h, b_h, a_n, b_n];

        % Evan's implementation
        a_n = F_PotassiumChannel.alpha(v + GlobalConstants.Vrest);
        b_n = F_PotassiumChannel.beta(v + GlobalConstants.Vrest);
        a_h = F_SodiumChannel.alpha_h(v + GlobalConstants.Vrest);
        b_h = F_SodiumChannel.beta_h(v + GlobalConstants.Vrest);
        a_m = F_SodiumChannel.alpha_m(v + GlobalConstants.Vrest);
        b_m = F_SodiumChannel.beta_m(v + GlobalConstants.Vrest);
        evan(i, :) = [a_m, b_m, a_h, b_h, a_n, b_n];
    end

    % Calculate the differences
    error = sum(abs(david - evan));
    fprintf('Error: %f\n', error);

    % Plot the results as a surface
    surf(david - evan);

    % print out the V values where there are NaN
    for i = 1:length(V)
        if any(isnan(david(i, :))) || any(isnan(evan(i, :)))
            fprintf('NaN at V = %f mV\n', V(i) * 1000);
        end
    end

end

%%

function [V, david, evan] = validate_potassium_current()

    V = -100:0.1:100; 
    V = V / 1000; % mV to V
    T = GlobalConstants.T;

    n = 0:0.01:1;

    % Add a small offset to avoid numerical errors
    V = V + 1e-6;
    n = n + 1e-6;
    n(n > 1) = 1;

    % Constants
    PNa = F_SodiumChannel.Pna;
    PK = F_PotassiumChannel.Pk;
    Na_in = GlobalConstants.Conc_NaI;
    Na_out = GlobalConstants.Conc_NaO;
    K_in = GlobalConstants.Conc_KI;
    K_out = GlobalConstants.Conc_KO;

    nd = 10e-6; % Diameter in m
    nl = 10e-6; % Length in m

    % Preallocate
    % david_sodium = zeros(length(V), length(m), length(h));
    david = zeros(length(V), length(n));
    % evan_sodium = zeros(length(V), length(m), length(h));
    evan = zeros(length(V), length(n));


    potassiumChannel = F_PotassiumChannel(nd, nl, 0);

    
    % Use steady state values for gating variables
    for i = 1:length(V)

        v = V(i);

        for j = 1:length(n)

            ni = n(j);
            mi = 0.5; % dummy
            hi = 0.5; % dummy

            % David's implementation
            [~, iK, ~, ~, ~] = SENNChannelStep(...
                v, -0.0846, 0, ...           % V, sec
                hi, mi, ni, ...     % h, m, n
                PNa, PK, ...        % PNa, PK
                Na_in, Na_out, ...  % Na_in, Na_out
                K_in, K_out, ...    % K_in, K_out
                T ...               % T
            );
            iK = (pi * nd * nl) * iK; % convert to total current
            david(i, j) = iK;

            % Evan's implementation
            potassiumChannel.setState(ni);
            iK = potassiumChannel.current(v + GlobalConstants.Vrest);
            evan(i, j) = iK;

        end
    end

    % Calculate the differences
    error_k = sum(abs(david - evan), 'all');
    fprintf('Potassium Error: %f\n', error_k);

end

function [V, david, evan] = validate_sodium_current()

    V = -100:0.5:100; 
    V = V / 1000; % mV to V
    T = GlobalConstants.T;

    m = 0:0.01:1;
    h = 0:0.01:1;

    % Add a small offset to avoid numerical errors
    V = V + 1e-6;
    m = m + 1e-6;
    h = h + 1e-6;
    m(m > 1) = 1;
    h(h > 1) = 1;

    % Constants
    PNa = F_SodiumChannel.Pna;
    PK = F_PotassiumChannel.Pk;
    Na_in = GlobalConstants.Conc_NaI;
    Na_out = GlobalConstants.Conc_NaO;
    K_in = GlobalConstants.Conc_KI;
    K_out = GlobalConstants.Conc_KO;

    nd = 10e-6; % Diameter in m
    nl = 10e-6; % Length in m

    % Preallocate
    david = zeros(length(V), length(m), length(h));
    evan = zeros(length(V), length(m), length(h));

    
    sodiumChannel = F_SodiumChannel(nd, nl, 0);

    
    % Use steady state values for gating variables
    for i = 1:length(V)
        fprintf('V = %f mV\n', V(i) * 1000);

        v = V(i);

        for j = 1:length(m)
            mi = m(j);

            for k = 1:length(h)
                hi = h(k);

                ni = 0.5; % dummy

                % David's implementation
                [iNa, ~, ~, ~, ~] = SENNChannelStep(...
                    v, -0.0846, 0, ...           % V, sec
                    hi, mi, ni, ...     % h, m, n
                    PNa, PK, ...        % PNa, PK
                    Na_in, Na_out, ...  % Na_in, Na_out
                    K_in, K_out, ...    % K_in, K_out
                    T ...               % T
                );
                iNa = (pi * nd * nl) * iNa; % convert to total current
                david(i, j, k) = iNa;

                % Evan's implementation
                sodiumChannel.setState([mi, hi]);
                iNa = sodiumChannel.current(v + GlobalConstants.Vrest);
                evan(i, j, k) = iNa;
            end
        end
    end

    % Calculate the differences
    error_na = sum(abs(david - evan), 'all');
    fprintf('Sodium Error: %f\n', error_na);
    
end

%%

function [t, V, stim, h, m, n, iNa, iK] = validate_eulers()

    % Evaluate one trajectory of the two implementations and 
    % identify if and when they diverge

    dt = 1e-6;
    t = 0:dt:0.002; % sec
    % t = 0:dt:0.0004; % sec
    V = zeros(2, length(t));
    h = zeros(2, length(t));
    m = zeros(2, length(t));
    n = zeros(2, length(t));
    iNa = zeros(2, length(t));
    iK = zeros(2, length(t));

    pa = 10e-9;                     % A
    pd = 0.5e-3;                      % sec
    pw = 0.04e-3;                    % sec
    stim = pa * (stepfun(t, pd) - stepfun(t, pd+pw));

    % Constants
    T = GlobalConstants.T;
    PNa = F_SodiumChannel.Pna;
    PK = F_PotassiumChannel.Pk;
    Na_in = GlobalConstants.Conc_NaI;
    Na_out = GlobalConstants.Conc_NaO;
    K_in = GlobalConstants.Conc_KI;
    K_out = GlobalConstants.Conc_KO;

    nd = 10e-6; % Diameter in m
    nl = 10e-6; % Length in m

    c_m = 0.02;                     % F/m^2     Membrane Capacitance per unit area
    C_m = c_m * pi * nd * nl;       % Farads

    gL = 728;       % S/m^2     Leak conductance per unit area
    GL = gL * pi * nd * nl; % Siemens

    VL = 0;         % V        Leak current equilibrium potential

    [alpha_h, beta_h] = dynamics.h_dynamics(V(1), T);
    [alpha_m, beta_m] = dynamics.m_dynamics(V(1), T);
    [alpha_n, beta_n] = dynamics.n_dynamics(V(1), T);

    h(:,1) = alpha_h / (alpha_h + beta_h);
    m(:,1) = alpha_m / (alpha_m + beta_m);
    n(:,1) = alpha_n / (alpha_n + beta_n);

    node = SENN_ActiveNode(nd, nl, dt);
    node.channel_Na.setState([m(2,1), h(2,1)]);
    node.channel_K.setState(n(2,1));

    for i = 2:length(t)

        % David's implementation
        [Nai, Ki, hi, mi, ni] = SENNChannelStep(...
            V(1,i-1), -0.0846, dt, ...
            h(1,i-1), m(1,i-1), n(1,i-1), ...
            PNa, PK, ...
            Na_in, Na_out, ...
            K_in, K_out, ...
            T ...
        );

        Nai = (pi * nd * nl) * Nai; % convert to total current
        Ki = (pi * nd * nl) * Ki;   % convert to total current

        iNa(1,i) = Nai;
        iK(1,i) = Ki;

        h(1,i) = hi;
        m(1,i) = mi;
        n(1,i) = ni;

        dVdt = (stim(i) - Nai - Ki - GL * (V(1,i-1) - VL)) / C_m;
        V(1,i) = V(1,i-1) + dVdt * dt;

        % Evan's implementation
        evan = node.compute(0, 0, 0, V(2,i-1), 0, 0, 1, 1, 0, 0, stim(i));
        V(2,i) = evan(1);
        iNa(2,i) = evan(2);
        iK(2,i) = evan(3);

        n(2,i) = node.channel_K.getState();
        state = node.channel_Na.getState();
        m(2,i) = state(1);
        h(2,i) = state(2);

    end

end

%%

% This finds issues, but they're only numerical stability issues
% with the Java implementation. The Matlab implementation is correct.
% Its unlikely that we'll encounter these issues in practice, because the
% hole in the curve is very small.
%   [x,y] = dynamics.n_dynamics(35, GlobalConstants.T)
%   F_PotassiumChannel.beta(0.035 + GlobalConstants.Vrest)
% validate_dynamics();

% [V, david, evan] = validate_potassium_current();
% [V, david, evan] = validate_sodium_current();
% error = abs(david - evan);
% v_error = max(error, [], [2 3]);

[t, V, stim, h, m, n, iNa, iK] = validate_eulers();
% Plot the results
figure()
plot(t*1e3, V(1,:)*1e3); hold on
plot(t*1e3, V(2,:)*1e3);
title('Membrane Potential vs Time')
ylabel('Membrane Potential (mV)')
xlabel('Time (msec)')
% yyaxis right
% plot(t*1e3, stim*1e6)
% ylabel('Stimulus (uA)')
% legend('David', 'Evan', 'Stimulus')
% zoom yon

% Make a figure comparing each of the gating variables
figure()
subplot(3,1,1)
plot(t*1e3, m(1,:)); hold on
plot(t*1e3, m(2,:));
title('m gating variable vs Time')
ylabel('m')
xlabel('Time (msec)')
legend('David', 'Evan')
subplot(3,1,2)
plot(t*1e3, h(1,:)); hold on
plot(t*1e3, h(2,:));    
title('h gating variable vs Time')
ylabel('h')
xlabel('Time (msec)')
legend('David', 'Evan')
subplot(3,1,3)
plot(t*1e3, n(1,:)); hold on
plot(t*1e3, n(2,:));
title('n gating variable vs Time')
ylabel('n')
xlabel('Time (msec)')
legend('David', 'Evan')

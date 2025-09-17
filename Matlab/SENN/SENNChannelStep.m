% Simulates the ion channel/current dynamics
% at a single node of Ravnier, and does not
% consider the effect of the activating function

% Expect T = 310.15 K
function [iNa, iK, h, m, n] = SENNChannelStep(...
    V, dt, ...
    h_prev, m_prev, n_prev, ...
    PNa, PK, ...
    Na_in, Na_out, ...
    K_in, K_out, ...
    T ...
)

    F = 96485;  % Coulombs / mol
    R = 8.314;  % J / (K * mol)

    % In terms of msec, mV, and K
    [ah, bh] = dynamics.h_dynamics(V, T);
    [am, bm] = dynamics.m_dynamics(V, T);
    [an, bn] = dynamics.n_dynamics(V, T);

    dhdt = ah - (ah + bh) * h_prev;
    dmdt = am - (am + bm) * m_prev;
    dndt = an - (an + bn) * n_prev;

    % disp([dhdt dmdt dndt])

    h = h_prev + dhdt * dt;
    m = m_prev + dmdt * dt;
    n = n_prev + dndt * dt;
    
    V = V * 1e3;    % mV -> V

    % (um/s) * (C^2 / mol^2)
    iNa = dynamics.vtrap(V, F/(R*T));
    iNa = iNa * PNa * h * m^3 * F^2 * (Na_out - Na_in * exp((V*F)/(R*T)));
    iNa = iNa / (R*T);

    iK = dynamics.vtrap(V, F/(R*T));
    iK = iK * PK * n^2 * F^2 * (K_out - K_in * exp((V*F)/(R*T)));
    iK = iK / (R*T);

    % disp([iNa iK])

end

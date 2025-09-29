% Simulates the ion channel/current dynamics
% at a single node of Ravnier, and does not
% consider the effect of the activating function

% Expect T = 310.15 K
function [iNa, iK, h, m, n] = SENNChannelStep(...
    V, Vrest, dt, ...
    h_prev, m_prev, n_prev, ...
    PNa, PK, ...
    Na_in, Na_out, ...
    K_in, K_out, ...
    T ...
)

    F = 96485;  % Coulombs / mol
    R = 8.314;  % J / (K * mol)

    % In terms of msec, mV, and K, gives mV/msec
    [ah, bh] = dynamics.h_dynamics(V * 1e3, T);
    [am, bm] = dynamics.m_dynamics(V * 1e3, T);
    [an, bn] = dynamics.n_dynamics(V * 1e3, T);

    dhdt = ah - (ah + bh) * h_prev;
    dmdt = am - (am + bm) * m_prev;
    dndt = an - (an + bn) * n_prev;

    % dhdt = ah * (1 - h_prev) - bh * h_prev;
    % dmdt = am * (1 - m_prev) - bm * m_prev;
    % dndt = an * (1 - n_prev) - bn * n_prev;

    h = h_prev + dhdt * dt * 1e3;
    m = m_prev + dmdt * dt * 1e3;
    n = n_prev + dndt * dt * 1e3;

    V = V + Vrest;

    % (um/s) * (C^2 / mol^2)
    iNa = dynamics.vtrap(V, F/(R*T));
    iNa = iNa * PNa * h * m^3 * F^2 * (Na_out - Na_in * exp((V*F)/(R*T)));
    iNa = iNa / (R*T);

    iK = dynamics.vtrap(V, F/(R*T));
    iK = iK * PK * n^2 * F^2 * (K_out - K_in * exp((V*F)/(R*T)));
    iK = iK / (R*T);

end

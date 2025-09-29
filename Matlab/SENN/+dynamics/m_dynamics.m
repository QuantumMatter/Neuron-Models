function [alpha, beta] = m_dynamics(V, T)

    % V is mV
    % K is Kelvin

    Aam = 0.49;     % msec^-1
    Abm = 1.04;     % msec^-1
    Bam = 25.41;    % mV
    Bbm = 21.00;    % mV
    Cam = 6.06;     % mV
    Cbm = 9.41;     % mV
    Qam = 2.2;
    Qbm = 2.2;
    T0 = 293.15;

    alpha = -Aam * power(Qam, (T-T0)/10);
    alpha = alpha .* dynamics.vtrap(Bam - V, 1/Cam);

    beta = Abm * power(Qbm, (T-T0)/10);
    beta = beta .* dynamics.vtrap(Bbm - V, -1/Cbm);

end
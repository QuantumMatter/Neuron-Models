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
    T0 = 37 + 273.15;

    alpha = power(Qam, (T-T0)/10) * ...
        (Aam * (V - Bam)) ...
        / ...
        (1 - exp((Bam - V) / Cam));

    beta = power(Qbm, (T-T0)/10) * ...
        (Abm * (Bbm - V)) ...
        / ...
        (1 - exp((V - Bbm) / Cbm));

end
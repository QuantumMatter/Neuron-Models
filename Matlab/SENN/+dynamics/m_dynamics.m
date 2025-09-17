function [alpha, beta] = m_dynamics(V, T)

    Aam = 0.49;
    Abm = 1.04;
    Bam = 25.41;
    Bbm = 21.00;
    Cam = 6.06;
    Cbm = 9.41;
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
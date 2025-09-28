function [alpha, beta] = h_dynamics(V, T)
    
    Aah = 0.09; % ms^-1
    Abh = 3.70;
    Bah = -27.74;
    Bbh = 56.00;
    Cah = 9.06;
    Cbh = 12.50;
    Qah = 2.9;
    Qbh = 2.9;
    T0 = 293.15;

    alpha = Aah * power(Qah, (T-T0)/10);
    alpha = alpha * dynamics.vtrap(Bah - V, -1/Cah);

    beta = power(Qbh, (T-T0)/10) * ...
        (Abh) ...
        / ...
        (1 + exp((Bbh - V) / Cbh));

end
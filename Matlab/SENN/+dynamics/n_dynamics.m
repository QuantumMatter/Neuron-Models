function [alpha, beta] = n_dynamics(V, T)

    Aan = 0.02;
    Abn = 0.05;
    Ban = 35.00;
    Bbn = 10.00;
    Can = 10.00;
    Cbn = 10.00;
    Qan = 3.0;
    Qbn = 3.0;
    T0 = 293.15;

    alpha = power(Qan, (T-T0)/10) * ...
        (Aan * (V - Ban)) ...
        / ...
        (1 - exp((Ban - V) / Can));

    beta = power(Qbn, (T-T0)/10) * ...
        (Abn * (Bbn - V)) ...
        / ...
        (1 - exp((V - Bbn) / Cbn));

end
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

    alpha = -Aan * power(Qan, (T-T0)/10);
    alpha = alpha .* dynamics.vtrap(Ban - V, 1/Can);

    beta = Abn * power(Qbn, (T-T0)/10);
    beta = beta .* dynamics.vtrap(Bbn - V, -1/Cbn);

end
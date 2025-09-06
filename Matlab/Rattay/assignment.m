%% Part 1
% In one figure, plot the activation function for a range of electrode currents

function [x, f] = activation_func(i_elec, N, rho_e, rho_i, z, D, L)
    dx = L / (N-1); % um
    
    d = 0.7 * D;    % axon diameter;        um
    Ga = (pi * power(d * 1e-4, 2)) / (4 * rho_i * (dx * 1e-4)); % axial conductivity, mS
    % uF/cm^2 * cm  * cm
    C_m = 1  * pi * (d * 1e-4) * (2.5 * 1e-4);                     % nodal membrane capacitance; uF

    x = dx * (-(N-1)/2 : 1 : (N-1)/2);      % x coordinates of nodes; um
    r = sqrt(x.*x + power(z * 1e3, 2));     % distance of each node from the electrode; um

    f = zeros(size(x));

    % mV = (kOhm * cm * uA) / (cm)
    V_e = (rho_e * i_elec) ./ (4 * pi * (r * 1e-4));
    f(1) = V_e(2) - V_e(1);
    f(2:N-1) = V_e(1:N-2) - 2*V_e(2 : N-1) + V_e(3:N);
    f(N) = V_e(N-1) - V_e(N);

    % mS / uF = 1 / (uF * kOhm) = ms^-1
    f = Ga / C_m * f;
end

figure()
for i_elec = [-500 500 -2500 2500]
% for i_elec = [2500]
    [x, f] = activation_func(i_elec, 51, 0.3, 0.055, 1, 10, 100*10*51);
    plot(x, f); hold on
end
hold off
% ylim([-100 100])

% R = rho * l / A
%   = rho * l / (pi * r ^ 2)
%   = rho * l / (pi * d ^ 2 / 4)
%   = 4 * rho * l / (pi * d ^ 2)
%   = (kOhm * cm) * cm / cm^2
%   = kOhm
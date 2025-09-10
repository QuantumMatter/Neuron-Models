%% Part 1
% In one figure, plot the activation function for a range of electrode currents

function [f] = second_derivative(x)
    N = length(x);
    f = zeros(size(x));

    f(1) = x(2) - x(1);
    f(2:N-1) = x(1:N-2) - 2.*x(2 : N-1) + x(3:N);
    f(N) = x(N-1) - x(N);
end

% Activation Function - Descibes the effect of the electrode on the nodes
% x - The x-coordinates of each of the nodes;           um
% z - The distance from the axon to the electrode;      mm
% i_elec - The electrode stimulus;                      uA
% rho_e - The resistivity of the extra-cellular space;  kOhm * cm
% rho_i - The resistivity of the intra-cellular space;  kOhm * cm
% d - The diameter of the axon;                         um
% f - The value of the activation function at each x;   mV / msec
function [f] = activation_func(x, z, i_elec, rho_e, rho_i, d)
    
    % Assume that all x is equally spaced
    dx = x(2) - x(1);
    
    Ga = (pi * power(d * 1e-4, 2)) / (4 * rho_i * (dx * 1e-4)); % axial conductivity, mS
    % uF/cm^2 * cm  * cm
    C_m = 1  * pi * (d * 1e-4) * (2.5 * 1e-4);                     % nodal membrane capacitance; uF

    r = sqrt(x.*x + power(z * 1e3, 2));     % distance of each node from the electrode; um

    % mV = (kOhm * cm * uA) / (cm)
    V_e = (rho_e * i_elec) ./ (4 * pi * (r * 1e-4));
    [f] = second_derivative(V_e);

    % mS / uF = 1 / (uF * kOhm) = mV * msec^-1
    f = Ga / C_m * f;
end

N = 51;
dx = 100 * 10;  % um
x = dx * (-(N-1)/2 : 1 : (N-1)/2);
z = 1;          % mm
D = 10;
d = 0.7 * D;

figure()
for i_elec = [-500 500 -2500 2500]
    [f] = activation_func(x, z, i_elec, 0.3, 0.055, d);
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
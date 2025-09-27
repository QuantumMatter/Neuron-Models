% Compares my Matlab implementation to Evan's Java implementation
% Specifically, it compares the following
%   - Channel dynamics, forward and backward rates vs Vm
%       - (a_m, b_m, ...) vs (Vm, T)
%   - Current for a given state
%       - (i_Na, i_K, i_L) vs (Vm, m, h, n)
%   - Response to constant stimulus

clear java;
javaaddpath('.\+ref\');

% function validate_dynamics()

    V = -100:0.1:100; 
    V = V / 1000; % mV to V
    T = GlobalConstants.T;

    % Preallocate
    david = zeros(length(V), 6);
    evan = zeros(length(V), 6);
    
    for i = 1:length(V)
        v = V(i);

        % David's implementation
        [a_m, b_m] = dynamics.m_dynamics(v * 1000, T);
        [a_h, b_h] = dynamics.h_dynamics(v * 1000, T);
        [a_n, b_n] = dynamics.n_dynamics(v * 1000, T);

        david(i, :) = [a_m, b_m, a_h, b_h, a_n, b_n];

        % Evan's implementation
        a_n = F_PotassiumChannel.alpha(v + GlobalConstants.Vrest);
        b_n = F_PotassiumChannel.beta(v + GlobalConstants.Vrest);
        a_h = F_SodiumChannel.alpha_h(v + GlobalConstants.Vrest);
        b_h = F_SodiumChannel.beta_h(v + GlobalConstants.Vrest);
        a_m = F_SodiumChannel.alpha_m(v + GlobalConstants.Vrest);
        b_m = F_SodiumChannel.beta_m(v + GlobalConstants.Vrest);
        evan(i, :) = [a_m, b_m, a_h, b_h, a_n, b_n];
    end

    % Calculate the differences
    error = sum((david - evan).^2);
    fprintf('Error: %f\n', error);

    % Plot the results as a surface
    surf(david - evan);

% print out the V values where there are NaN
for i = 1:length(V)
    if any(isnan(david(i, :))) || any(isnan(evan(i, :)))
        fprintf('NaN at V = %f mV\n', V(i) * 1000);
    end
end

% end

% validate_dynamics();
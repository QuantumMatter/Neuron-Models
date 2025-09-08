c_m = 1;        % uF/cm^2; membrane capacitance per unit area
gNa_bar = 120;  % mS/cm^2; maximum sodium conductance
gK_bar = 36;    % mS/cm^2; maximum potassium conductance
gL_bar = 0.3;   % mS/cm^2; maximum leak conductance

na_conc_in = 12;    % mM; intra-cellular sodium concentration
na_conc_out = 145;  % mM; extra-cellular sodium concentration
k_conc_in = 140;    % mM; intra-cellular potassium concentration
k_conc_out = 5;     % mM; extra-cellular potassium concentration
cl_conc_in = 4;     % mM; intra-cellular chloride concentration
cl_conc_out = 110;  % mM; extra-cellular chloride concentration

dt = 1e-3;      % time step;                        msec
t = 0:dt:15;  % time vector;                        msec
i_start = 1;    % stimulation start;                msec
i_duration = 500 * 1e-3; % stimulation duration;    msec
i_mask = stepfun(t, i_start) - stepfun(t, i_start + i_duration);
i_elec = -500 * i_mask;

z = 1;
N = 51;
L = 2.5; % Node length;     um
D = 10; % Neuron diameter;  um
rho_e = 0.3;    % Extra-cellular resistivity; kOhm * cm
rho_i = 0.055;  % Intra-cellular resistivity; kOhm * cm

[V, f, Iinj] = Rattay( ...
    t, i_elec, z, ...
    N, L, D, rho_e, rho_i, ...
    gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out ...
);

x = -((N-1)/2):1:(N-1)/2;

figure('Name', 'Stimulus')
plot(t, i_elec)

% Activation function during stimulation
figure()
plot(x, f(:, (i_start / dt)+2))

figure('Name', 'Activation Function')
mesh(t, x, f)

% disp(size(x))
% disp(size(t))
% disp(size(V))

figure('Name', 'Membrane Potential for Central Neuron')
plot(t, V((N+1)/2,:))

figure('Name', 'Membrane Potential for Full Axon')
mesh(t, x, V)
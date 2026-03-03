% Find the threshold (minimum) current density amplitude that will elicit
% an AP, w/ stim parameters pw = 500us, pd = 2ms
% Resolution of the search should be 1 uA

% Provide a plot of the resulting membrane potential as a function of
% time from 0 to 20 ms for the following conditions
% 1 - The stimulus just below threshold
% 2 - The stimulus at threshold

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

t = 0:1e-3:20;      % time vector;      msec
pd = 2;             % pulse delay;      usec
pw = 0.5;           % pulse width;      usec

pa_space = 1:100;           % pluse amplitude;  uA/cm^2

v_traces = zeros([length(pa_space) length(t)]);
indicator = zeros([length(pa_space) 1]);

for ii = 1:length(pa_space)
    pa = pa_space(ii);
    [V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa, pw, gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out);
    v_traces(ii,:) = V;

    if any(V > 0)
        indicator(ii) = 1;
    end
end

%%

sub_idx = find(diff(indicator), 1);
thresh_idx = sub_idx + 1;

fprintf("Sub-Treshold: %.0f uA/cm^2\r\n", pa_space(sub_idx));
fprintf("Treshold: %.0f uA/cm^2\r\n", pa_space(thresh_idx));

%%

figure('Name', 'Thresholding')

plot(t, v_traces(sub_idx, :))
hold on
plot(t, v_traces(thresh_idx, :))

legend(sprintf('Sub-Threshold: %0.f uA/cm^2', pa_space(sub_idx)), sprintf('Threshold: %.0f uA/cm^2', pa_space(thresh_idx)))
ylabel('V (mV)')
xlabel('Time (ms)')
title('Hodgkin–Huxley: V(t) and I_{inj}(t)')
grid on
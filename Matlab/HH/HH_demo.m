c_m = 1;        % uF/cm^2; membrane capacitance per unit area
gNa_bar = 120;  % mS/cm^2; maximum sodium conductance
gK_bar = 36;    % mS/cm^2; maximum potassium conductance
gL_bar = 0.3;   % mS/cm^2; maximum leak conductance

% Problem set values
na_conc_in = 12;    % mM; intra-cellular sodium concentration
na_conc_out = 145;  % mM; extra-cellular sodium concentration
k_conc_in = 140;    % mM; intra-cellular potassium concentration
k_conc_out = 5;     % mM; extra-cellular potassium concentration
cl_conc_in = 4;     % mM; intra-cellular chloride concentration
cl_conc_out = 110;  % mM; extra-cellular chloride concentration

% % Lecture values
% na_conc_in = 5;    % mM; intra-cellular sodium concentration
% na_conc_out = 145;  % mM; extra-cellular sodium concentration
% k_conc_in = 140;    % mM; intra-cellular potassium concentration
% k_conc_out = 5;     % mM; extra-cellular potassium concentration
% cl_conc_in = 4;     % mM; intra-cellular chloride concentration
% cl_conc_out = 110;  % mM; extra-cellular chloride concentration

% R = 8.314;      % mJ/K*mmol;    gas constant
% T = 310;        % K;            temperature
% F = 96.4853;    % C/mmol;       Faraday's constant

t = 0:1e-3:20;    % time vector;      msec
pd = 2;          % pulse delay;      usec
pa = 100;           % pluse amplitude;  uA/cm^2
pw = 0.5;           % pulse width;      usec

[V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa, pw, gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out);

% % plot(t, V)
% % hold on
% plot(t, Iinj)
% xlabel('Time (ms)');
% ylabel('Membrane Potential (mV)');
% title('Hodgkin-Huxley Model Simulation');
% grid on;
% % legend('Membrane Potential', 'Injected Current');

%% FIGURE 1 — Membrane potential and injected current
figure('Name','Membrane Potential & Injected Current');
yyaxis left
plot(t, V, 'LineWidth', 1.6)
ylabel('V (mV)')
yyaxis right
plot(t, Iinj, '--', 'LineWidth', 1.4)
ylabel('I_{inj} (\muA/cm^2)')
xlabel('Time (ms)')
title('Hodgkin–Huxley: V(t) and I_{inj}(t)')
grid on
legend('V', 'I_{inj}', 'Location', 'best')

%% FIGURE 2 — Ionic currents
figure('Name','Ionic Currents');
plot(t, INa, 'LineWidth', 1.6); hold on
plot(t, IK,  'LineWidth', 1.6);
plot(t, IL,  'LineWidth', 1.6); hold off
xlabel('Time (ms)')
ylabel('Current (\muA/cm^2)')
title('Ionic Currents: I_{Na}, I_{K}, I_{L}')
legend('I_{Na}', 'I_{K}', 'I_{L}', 'Location', 'best')
grid on

%% FIGURE 3 — Gate probabilities
figure('Name','Gate Probabilities');
plot(t, n, 'LineWidth', 1.6); hold on
plot(t, m, 'LineWidth', 1.6); 
plot(t, h, 'LineWidth', 1.6); hold off
xlabel('Time (ms)')
ylabel('Probability')
title('Gates: n, m, h')
legend('n', 'm', 'h', 'Location', 'best')
ylim([0 1])
grid on
 
% P_Na = 0.03;    % Resting permeability of sodium
% P_K = 1;        % Resting permeability of potassium
% P_Cl = 0.1;     % Resting permeability of chloride
% 
% % Resting potential
% V_rest = (R*T/F) * log( ...
%     (P_Na*na_conc_out + P_K*k_conc_out + P_Cl*cl_conc_in ) ...
%     / ...
%     (P_Na*na_conc_in  + P_K*k_conc_in  + P_Cl*cl_conc_out) ...
% );
% 
% E_Na = (R*T/F) * log(na_conc_out / na_conc_in);
% E_K = (R*T/F) * log(k_conc_out / k_conc_in);
% 
% % % Lecture values
% % E_Na = 50;
% % E_K = -77;
% 
% 
% % Calculates the forward rate constant for the 'm' gate
% % Units are ms^-1
% % V - membrane potential in millivolts
% function [alpha_m] = calc_alpha_m(V)
%     alpha_m = (0.1 * (V + 40)) / (1 - exp(-(V+40)/10));
% end
% 
% % Calculates the reverse rate constant for the 'm' gate
% % Units are ms^-1
% % V - membrane potentional in millivolts
% function [beta_m] = calc_beta_m(V)
%     beta_m = 4 * exp(-(V+65)/18);
% end
% 
% function [m] = calc_m(V)
%     alpha_m = calc_alpha_m(V);
%     beta_m = calc_beta_m(V);
%     m = alpha_m / (alpha_m + beta_m);
% end
% 
% function [alpha_h] = calc_alpha_h(V)
%     alpha_h = 0.07 * exp(-(V+65)/20);
% end
% 
% function [beta_h] = calc_beta_h(V)
%     beta_h = 1 / (1 + exp(-(V+35)/10));
% end
% 
% function [h] = calc_h(V)
%     alpha_h = calc_alpha_h(V);
%     beta_h = calc_beta_h(V);
%     h = alpha_h / (alpha_h + beta_h);
% end
% 
% function [alpha_n] = calc_alpha_n(V)
%     alpha_n = (0.01 * (V + 55)) / (1 - exp(-(V+55)/10));
% end
% 
% function [beta_n] = calc_beta_n(V)
%     beta_n = 0.125 * exp(-(V + 65) / 80);
% end
% 
% function [n] = calc_n(V)
%     alpha_n = calc_alpha_n(V);
%     beta_n = calc_beta_n(V);
%     n = alpha_n / (alpha_n + beta_n);
% end
% 
% m_0 = calc_m(V_rest);
% h_0 = calc_h(V_rest);
% n_0 = calc_n(V_rest);
% 
% function [gNa] = calc_gNa(gNa_bar, m, h)
%     gNa = gNa_bar * power(m, 4) * h;
% end
% 
% function [gK] = calc_gK(gK_bar, n)
%     gK = gK_bar * power(n, 4);
% end
% 
% gNa_0 = calc_gNa(gNa_bar, m_0, h_0);
% gK_0 = calc_gK(gK_bar, n_0);
% 
% E_l = V_rest + ((gK_0 * (V_rest - E_K)) + (gNa_0 * (V_rest - E_Na)) / gL_bar)
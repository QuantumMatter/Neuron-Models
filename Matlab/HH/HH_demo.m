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
pa = 100;           % pluse amplitude;  uA/cm^2
pw = 0.5;           % pulse width;      usec

[V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa, pw, gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out);

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


%% FIGURE 4 - Time series response for increasing stimulus

t = 0:1e-3:50;
pd = 1;
pa = -10:0.1:20;
pw = 50000;

V_multiple = zeros(length(t), length(pa));

for i = 1:length(pa)

    [V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa(i), pw, gNa_bar, gK_bar, gL_bar, k_conc_out, na_conc_out);
    V_multiple(:,i) = V;

end

fig = figure();
fig.set('Name', 'Animation')
ax1 = nexttile;
hlv = plot(ax1, t, V_multiple(:,1), 'LineWidth', 1.4);
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Membrane Potential Response to Constant Stimulus Step')
ylim([-110 80])
grid on
sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.96 0.7 0.03], 'Min',min(pa), 'Max',max(pa), 'Value',0);
lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.96 0.12 0.03], 'String','init', 'BackgroundColor','w');

function UpdateStim(k, V_multiple, hlv, lbl, pa)

    i = round((k - pa(1)) / (pa(2) - pa(1))) + 1;

    set(hlv, 'YData', V_multiple(:,i));
    set(lbl, 'String', sprintf("Stimulus = %.01f uA/cm^2", pa(i)));

end

sld.addlistener('ContinuousValueChange', @(src,~) UpdateStim(src.Value, V_multiple, hlv, lbl, pa));
UpdateStim(0, V_multiple, hlv, lbl, pa);

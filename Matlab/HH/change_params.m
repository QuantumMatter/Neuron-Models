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

% V_multiple = zeros(length(t), length(pa));

gNa_bar = 100:5:140;
gL_bar = 0.15:0.1:0.45;
gK_bar = 20:4:52;
na_conc_out = 120:5:160;  % mM; extra-cellular sodium concentration
k_conc_out = 2:1:10;     % mM; extra-cellular potassium concentration

traces = zeros(length(t), length(gNa_bar), length(gL_bar), length(gK_bar), length(na_conc_out), length(k_conc_out));

for ii = 1:length(gNa_bar)
    for jj = 1:length(gL_bar)
        for kk = 1:length(gK_bar)
            for xx = 1:length(na_conc_out)
                for yy = 1:length(k_conc_out)

                    [V, n, m, h, INa, IK, IL, Iinj] = HH(t, pd, pa, pw, gNa_bar(ii), gK_bar(kk), gL_bar(jj), k_conc_out(yy), na_conc_out(xx));
                    traces(:,ii,jj,kk,xx,yy) = V;

                end
            end
        end
    end
end

%%
clc
clear state

fig = figure();
fig.set('Name', 'Animation')
ax1 = gca;
ax1.Position = [0.05 0.05 0.9 0.7];

hlv = plot(ax1, t, traces(:,1,1,1,1,1), 'LineWidth', 1.4);

xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Membrane Potential Response to Constant Stimulus Step')
ylim([-110 80])
grid on

gNa_sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.96 0.7 0.03], 'Min',min(gNa_bar), 'Max',max(gNa_bar), 'Value',gNa_bar(1));
gNa_lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.96 0.12 0.03], 'String','init', 'BackgroundColor','w');
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.03, 0.96, 0.10, 0.03], 'String', 'gNa');

gK_sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.92 0.7 0.03], 'Min',min(gK_bar), 'Max',max(gK_bar), 'Value',mean(gK_bar));
gK_lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.92 0.12 0.03], 'String','init', 'BackgroundColor','w');
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.03, 0.92, 0.10, 0.03], 'String', 'gK');

gL_sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.88 0.7 0.03], 'Min',min(gL_bar), 'Max',max(gL_bar), 'Value',mean(gL_bar));
gL_lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.88 0.12 0.03], 'String','init', 'BackgroundColor','w');
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.03, 0.88, 0.10, 0.03], 'String', 'gL');

whos na_conc_out
disp([min(na_conc_out) max(na_conc_out) mean(na_conc_out)])
naout_sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.84 0.7 0.03], 'Min',min(na_conc_out), 'Max',max(na_conc_out), 'Value',mean(na_conc_out));
naout_lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.84 0.12 0.03], 'String','init', 'BackgroundColor','w');
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.03 0.84 0.10 0.03], 'String', '[Na] out');

kout_sld = uicontrol('Style','slider', 'Units','normalized', ...
    'Position',[0.15 0.8 0.7 0.03], 'Min',min(k_conc_out), 'Max',max(k_conc_out), 'Value',mean(k_conc_out));
kout_lbl = uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.87 0.8 0.12 0.03], 'String','init', 'BackgroundColor','w');
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.03 0.8 0.10 0.03], 'String', '[K] out');

state = SimState(1, 1, 1, 1, 1);

function UpdateStim(k, which, traces, hlv, lbl, data, state)

    display(k);

    i = round((k - data(1)) / (data(2) - data(1))) + 1;
    i = max(1, min(numel(data), i));

    if which == 1
        state.gNa_ii = i;
    elseif which == 2
        state.gL_jj = i;
    elseif which == 3
        state.gK_kk = i;
    elseif which == 4
        state.naout_xx = i;
    elseif which == 5
        state.kout_yy = i;
    end

    % display(which);
    % display(state);
    
    set(hlv, 'YData', traces(:, state.gNa_ii, state.gL_jj, state.gK_kk, state.naout_xx, state.kout_yy));
    set(lbl, 'String', sprintf("%.02f", data(i)));

end

L = cell(5,1);

L{1} = addlistener(gNa_sld,   'Value', 'PostSet', @(~,evt) UpdateStim(evt.AffectedObject.Value, 1, traces, hlv, gNa_lbl,   gNa_bar,     state));
L{2} = addlistener(gL_sld,    'Value', 'PostSet', @(~,evt) UpdateStim(evt.AffectedObject.Value, 2, traces, hlv, gL_lbl,    gL_bar,      state));
L{3} = addlistener(gK_sld,    'Value', 'PostSet', @(~,evt) UpdateStim(evt.AffectedObject.Value, 3, traces, hlv, gK_lbl,    gK_bar,      state));
L{4} = addlistener(naout_sld, 'Value', 'PostSet', @(~,evt) UpdateStim(evt.AffectedObject.Value, 4, traces, hlv, naout_lbl, na_conc_out, state));
L{5} = addlistener(kout_sld,  'Value', 'PostSet', @(~,evt) UpdateStim(evt.AffectedObject.Value, 5, traces, hlv, kout_lbl,  k_conc_out,  state));

fig.UserData.listeners = L;

UpdateStim(mean(gNa_bar), 1, traces, hlv, gNa_lbl, gNa_bar, state);
UpdateStim(mean(gL_bar), 2, traces, hlv, gL_lbl, gL_bar, state);
UpdateStim(mean(gK_bar), 3, traces, hlv, gK_lbl, gK_bar, state);
UpdateStim(mean(na_conc_out), 4, traces, hlv, naout_lbl, na_conc_out, state);
UpdateStim(mean(k_conc_out), 5, traces, hlv, kout_lbl, k_conc_out, state);

test = min(na_conc_out);
naout_sld.Min = test;
display(test);
display(naout_sld.Min);
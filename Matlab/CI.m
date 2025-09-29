function [] = CI(p)

    close all;
    clearvars -except p path

    cwd = fileparts(mfilename('fullpath'));

    addpath(fullfile(cwd, 'HH'));
    addpath(fullfile(cwd, 'Rattay'));
    addpath(fullfile(cwd, 'SENN'));
    addpath(fullfile(cwd, 'utils'));

    HH_demo;
    save_all_open_figs(fullfile(p, 'HH'), '')

    close all;
    clearvars -except p path
    Rattay_demo;
    save_all_open_figs(fullfile(p, 'Rattay'), '')

    close all;
    clearvars -except p path
    demo_propagate;
    save_all_open_figs(fullfile(p, 'SENN'), '')

    close all
    clearvars

end
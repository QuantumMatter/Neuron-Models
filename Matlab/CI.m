function [] = CI(p)

    close all;
    clearvars -except p

    cwd = fileparts(mfilename('fullpath'));

    addpath(fullfile(cwd, 'HH'));
    addpath(fullfile(cwd, 'Rattay'));
    addpath(fullfile(cwd, 'utils'));

    HH_demo;
    save_all_open_figs(fullfile(p, 'HH'), '')

    close all;
    clearvars -except p
    Rattay_demo;
    save_all_open_figs(fullfile(p, 'Rattay'), '')

    close all
    clearvars

end
function [] = CI(path)

    close all;
    clearvars -except path

    addpath('./HH');
    addpath('./Rattay');
    addpath('./utils');

    HH_demo;
    save_all_open_figs(sprintf('%s/HH', path), '')

    close all;
    clearvars -except path
    Rattay_demo;
    save_all_open_figs(sprintf('%s/Rattay', path), '')

    close all
    clearvars

end
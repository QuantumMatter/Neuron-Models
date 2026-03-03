classdef SimState < handle
    properties
        gNa_ii = 1;
        gL_jj = 1;
        gK_kk = 1;
        naout_xx = 1;
        kout_yy = 1;
    end
    methods
        function obj = SimState(ii, jj, kk, xx, yy)
            obj.gNa_ii = ii;
            obj.gL_jj = jj;
            obj.gK_kk = kk;
            obj.naout_xx = xx;
            obj.kout_yy = yy;
        end
    end
end
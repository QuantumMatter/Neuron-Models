% Indices go from distal to medial
%   Index 1 is the distal node
%   Index N is the medial node
function [V] = Hayden(...
    t, i_elec, z, ...
    N, ...
    cm, dk, lk, rhok, L, ...
    gL)

    arguments (Input)
        t (1,:) % test
        i_elec (1,:) {mustBeEqmustBeEqualSize(t, i_elec)}  %
        z % Distance to electrode
    end

    arguments (Output)
        V
    end

end


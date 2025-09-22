
% Approximates the function
%      x
% ------------
% 1 - exp(x*a)
%
% around x=0 to fill in the discontinuity
% This was done by finding the maclaurin
% expansion and then simplifying
function [y] = vtrap(x, alpha)
    if abs(x*alpha) > 1e-6
        y = x / (1 - exp(-x*alpha));
        return;
    end
    y = alpha + x*alpha^2/2 + x.^2*alpha^3/6;
    y = -1./y;
end

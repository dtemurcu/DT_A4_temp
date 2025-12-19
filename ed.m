function [g, C, lambda] = ed(co, a, b, gmin, gmax, d, toler)
% Deniz Temurcu 261089503
% This function performs economic dispatch using the lambda-iteration method

% Our inputs:
% co     is the vector of constant cost coefficients ($/h)
% a      is the vector of linear cost coefficients ($/MWh)
% b      is the vector of quadratic cost coefficients ($/MW^2h)
% gmin   is the vector of minimum generation limits (MW)
% gmax   is the vector of maximum generation limits (MW)
% d      is the total power demand (MW)
% toler  is the power balancing tolerance (MW)

% Our outputs:
% g      is the vector of generation outputs (MW)
% C      is the total generation cost ($/h)
% lambda is the power balance lagrange multiplier ($/MWh)

% make sure the sizes match up
n = length(co);
if any([numel(a),numel(b),numel(gmin),numel(gmax)] ~= n)
    error('co, a, b, gmin, gmax must all be length n.');
end
co = co(:); a = a(:); b = b(:); gmin = gmin(:); gmax = gmax(:); % ensure column vectors
% feasibility check
if d > sum(gmax) || d < sum(gmin)
    error('infeasible demand: d is outside global generation limits.');
end
% initialize
lambda = mean(a);    % initial guess based on average linear costs
iter = 0;
maxiter = 1000;
% lambda iteration loop
while iter < maxiter
    iter = iter + 1;
    % calculate tentative generation
    g = (lambda - a) ./ b;              % unconstrained dispatch
    g = max(g, gmin);                   % apply lower limits
    g = min(g, gmax);                   % apply upper limits
    % check mismatch
    mismatch = d - sum(g);
    if abs(mismatch) <= toler
        break
    end
    % update lambda using gradient method from slides
    unconstrained = (g > gmin) & (g < gmax);
    if any(unconstrained)
        % slope = sum(1/b_i) for active marginal units
        slope = sum(1 ./ b(unconstrained));
        dlambda = mismatch / slope;
    else
        % fallback: use aggregate slope if all units hit limits to avoid stall
        slope_all = sum(1 ./ b);
        dlambda = mismatch / slope_all;
    end
    lambda = lambda + dlambda;
end
if iter >= maxiter
    error('economic dispatch failed to converge.');
end
% calculate final cost
C = sum(co + a .* g + 0.5 * b .* g.^2);
end
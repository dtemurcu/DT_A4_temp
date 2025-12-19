function [u, g, C, lambda] = uc(co, a, b, gmin, gmax, d, toler)
% Deniz Temurcu 261089503
% This function performs static unit commitment using full enumeration

% Our inputs:
% co     is the vector of constant cost coefficients ($/h)
% a      is the vector of linear cost coefficients ($/MWh)
% b      is the vector of quadratic cost coefficients ($/MW^2h)
% gmin   is the vector of minimum generation limits (MW)
% gmax   is the vector of maximum generation limits (MW)
% d      is the total power demand (MW)
% toler  is the power balancing tolerance (MW)

% Our outputs:
% u      is the vector of generation statuses (binary)
% g      is the vector of generation outputs (MW)
% C      is the total generation cost ($/h)
% lambda is the power balance lagrange multiplier ($/MWh)

% make sure the sizes match up
n = length(co);
co = co(:); a = a(:); b = b(:); gmin = gmin(:); gmax = gmax(:);
% initialize
min_C = Inf;
best_u = zeros(n, 1);
best_g = nan(n, 1);
best_lambda = nan;
num_comb = 2^n - 1; % total combinations excluding all-off
% enumeration loop
for k = 1:num_comb
    % generate binary status vector
    bin_str = dec2bin(k, n);
    u_trial = zeros(n, 1);
    for i = 1:n
        u_trial(i) = str2double(bin_str(n - i + 1));
    end
    idx = find(u_trial == 1); % active unit indices
    if isempty(idx), continue; end
    % capacity check
    if d > sum(gmax(idx)) || d < sum(gmin(idx))
        continue; % skip if infeasible capacity
    end
    % run economic dispatch on subset
    try
        [g_sub, var_C, lam_sub] = ed(co(idx), a(idx), b(idx), gmin(idx), gmax(idx), d, toler);
        total_C = var_C + sum(co(idx)); % add fixed costs
        if total_C < min_C
            min_C = total_C;
            best_u = u_trial;
            best_g = zeros(n, 1);
            best_g(idx) = g_sub;
            best_lambda = lam_sub;
        end
    catch
        continue; % skip if ed fails
    end
end
% assign outputs
u = best_u;
g = best_g;
C = min_C;
lambda = best_lambda;
if isinf(min_C)
    warning('no feasible unit commitment found.');
end
end
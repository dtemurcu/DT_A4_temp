function [g, C, lambda, LMP, CongestionSurplus] = dcscopf(line_data, cost_data, gmin, gmax, d, ngen)
% Deniz Temurcu 261089503
% This function performs dc-scopf using quadratic programming

% Our inputs:
% line_data : matrix [from, to, x, fmax]
% cost_data : matrix [co, a, b]
% gmin, gmax : generation limits (MW)
% d         : nodal demand vector (MW)

% ngen      : generator bus indices
% Our outputs:
% g         : optimal generation (MW)
% C         : total cost ($/h)
% lambda    : system lambda ($/MWh)
% LMP       : locational marginal prices ($/MWh)
% CongestionSurplus : total congestion rent ($/h)

% parse inputs
nfrom = line_data(:,1); nto = line_data(:,2);
x_line = line_data(:,3); fmax = line_data(:,4);
co = cost_data(:,1); a = cost_data(:,2); b = cost_data(:,3);
% setup dimensions and b-matrix
nb = length(d); ng = length(co); nl = length(nfrom); ref_bus = 1;
B = zeros(nb, nb);
for k = 1:nl
    i = nfrom(k); j = nto(k); y = 1/x_line(k);
    B(i,i) = B(i,i) + y; B(j,j) = B(j,j) + y;
    B(i,j) = B(i,j) - y; B(j,i) = B(j,i) - y;
end
% optimization setup (x = [g; theta])
H = blkdiag(diag(b), zeros(nb, nb));
f = [a; zeros(nb, 1)];
% equality constraints: P_gen - B*theta = P_load
Cg = zeros(nb, ng);
for k = 1:ng, Cg(ngen(k), k) = 1; end
Aeq = [Cg, -B; zeros(1, ng+nb)];
Aeq(end, ng+ref_bus) = 1; % ref bus angle = 0
beq = [d; 0];
% inequality constraints: flow <= fmax
A = zeros(2*nl, ng+nb); bvec = zeros(2*nl, 1);
for k = 1:nl
    i = nfrom(k); j = nto(k); y = 1/x_line(k);
    % flow_ij = (theta_i - theta_j) / x
    A(k, ng+i) = y; A(k, ng+j) = -y; bvec(k) = fmax(k);       % +flow limit
    A(k+nl, ng+i) = -y; A(k+nl, ng+j) = y; bvec(k+nl) = fmax(k); % -flow limit
end
% solve qp
lb = [gmin; -Inf(nb, 1)]; ub = [gmax; Inf(nb, 1)];
opts = optimset('Algorithm','interior-point-convex','Display','off');
[x, ~, exitflag, ~, lam] = quadprog(H, f, A, bvec, Aeq, beq, lb, ub, [], opts);
if exitflag <= 0, error('optimization failed.'); end
% extract results
g = x(1:ng);
C = sum(co + a.*g + 0.5*b.*g.^2);
% lmp calculation (lmp = dCost/dLoad)
% matlab returns lagrange multipliers for Aeq*x=beq such that grad(L)=0
% standard interpretation: lmp is negative of lambda_eq for P-B*theta=d formulation
LMP = -lam.eqlin(1:nb);
lambda = LMP(ref_bus);
% congestion surplus
P_gen_bus = Cg * g;
CongestionSurplus = sum(LMP .* d) - sum(LMP .* P_gen_bus);
end
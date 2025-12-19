function [delta, V, N, time] = fdwlsse(nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)
% Deniz Temurcu 261089503
% This function performs fast-decoupled weighted least squares state estimation

% Our inputs:
% nfrom, nto : vectors of line connectivity
% r, x, b    : vectors of line resistance, reactance, and charging susceptance (p.u.)
% Pinj, Qinj : matrices [bus_idx, val, sigma]
% Pflow, Qflow: matrices [from, to, val, sigma]
% Vnode      : matrix [bus_idx, val, sigma]
% toler      : convergence tolerance
% maxiter    : maximum iterations

% Our outputs:
% delta     : estimated voltage angles (rad)
% V         : estimated voltage magnitudes (p.u.)
% N         : iteration count
% time      : cpu time (s)

tstart = tic;

% Ensure column vectors for network data
nfrom = nfrom(:); nto = nto(:); r = r(:); x = x(:); b_ch = b(:);

nb = max([nfrom; nto]); nl = length(nfrom); ref_bus = 1;

% build exact y-bus for mismatch calc
Y = zeros(nb, nb);
for k = 1:nl
    i = nfrom(k); j = nto(k); z = r(k) + 1j*x(k); y = 1/z; bc = 1j*b_ch(k)/2;
    Y(i,i) = Y(i,i) + y + bc; Y(j,j) = Y(j,j) + y + bc;
    Y(i,j) = Y(i,j) - y; Y(j,i) = Y(j,i) - y;
end
G_mat = real(Y); B_mat = imag(Y);

% build measurement vectors and weights
% Pinj: [bus, val, sig] -> zP includes Pinj val
% Pflow: [from, to, val, sig] -> zP includes Pflow val
zP = [Pinj(:,2); Pflow(:,3)];
sigP = [Pinj(:,3); Pflow(:,4)];
WP = diag(sigP.^-2);

% Qinj: [bus, val, sig]
% Qflow: [from, to, val, sig]
% Vnode: [bus, val, sig]
zQ = [Qinj(:,2); Qflow(:,3); Vnode(:,2)];
sigQ = [Qinj(:,3); Qflow(:,4); Vnode(:,3)];
WQ = diag(sigQ.^-2);

npi = size(Pinj,1); npf = size(Pflow,1);
nqi = size(Qinj,1); nqf = size(Qflow,1);
nvi = size(Vnode,1);

% build constant jacobians HP and HQ using 1/x approx
HP = zeros(npi+npf, nb); HQ = zeros(nqi+nqf+nvi, nb);

% P-inj: sum(1/x)
for k = 1:npi 
    i = Pinj(k,1); idx = (nfrom==i | nto==i);
    HP(k, i) = sum(1./x(idx));
    for l = find(idx)', nbr = setdiff([nfrom(l), nto(l)], i); HP(k, nbr) = -1/x(l); end
end

% P-flow: 1/x
for k = 1:npf 
    i = Pflow(k,1); j = Pflow(k,2); 
    l = find((nfrom==i & nto==j) | (nfrom==j & nto==i), 1);
    HP(npi+k, i) = 1/x(l); HP(npi+k, j) = -1/x(l);
end

% Q-inj: same as P-inj
for k = 1:nqi 
    i = Qinj(k,1); idx = (nfrom==i | nto==i);
    HQ(k, i) = sum(1./x(idx));
    for l = find(idx)', nbr = setdiff([nfrom(l), nto(l)], i); HQ(k, nbr) = -1/x(l); end
end

% Q-flow: 1/x
for k = 1:nqf 
    i = Qflow(k,1); j = Qflow(k,2);
    l = find((nfrom==i & nto==j) | (nfrom==j & nto==i), 1);
    HQ(nqi+k, i) = 1/x(l); HQ(nqi+k, j) = -1/x(l);
end

% V-mag
for k = 1:nvi, HQ(nqi+nqf+k, Vnode(k,1)) = 1; end 

HP(:, ref_bus) = []; % remove ref bus col

% gain matrices
GP = HP'*WP*HP; GQ = HQ'*WQ*HQ;
if rcond(GP)<1e-12, GP = GP + 1e-4*eye(size(GP)); warning('regularizing GP'); end
if rcond(GQ)<1e-12, GQ = GQ + 1e-4*eye(size(GQ)); warning('regularizing GQ'); end

% estimation loop
V = ones(nb,1); delta = zeros(nb,1); N = 0;
for iter = 1:maxiter
    % calculate exact mismatches h(x)
    Pcal = zeros(nb,1); Qcal = zeros(nb,1);
    for i=1:nb, for j=1:nb, if Y(i,j)~=0
        t = delta(i)-delta(j); m = V(i)*V(j);
        Pcal(i) = Pcal(i) + m*(G_mat(i,j)*cos(t)+B_mat(i,j)*sin(t));
        Qcal(i) = Qcal(i) + m*(G_mat(i,j)*sin(t)-B_mat(i,j)*cos(t));
    end, end, end

    hP = zeros(npi+npf, 1); hQ = zeros(nqi+nqf+nvi, 1);
    
    % hP for Injections
    for k=1:npi, hP(k) = Pcal(Pinj(k,1)); end
    
    % hP for Flows
    for k=1:npf
        i=Pflow(k,1); j=Pflow(k,2); l=find((nfrom==i&nto==j)|(nfrom==j&nto==i),1);
        y=1/(r(l)+1j*x(l)); t=delta(i)-delta(j);
        hP(npi+k) = real(y)*V(i)^2 - V(i)*V(j)*(real(y)*cos(t)+imag(y)*sin(t));
    end

    % hQ for Injections
    for k=1:nqi, hQ(k) = Qcal(Qinj(k,1)); end
    
    % hQ for Flows
    for k=1:nqf
        i=Qflow(k,1); j=Qflow(k,2); l=find((nfrom==i&nto==j)|(nfrom==j&nto==i),1);
        y=1/(r(l)+1j*x(l)); t=delta(i)-delta(j);
        hQ(nqi+k) = -V(i)^2*(imag(y)+b_ch(l)/2) - V(i)*V(j)*(real(y)*sin(t)-imag(y)*cos(t));
    end
    
    % hQ for Vnode
    for k=1:nvi, hQ(nqi+nqf+k) = V(Vnode(k,1)); end

    % residuals
    rP = zP - hP; rQ = zQ - hQ;
    if max(abs([rP; rQ])) < toler, N=iter; break; end

    % solve decoupled
    dd = zeros(nb,1);
    dd(2:end) = GP \ (HP'*WP*rP);
    delta = delta + dd;
    V = V + (GQ \ (HQ'*WQ*rQ));
end
if N==0, N=maxiter; end
time = toc(tstart);
end
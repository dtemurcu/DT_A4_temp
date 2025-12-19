function [delta, V, N, time] = fdwlsse(line_data, z_inj, z_flow, z_node, toler, maxiter)
% Deniz Temurcu 261089503
% This function performs fast-decoupled weighted least squares state estimation

% Our inputs:
% line_data : matrix [from, to, r, x, b]
% z_inj     : matrix [bus_idx, P_val, sigP, Q_val, sigQ]
% z_flow    : matrix [from, to, P_val, sigP, Q_val, sigQ]
% z_node    : matrix [bus_idx, V_val, sigV]
% toler     : convergence tolerance

% maxiter   : maximum iterations
% Our outputs:
% delta     : estimated voltage angles (rad)
% V         : estimated voltage magnitudes (p.u.)
% N         : iteration count
% time      : cpu time (s)

% parse network data
tstart = tic;
nfrom = line_data(:,1); nto = line_data(:,2);
r = line_data(:,3); x = line_data(:,4); b_ch = line_data(:,5);
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
Pinj = z_inj(:,1:2); sigPinj = z_inj(:,3);
Qinj = [z_inj(:,1), z_inj(:,4)]; sigQinj = z_inj(:,5);
Pflow = z_flow(:,1:3); sigPflow = z_flow(:,4);
Qflow = [z_flow(:,1:2), z_flow(:,5)]; sigQflow = z_flow(:,6);
Vnode = z_node(:,1:2); sigVnode = z_node(:,3);
zP = [Pinj(:,2); Pflow(:,3)]; sigP = [sigPinj; sigPflow]; WP = diag(sigP.^-2);
zQ = [Qinj(:,2); Qflow(:,3); Vnode(:,2)]; sigQ = [sigQinj; sigQflow; sigVnode]; WQ = diag(sigQ.^-2);
npi = size(Pinj,1); npf = size(Pflow,1); nqi = size(Qinj,1); nqf = size(Qflow,1); nvi = size(Vnode,1);
% build constant jacobians HP and HQ using 1/x approx
HP = zeros(npi+npf, nb); HQ = zeros(nqi+nqf+nvi, nb);
for k = 1:npi % P-inj: sum(1/x)
    i = Pinj(k,1); idx = (nfrom==i | nto==i);
    HP(k, i) = sum(1./x(idx));
    for l = find(idx)', nbr = setdiff([nfrom(l), nto(l)], i); HP(k, nbr) = -1/x(l); end
end
for k = 1:npf % P-flow: 1/x
    i = Pflow(k,1); j = Pflow(k,2); l = find((nfrom==i & nto==j) | (nfrom==j & nto==i), 1);
    HP(npi+k, i) = 1/x(l); HP(npi+k, j) = -1/x(l);
end
for k = 1:nqi % Q-inj: same as P-inj
    i = Qinj(k,1); idx = (nfrom==i | nto==i);
    HQ(k, i) = sum(1./x(idx));
    for l = find(idx)', nbr = setdiff([nfrom(l), nto(l)], i); HQ(k, nbr) = -1/x(l); end
end
for k = 1:nqf % Q-flow: 1/x
    l = find((nfrom==Qflow(k,1) & nto==Qflow(k,2)) | (nfrom==Qflow(k,2) & nto==Qflow(k,1)), 1);
    HQ(nqi+k, Qflow(k,1)) = 1/x(l); HQ(nqi+k, Qflow(k,2)) = -1/x(l);
end
for k = 1:nvi, HQ(nqi+nqf+k, Vnode(k,1)) = 1; end % V-mag
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
    for k=1:npi, hP(k) = Pcal(Pinj(k,1)); end
    for k=1:npf
        i=Pflow(k,1); j=Pflow(k,2); l=find((nfrom==i&nto==j)|(nfrom==j&nto==i),1);
        y=1/(r(l)+1j*x(l)); t=delta(i)-delta(j);
        hP(npi+k) = real(y)*V(i)^2 - V(i)*V(j)*(real(y)*cos(t)+imag(y)*sin(t));
    end
    for k=1:nqi, hQ(k) = Qcal(Qinj(k,1)); end
    for k=1:nqf
        i=Qflow(k,1); j=Qflow(k,2); l=find((nfrom==i&nto==j)|(nfrom==j&nto==i),1);
        y=1/(r(l)+1j*x(l)); t=delta(i)-delta(j);
        hQ(nqi+k) = -V(i)^2*(imag(y)+b_ch(l)/2) - V(i)*V(j)*(real(y)*sin(t)-imag(y)*cos(t));
    end
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
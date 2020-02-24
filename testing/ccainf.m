function ccainf(N,Ny,Nx,Nz,Nw,Npca,dist_parm,sig,nR,nP,deconf_meth,use_nullspace,seed0)
% Monte Carlo evaluation of CCA inference
% ccainf(N,Ny,Nx,Nz,Nw,Npca,dist_parm,sig,nR,nP,deconf_meth,seed0)
%
% N           : Number of subjects
% Ny          : Number of variables in Y
% Nx          : Number of variables in X
% Nz          : Number of nuisance variables (in Z)
% Nw          : Number of nuisance variables (in W). If NaN, then W=Z.
% Npca        : Number of PCA dims to retain (empty or 0 for no PCA)
% dist_parm   : Parameter that determines the distribution of X, Y and Z.
%               - dparm = 0 : normal distribution
%               - 0 < dparm < 1 : categorical (0/1), with the E(x) = dparm
%               - t >= 1 : t distribution with df = dparm
% sig         : Add signal?
%               If univariate, standard deviation of added common part
%               If multivariate: standard deviation separately per
%               column. Should be of length <= min(Nx, Ny)
% nR          : Number of realizations
% nP          : Number of permutations per realization
% deconf_meth : Deconfound-fix method:
%               - 0 : none
%               - 1 : HJ w/ SVD
%               - 2 : HJ w/ Schur
%               - 3 : Theil's BLUSs
%               - 4 : Random deletion of Nz+1 rows
% use_null    : If true, includes the complement of the canonical
%               coefficients (i.e., it uses the null space).
% seed0       : Initial seed for the random number generator (gets
%               incremented at each realization)
%
% All arguments required. Summary of results printed; nothing saved.
%
% An interecept is always used; Nz = 1 means one nuisance variable in
% addition to the intercept. Thus, if e.g., Nz = 3, Z will have 4 cols,
% one being the intercept.
%
%----------------------------------------------------------------------------
% Anderson Winkler, Olivier Renaud, Thomas Nichols
% August 2019

if numel(sig) > 1
    sigstr = numel(sig);
else
    sigstr = sig;
end

fprintf('Simulation parameters: ')
fprintf('N: %d, Ny: %d, Nx: %d, Nz: %d, Nw: %d, Npca: %d, dparm: %d, sig: %f, nR: %d, nP: %d, dmeth: %d, use_null: %d, sd0: %d\n',...
    N,Ny,Nx,Nz,Nw,Npca,dist_parm,sigstr,nR,nP,deconf_meth,use_nullspace,seed0);

if Npca == 0
    Npca = [];
end
if isnan(Nw)
    WeqZ = true;
else
    WeqZ = false;
end

alpha = 0.05;
ppara_rep  = zeros(nR,min([Ny Nx Npca]));
puncRo_rep = zeros(nR,min([Ny Nx Npca]));
puncWo_rep = zeros(nR,min([Ny Nx Npca]));
puncRs_rep = zeros(nR,min([Ny Nx Npca]));
puncWs_rep = zeros(nR,min([Ny Nx Npca]));
pcloRo_rep = zeros(nR,min([Ny Nx Npca]));
pcloWo_rep = zeros(nR,min([Ny Nx Npca]));
pcloRs_rep = zeros(nR,min([Ny Nx Npca]));
pcloWs_rep = zeros(nR,min([Ny Nx Npca]));
pmaxRo_rep = zeros(nR,min([Ny Nx Npca]));
pmaxWo_rep = zeros(nR,min([Ny Nx Npca]));
pmaxRs_rep = zeros(nR,min([Ny Nx Npca]));
pmaxWs_rep = zeros(nR,min([Ny Nx Npca]));
cco1_rep   = zeros(nR,min([Ny Nx Npca]));
ccor_rep   = zeros(nR,min([Ny Nx Npca]));
ccs1_rep   = zeros(nR,min([Ny Nx Npca]));
ccsr_rep   = zeros(nR,min([Ny Nx Npca]));

% For each realization:
for r = 1:nR
    fprintf('Realization %d:\n',r)
    
    % Use rng for repeatability:
    if isoctave
        rand('state',r+seed0); %#ok
    else
        rng(r+seed0);
    end
    
    % Create some random data.
    if dist_parm == 0
        % Normal
        Y = randn(N,Ny);
        X = randn(N,Nx);
        Z = [randn(N,Nz) ones(N,1)];
        if WeqZ
            W  = Z;
            Nw = Nz;
        else
            W  = [randn(N,abs(Nw)) ones(N,1)];
        end
    elseif dist_parm > 0 && dist_parm < 1
        % Bernoulli (binary)
        Y = double(rand(N,Ny) > dist_parm);
        X = double(rand(N,Nx) > dist_parm);
        Z = [(randn(N,Nz) > dist_parm) ones(N,1)];
        if WeqZ
            W  = Z;
            Nw = Nz;
        else
            W  = [(randn(N,abs(Nw)) > dist_parm) ones(N,1)];
        end
    elseif dist_parm >= 1
        % Student's t (kurtotic)
        Y = trnd(dist_parm,N,Ny);
        X = trnd(dist_parm,N,Nx);
        Z = [trnd(dist_parm,N,Nz) ones(N,1)];
        if WeqZ
            W  = Z;
            Nw = Nz; 
        else
            W  = [trnd(dist_parm,N,abs(Nw)) ones(N,1)];
        end
    end
    
    % % ========
    % % This block is to introduce correlations between Z and W.
    % % Comment out to have these two random.
    % Nc = min(Nz,Nw)-1;
    % C = randn(N,Nc);
    % Z(:,1:Nc) = Z(:,1:Nc)/2 + C;
    % W(:,1:Nc) = W(:,1:Nc)/2 + C;
    % [~,~,rtmp,~,~] = canoncorr(Z(:,1:Nz),W(:,1:Nw));
    % disp(rtmp);
    % %disp(corr(Z(:,1:Nz),W(:,1:Nw)));
    % % ========
    
    if Nw < 0
        Nc = min(Nz,abs(Nw))-1;
        C  = randn(N,Nc);
        Z(:,1:Nc) = Z(:,1:Nc)/2 + C;
        W(:,1:Nc) = W(:,1:Nc)/2 + C;
    end
    
    % Add common signal
    lsig = length(sig);
    if lsig > min(Nx,Ny)
        error('Length of sig must be <= than Nx and Ny.');
    end
    if lsig > 1
        sigstr = sprintf('Power for k=(1:%d), FPR for the rest.',lsig);
        for k = 1:lsig
            tmp = randn(N,1)*sig(k);
            Y(:,k) = Y(:,k) + tmp;
            X(:,k) = X(:,k) + tmp;
        end
    elseif sig > 0
        sigstr = 'Power';
        tmp = randn(N,1)*sig;
        Y = Y + tmp;
        X = X + tmp;
    else
        sigstr = 'Error rate';
    end
    
    % Deconfound using the specified method
    Rz = eye(N) - Z*pinv(Z);
    Rw = eye(N) - W*pinv(W);
    switch deconf_meth
        case {0,4}
            % None
            Qz      = eye(N);
            Qw      = eye(N);
            
        case 1
            % Tom's version (SVD)
            [Qz,S]  = svd(Rz);
            o       = abs(diag(S)) < N*eps(S(1));
            Qz(:,o) = [];
            [Qw,S]  = svd(Rw);
            o       = abs(diag(S)) < N*eps(S(1));
            Qw(:,o) = [];
            
        case 2
            % Huh Juhn's / Anderson's version (Schur)
            [Qz,D]  = schur(Rz);
            o       = abs(diag(D)) < N*eps(D(end:end));
            Qz(:,o) = [];
            [Qw,D]  = schur(Rw);
            o       = abs(diag(D)) < N*eps(D(end:end));
            Qw(:,o) = [];
            
        case 3
            % Theil's BLUS residuals
            Qz      = theil(Z,Rz);
            Qw      = theil(W,Rw);
    end
    NnewY = size(Qz,2);
    NnewX = size(Qw,2);
    
    % Remove nuisance, deconfound as defined above, make sure nothing is rank deficient:
    Y = Qz'*Rz*Y;
    X = Qw'*Rw*X;
    %    Y = Y - mean(Y,1);
    %    X = X - mean(X,1);
    if ~ isempty(Npca) && Npca > 0
        Y = epca(Y,Npca);
        X = epca(X,Npca);
    end
    
    % For each permutation
    for p = 1:nP
        if p == 1
            % First permutation is no permutation
            idxY = (1:NnewY);
            idxX = (1:NnewX);
            
            % Initial CCA, not used for the null distibution.
            % Later we'll use simply these initial U and V, thus
            % dropping Y and X. To make sure U and V span the same
            % spaces as Y and X, include in the loadings the null
            % spaces of A and B into "new" A and B.
            [A0,B0,cc0] = cca(Qz*Y,Qw*X,Nz,abs(Nw));
            %[A0,B0,cc0,~,~] = canoncorr(Y,X);
            if use_nullspace
                nA0 = null(A0');
                nB0 = null(B0');
                U0  = Y*[A0 nA0];
                V0  = X*[B0 nB0];
            else
                U0  = Y*A0;
                V0  = X*B0;
            end
            
            % Since this is for the first perm, initialize
            % variables for later.
            ccs    = zeros(size(cc0));
            lWs    = zeros(size(cc0));
            puncRo = zeros(size(cc0));
            puncWo = zeros(size(cc0));
            puncRs = zeros(size(cc0));
            puncWs = zeros(size(cc0));
            pmaxRo = zeros(size(cc0));
            pmaxWo = zeros(size(cc0));
            pmaxRs = zeros(size(cc0));
            pmaxWs = zeros(size(cc0));
        else
            % If not first permutation, permute randomly
            if deconf_meth == 4
                idxY = randperm(N);
                idxY = idxY(1:(N-size(Z,2)));
                idxX = randperm(N);
                idxX = idxX(1:(N-size(W,2)));
            else
                idxY = randperm(NnewY);
                idxX = randperm(NnewX);
            end
        end
        
        % This is the "one-step" way, all variables
        [~,~,cco,~,~,lWo,ppara] = cca(Qz*U0(idxY,:),Qw*V0(idxX,:),Nz,abs(Nw));
        %[~,~,cco,~,~,stats] = canoncorr(U0(idxY,:),V0(idxX,:));
        %lWo   = -log(stats.Wilks);
        %ppara = stats.pChisq;
        
        % This is the "sequential" way
        % For each canonical component, do the CCA from k onwards
        % This will not affect the first permutation, but will effectively
        % remove from later k the earlier ones. Since these later ones will
        % only be tested if the earlier ones were significant (step-down)
        % this means the earlier ones will act as nuisance, and will not be
        % permuted together.
        for k = 1:numel(cc0)
            [~,~,cctmp,~,~,lWtmp,~] = cca(Qz*U0(idxY,k:end),Qw*V0(idxX,k:end),Nz,abs(Nw));
            ccs(k) = cctmp(1);
            lWs(k) = lWtmp(1);
            %[~,~,cctmp,~,~,statstmp] = canoncorr(U0(idxY,k:end),V0(idxX,k:end));
            %ccs(k) = cctmp(1);
            %lWs(k) = -log(statstmp.Wilks(1));
        end
        
        % Keep the results of the 1st permutation
        if p == 1
            cco1   = cco;
            lWo1   = lWo;
            ccs1   = ccs;
            lWs1   = lWs;
            ppara1 = ppara;
        end
        
        % Increment counters
        puncRo = puncRo + (cco1 <= cco);
        puncWo = puncWo + (lWo1 <= lWo);
        puncRs = puncRs + (ccs1 <= ccs);
        puncWs = puncWs + (lWs1 <= lWs);
        pmaxRo = pmaxRo + (cco1 <= max(cco));
        pmaxWo = pmaxWo + (lWo1 <= max(lWo));
        pmaxRs = pmaxRs + (ccs1 <= max(ccs));
        pmaxWs = pmaxWs + (lWs1 <= max(lWs));
    end
    
    % Compute the p-values
    puncRo = puncRo/nP;
    puncWo = puncWo/nP;
    puncRs = puncRs/nP;
    puncWs = puncWs/nP;
    pmaxRo = pmaxRo/nP;
    pmaxWo = pmaxWo/nP;
    pmaxRs = pmaxRs/nP;
    pmaxWs = pmaxWs/nP;
    
    % Show results per realization
    fprintf('- p-para (parametric, Bartlett):');   fprintf(' %g',ppara1); fprintf('\n');
    fprintf('- p-unc (permutation, Roy,   one):'); fprintf(' %g',puncRo); fprintf('\n');
    fprintf('- p-unc (permutation, Wilks, one):'); fprintf(' %g',puncWo); fprintf('\n');
    fprintf('- p-unc (permutation, Roy,   seq):'); fprintf(' %g',puncRs); fprintf('\n');
    fprintf('- p-unc (permutation, Wilks, seq):'); fprintf(' %g',puncWs); fprintf('\n');
    fprintf('- p-clo (permutation, Roy,   one):'); fprintf(' %g',cummax(puncRo,2)); fprintf('\n');
    fprintf('- p-clo (permutation, Wilks, one):'); fprintf(' %g',cummax(puncWo,2)); fprintf('\n');
    fprintf('- p-clo (permutation, Roy,   seq):'); fprintf(' %g',cummax(puncRs,2)); fprintf('\n');
    fprintf('- p-clo (permutation, Wilks, seq):'); fprintf(' %g',cummax(puncWs,2)); fprintf('\n');
    fprintf('- p-max (permutation, Roy,   one):'); fprintf(' %g',pmaxRo); fprintf('\n');
    fprintf('- p-max (permutation, Wilks, one):'); fprintf(' %g',pmaxWo); fprintf('\n');
    fprintf('- p-max (permutation, Roy,   seq):'); fprintf(' %g',pmaxRs); fprintf('\n');
    fprintf('- p-max (permutation, Wilks, seq):'); fprintf(' %g',pmaxWs); fprintf('\n');
    fprintf('- CCs, one (not permuted):');         fprintf(' %g',cco1);   fprintf('\n');
    fprintf('- CCs, one (random perm):' );         fprintf(' %g',cco);    fprintf('\n');
    fprintf('- CCs, seq (not permuted):');         fprintf(' %g',ccs1);   fprintf('\n');
    fprintf('- CCs, seq (random perm):' );         fprintf(' %g',ccs);    fprintf('\n');
    
    % Store to average over realizations
    ppara_rep(r,:)  = ppara1;
    puncRo_rep(r,:) = puncRo;
    puncWo_rep(r,:) = puncWo;
    puncRs_rep(r,:) = puncRs;
    puncWs_rep(r,:) = puncWs;
    pcloRo_rep(r,:) = cummax(puncRo,2);
    pcloWo_rep(r,:) = cummax(puncWo,2);
    pcloRs_rep(r,:) = cummax(puncRs,2);
    pcloWs_rep(r,:) = cummax(puncWs,2);
    pmaxRo_rep(r,:) = pmaxRo;
    pmaxWo_rep(r,:) = pmaxWo;
    pmaxRs_rep(r,:) = pmaxRs;
    pmaxWs_rep(r,:) = pmaxWs;
    cco1_rep(r,:)   = cco1;
    ccor_rep(r,:)   = cco;
    ccs1_rep(r,:)   = ccs1;
    ccsr_rep(r,:)   = ccs;
end

fprintf('Results:\n')
fprintf('- %s (ppara, Bartlett):',sigstr); fprintf(' %g',mean(ppara_rep  <= alpha)); fprintf('\n');
fprintf('- %s (unc, Roy,   one):',sigstr); fprintf(' %g',mean(puncRo_rep <= alpha)); fprintf('\n');
fprintf('- %s (unc, Wilks, one):',sigstr); fprintf(' %g',mean(puncWo_rep <= alpha)); fprintf('\n');
fprintf('- %s (unc, Roy,   seq):',sigstr); fprintf(' %g',mean(puncRs_rep <= alpha)); fprintf('\n');
fprintf('- %s (unc, Wilks, seq):',sigstr); fprintf(' %g',mean(puncWs_rep <= alpha)); fprintf('\n');
fprintf('- %s (clo, Roy,   one):',sigstr); fprintf(' %g',mean(pcloRo_rep <= alpha)); fprintf('\n');
fprintf('- %s (clo, Wilks, one):',sigstr); fprintf(' %g',mean(pcloWo_rep <= alpha)); fprintf('\n');
fprintf('- %s (clo, Roy,   seq):',sigstr); fprintf(' %g',mean(pcloRs_rep <= alpha)); fprintf('\n');
fprintf('- %s (clo, Wilks, seq):',sigstr); fprintf(' %g',mean(pcloWs_rep <= alpha)); fprintf('\n');
fprintf('- %s (max, Roy,   one):',sigstr); fprintf(' %g',mean(pmaxRo_rep <= alpha)); fprintf('\n');
fprintf('- %s (max, Wilks, one):',sigstr); fprintf(' %g',mean(pmaxWo_rep <= alpha)); fprintf('\n');
fprintf('- %s (max, Roy,   seq):',sigstr); fprintf(' %g',mean(pmaxRs_rep <= alpha)); fprintf('\n');
fprintf('- %s (max, Wilks, seq):',sigstr); fprintf(' %g',mean(pmaxWs_rep <= alpha)); fprintf('\n');
fprintf('- CCs, one (not permuted):');     fprintf(' %g',mean(cco1_rep           )); fprintf('\n');
fprintf('- CCs, one (random perm):');      fprintf(' %g',mean(ccor_rep           )); fprintf('\n');
fprintf('- CCs, seq (not permuted):');     fprintf(' %g',mean(ccs1_rep           )); fprintf('\n');
fprintf('- CCs, seq (random perm):');      fprintf(' %g',mean(ccsr_rep           )); fprintf('\n');

% =================================================================
function [A,B,cc,U,V,lW,ppara] = cca(X,Y,Nz,Nw)
N = size(X,1);
[Qx,Rx,iX] = qr(X,0);
[Qy,Ry,iY] = qr(Y,0);
K  = min(rank(X),rank(Y));
[L,D,M] = svds(Qx'*Qy,K);
cc = min(max(diag(D(:,1:K))',0),1);
A  = Rx\L(:,1:K)*sqrt(N-Nz-1);
B  = Ry\M(:,1:K)*sqrt(N-Nw-1);
A(iX,:) = A;
B(iY,:) = B;
U  = X*A;
V  = Y*B;
lW = -fliplr(cumsum(fliplr(log(1-cc.^2))));
% - - - - - - - - - - - - - - - - - - - - - -
% Parametric p-val, Bartlett (1938)
[~,Nx] = size(X);
[~,Ny] = size(Y);
lW = (N - max(Nz,Nw) - (Nx+Ny+3)/2)*lW;
nu = (Nx-(1:K)+1).*(Ny-(1:K)+1);
ppara = gammaincfix(lW/2,nu/2,'upper');

% =================================================================
function Q = theil(Z,Rz)
% Create (n-k) x n matrix that generates BLUS residuals
% * Magnus JR, Sinha AK. On Theil's errors.
%   The Econometrics Journal. 2005;8(1):39-54.

[n,k] = size(Z);

% Find safe subset to drop
Done = false;
while ~ Done
    i = randperm(n,n-k);
    o = setdiff(1:n,i);
    if rank(Z(o,:)) == k
        Done = true;
    end
end
S = eye(n);
S = S(:,i);
Q = Rz*S*sqrtm(inv(S'*Rz*S));

% =================================================================
function Y = gammaincfix(X,A,tail)
% Fixes an issue with Octave's native gammainc.

% For most cases, this will be enough:
Y = gammainc(X,A,tail);

% This patch is really only for Octave:
if isoctave,
    idx = Y == 0;
    if any(idx),
        gfun = @(t)exp(-t).*t.^(A-1);
        if strcmpi(tail,'upper'),
            Y(idx) = arrayfun(@(x)quad(gfun,x,Inf),X(idx));
        else
            Y(idx) = arrayfun(@(x)quad(gfun,0,x),X(idx));
        end
        Y(idx) = bsxfun(@rdivide,Y(idx),gamma(A));
    end
end

% =================================================================
function y = isoctave
persistent isoct;
if isempty(isoct),
    isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end
y = isoct;

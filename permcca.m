function [pfwer,r,A,B,U,V] = permcca(varargin)
% Permutation inference for canonical correlation
% analysis (CCA).
%
% Inputs:
% - Y        : Left set of variables, size N by P.
% - X        : Left set of variables, size N by Q.
% - nP       : An integer representing the number
%              of permutations.
%              Default is 1000 permutations.
% - Z        : (Optional) Nuisance variables for
%              both (partial CCA) or left side
%              (part CCA) only.
% - W        : (Optional) Nuisance variables for the
%              right side only (bipartial CCA).
% - Sel      : (Optional) Selection matrix.
% - partial  : (Optional) Boolean indicating whether
%              this is partial (true) or part (false) CCA.
%              Default is true, i.e., partial CCA.
%
% Outputs:
% - p   : p-values, FWER corrected via closure.
% - r   : Canonical correlations.
% - A   : Canonical coefficients, left side.
% - B   : Canonical coefficients, right side.
% - U   : Canonical variables, left side.
% - V   : Canonical variables, right side.
%
% ___________________________________________
% AM Winkler, O Renaud, SM Smith, TE Nichols
% NIH - Univ. of Geneva - Univ. of Oxford
% Feb/2020

% Read input arguments
narginchk(2,7)
Y = varargin{1};
Y = center(Y);
X = varargin{2};
X = center(X);
if nargin >= 3
    nP = varargin{3};
end
if nargin >= 4
    Z = varargin{4};
    Z = center(Z);
else
    Z = [];
end
if nargin >= 5
    W = varargin{5};
    W = center(W);
else
    W = [];
end
if nargin >= 6
    Sel = varargin{6};
else
    Sel = [];
end
if nargin >= 7
    partial = varargin{7};
else
    partial = true;
end
[Ny,P] = size(Y);
[Nx,Q] = size(X);
K = min(P,Q);
if Ny ~= Nx
    error('Y and X do not have same number of rows.')
end
N = Ny; clear Ny Nx
I = eye(N);

% Residualise Y wrt Z
if ~ isempty(Z)
    Rz = I - Z*pinv(Z);
    Qz = semiortho(Rz,Sel);
else
    Rz = I;
    Qz = I;
end
Y    = Qz'*Rz*Y;
Pnew = size(Y,1);
R    = size(Z,2);

% Residualise X wrt W
if isempty(W) && partial
    W = Z;
end
if ~ isempty(W)
    if partial
        Rw = Rz;
        Qw = Qz;
    else % bipartial
        Rw = I - W*pinv(W);
        Qw = semiortho(Rw,Sel);
    end
else
    Rw = I;
    Qw = I;
end
X    = Qw'*Rw*X;
Qnew = size(X,1);
S    = size(W,2);

% Initial CCA
[A,B,~] = cca(Y,X,R,S);
U = Y*[A null(A')];
V = X*[B null(B')];

% Initialise counter
cnt = zeros(1,K);
lW  = zeros(1,K);

% For each permutation
for p = 1:nP
    
    % First permutation is no permutation
    if p == 1
        idxY = (1:Pnew);
        idxX = (1:Qnew);
    else
        idxY = randperm(Pnew);
        idxX = randperm(Qnew);
    end
    
    % For each canonical variable
    for k = 1:K
        [~,~,r] = cca(Qz*U(idxY,k:end),Qw*V(idxX,k:end),R,S);
        lWtmp = -fliplr(cumsum(fliplr(log(1-r.^2))));
        lW(k) = lWtmp(1);
    end
    if p == 1
        lW1 = lW;
    end
    
    cnt = cnt + (lW >= lW1);
end
punc  = cnt/nP;
pfwer = cummax(punc);

% =================================================================
function Q = semiortho(R,Sel)
% Compute a semi-orthogonal matrix according to
% the Huh-Jhun or Theil methods.
if isempty(Sel)
    % Huh-Jhun
    [Q,L]  = schur(R);
    o      = abs(diag(L)) < size(R,1)*eps(L(end:end));
    Q(:,o) = [];
else
    % Theil
    Q = R*Sel*sqrtm(inv(Sel'*R*Sel));
end

% =================================================================
function [A,B,cc] = cca(X,Y,R,S)
% Compute CCA.
N = size(X,1);
[Qx,Rx,iX] = qr(X,0);
[Qy,Ry,iY] = qr(Y,0);
K  = min(rank(X),rank(Y));
[L,D,M] = svds(Qx'*Qy,K);
cc = min(max(diag(D(:,1:K))',0),1);
A  = Rx\L(:,1:K)*sqrt(N-R);
B  = Ry\M(:,1:K)*sqrt(N-S);
A(iX,:) = A;
B(iY,:) = B;

% =================================================================
function X = center(X)
% Mean center a matrix and remove constant columns.
icte = sum(diff(X,1,1).^2,1) == 0;
X = bsxfun(@minus,X,mean(X,1));
X(:,icte) = [];

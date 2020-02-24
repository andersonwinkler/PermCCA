function [lEvec,rEvec,Evals,xtras] = epca(varargin)
% Return the first p eigenvectors and eigenvalues.
%
% Usage:
% [lEvec,rEvec,Evals,xtras] = epca(X,p,dim,norm)
%
% Inputs:
% X     : 2D array
% p     : Number of eigenvectors and eigenvalues to be
%         returned.
% dim   : Dimension used for mean-centering.
%         Valid inputs are 1 (or 'rows') or 2 (or 'cols').
%         If omitted, no mean-centering is performed.
% norm  : True/False, indicationg whether, in addition to
%         mean-centering, the data should also have the
%         variance normalised to unity.
%
% Outputs:
% lEvec : Left singular vectors (= loadings for the rEvec)
% rEvec : Right singular vectors (= loadings for the lEvec).
% Evals : Eigenvalues.
% xtras : A struct containing the eigenvectors scaled
%         by their eigenvalues, and the data recovered
%         after the p+1...rank(X) dimensions were removed.
%
% This function was called "pca", but after Matlab created its
% own similarly-named "pca" function, it was renamed to "epca",
% for "enhanced pca".
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Nov/2010 (first version)
% May/2015 (this version)
% http://brainder.org

% Accept and check arguments
p  = 1;
dm = false;
nm = false;
if nargin < 1 || nargin > 4,
    error('Incorrect number of arguments');
elseif nargin == 2,
    p  = varargin{2};
elseif nargin == 3,
    p  = varargin{2};
    dm = varargin{3};
elseif nargin == 4,
    p  = varargin{2};
    dm = varargin{3};
    nm = varargin{4};
end
X = varargin{1};
if ischar(dm),
    if strcmpi(dm,'rows'),
        dm = 1;
    elseif strcmpi(dm,'cols'),
        dm = 2;
    else
        error('Unknown option "%s"',dm);
    end
end
[nR,nC] = size(X);
if p > max(nR,nC),
    error('Cannot extract more eigenvalues than rows or columns.');
end

% Demean and/or normalise along columns or rows
if dm,
    Xm = mean(X,dm);
    X  = bsxfun(@minus,X,Xm);
    if nm,
        Xs = std(X,[],dm);
        X  = bsxfun(@rdivide,X,Xs);
    end
end

% Save some memory by working with the
% smallest possible square of X
if nR >= nC,
    [~,SS,V] = svd(X'*X);
    Vp  = V(:,1:p);
    SSp = SS(1:p,1:p);
    Up  = X*Vp;
else
    [U,SS,~] = svd(X*X');
    Up  = U(:,1:p);
    SSp = SS(1:p,1:p);
    Vp  = X'*Up;
end

% Pick one sign
s = diag(sign(Up(1,:)));

% Eigenvectors and eigenvalues
lEvec = Up*s;
rEvec = Vp*s;
Evals = diag(SSp)./(nR-1);

% Some extra outputs, not normally needed
if nargout == 4,
    
    % Scaled eigenvectors
    xtras.lScal = X*Vp;
    xtras.rScal = Up'*X;
    
    % Unit norm eigenvectors
    if nR >= nC,
        xtras.lEvec1 = lEvec/sqrt(SSp);
        xtras.rEvec1 = rEvec;
    else
        xtras.lEvec1 = lEvec;
        xtras.rEvec1 = (rEvec/sqrt(SSp))';
    end
    
    % Recovered data using the p eigenvectors
    xtras.lRec = Up*xtras.rScal;
    xtras.rRec = xtras.lScal*Vp';
end

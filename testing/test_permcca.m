% This is just a quick test for sanity of "permcca.m".
% For the main simulation function, see "ccainf.m".
N = 40;
Ny = 8;
Nx = 12;
Nz = 15;
nR = 100;
nP = 100;
pfwer = zeros(nR,8);
for r = 1:nR,
    fprintf('%d\n',r);
    Y = randn(N,Ny);
    X = randn(N,Nx);
    Z = randn(N,Nz);
    pfwer(r,:) = permcca(Y,X,nP,Z,Z);
end
fwer = mean(pfwer<=.05,1);


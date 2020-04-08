% This is just a quick test for sanity of "permcca.m".
% For the main simulation function, see "ccainf.m".
N = 50;
Ny = 8;
Nx = 12;
Nz = 15;
Nw = 3;
nR = 200;
nP = 200;
pfwer = zeros(nR,8);
for rep = 1:nR,
    fprintf('%d\n',rep);
    Y = randn(N,Ny);
    X = randn(N,Nx);
    Z = randn(N,Nz);
    W = randn(N,Nw);
    [pfwer(rep,:),r,A,B,U,V] = permcca(Y,X,nP,Z,Z);
end
fwer = mean(pfwer<=.05,1);


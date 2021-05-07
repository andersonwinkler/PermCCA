rng(0)
N = 100;
Ny = 15;
Nx = 12;
Nz = 10;
Nw = 10;
nP = 100;

Y = randn(N,Ny);
X = randn(N,Nx);
Z = randn(N,Nz);

%[pfwer,r,A,B,U,V] = permcca(Y,X,nP,Z);

[pfwer,r,A,B,U,V,pA,pB] = permloads(Y,X,nP,Z,[],[],true,false,10);





rng(1)
N = 100;
Ny = 15;
Nx = 12;
Nz = 10;
Nw = 10;
nP = 100;
nRlz = 500;

pAr = zeros(Ny,min(Ny,Nx),nRlz);
pBr = zeros(Nx,min(Ny,Nx),nRlz);
for rlz = 1:nRlz
Y = randn(N,Ny);
X = randn(N,Nx);
Z = randn(N,Nz);

[pfwer,r,A,B,U,V,pA,pB] = permloads(Y,X,nP,Z,[],[],true,false);
pAr(:,:,rlz) = pA;
pBr(:,:,rlz) = pB;
end





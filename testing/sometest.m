rng(0)
N = 100;
Ny = 28;
Nx = 21;
Y = randn(N,Ny)*10;
X = randn(N,Nx)*5;
Y = Y-mean(Y);
X = X-mean(X);

Y = bsxfun(@rdivide,Y,std(Y));
X = bsxfun(@rdivide,X,std(X));


[A,B,cc] = cca(Y,X,1,1);
U = Y*A;
V = X*B;

lA = corr(Y,U);
lB = corr(X,V);


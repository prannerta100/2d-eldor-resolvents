figure;
rng('shuffle');
sz = 100;
bnd = 4;
a = diag(randn(sz-bnd,1),-bnd) + diag(randn(sz,1)) + diag(randn(sz-bnd,1), bnd);
invA = inv(a);
imagesc(log10(abs(invA(:,:))));colorbar;
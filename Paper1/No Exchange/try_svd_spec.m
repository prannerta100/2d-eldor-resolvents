matxr = spconvert(load('matlab_realh_34050_9.mtxx'));
matxi = spconvert(load('matlab_imagh_34050_9.mtxx'));
shiftroffdiag = 1.0;
matx = matxr + 1i* matxi + shiftroffdiag * speye(size(matxr,1));
omarr = linspace(-333,333,128);
cfact = 17.6372;
omarrG = omarr/(cfact/(2*pi));
spec1 = zeros(128,1);
spec2 = zeros(128,1);
stvx = [1/sqrt(3);1/sqrt(3);1/sqrt(3);zeros(size(matx,1)-3,1)];
tic
parfor i=1:128
    %[u,s,v] = svds(matx+1i*omarrG(i)*speye(size(matx,1)),100,'smallest');
    [u,s,v] = svd(full(matx+1i*omarrG(i)*speye(size(matx,1))));
    spec1(i) = (stvx'*v)*(s\(u'*stvx));
    spec2(i) = stvx'*((matx+1i*omarrG(i)*speye(size(matx,1)))\stvx);
end
toc
%%
plot(omarrG(1:end-1),diff(real(spec1)));
hold on;
plot(omarrG(1:end-1),diff(real(spec2)));
legend('SVD','spLU');
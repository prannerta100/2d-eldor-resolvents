matxr = spconvert(load('matlab_realh_34050_5.mtxx'));
matxi = spconvert(load('matlab_imagh_34050_5.mtxx'));
stvx = load('stvech_34050_5.stvx');
stvx = stvx(:,1) + 1i* stvx(:,2);
size(stvx)
matx = matxr + 1i* matxi;

[Tcsym, ~] = myLanczosCSym(matx, 1000, stvx);
plot(abs(normlist));
%%
[Therm, ~] = myLanczosHerm(matx, 1000, stvx);
figure;
plot(abs(diag(Therm,1)));
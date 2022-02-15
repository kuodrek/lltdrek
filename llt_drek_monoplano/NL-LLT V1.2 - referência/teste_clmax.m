clear
clc
b = [0.630 0.431 0.132];
c = [0.525 0.525 0.26 0.13];
particoes = 3;
diedro = [0 0 0];
offset = [0 0.0662 0.0988];
twist = 0;
h = 0.1;
nperfil = [2 2 2 2];
Vinf  = 15;
[~,~,~,~,~,~,~,~,~,~,asa_geometria,airfoil_data] = llt_principal(particoes,b,c,diedro,offset,nperfil,twist,h,Vinf);

tic()
[CLmax,alfa_estol,Cl] = llt_CLmax(asa_geometria,airfoil_data,Vinf)
toc()
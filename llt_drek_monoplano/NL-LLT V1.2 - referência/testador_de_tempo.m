% Testador de tempo
clear
clc

particoes = 1;
b = 1;
c = [0.4 0.4];
diedro = 0;
offset = 0;
nperfil = [3 3];
twist = 0;
h = 0;
Vinf = 15;
damping = linspace(0.02,0.1,10);
for i=1:size(damping,2)
    damping(i)
[Sw,ARw,CLalphaw,CL0w,MACw,bw,CDi_coeff,Cmacw,xACw,Cl] = llt_principal(particoes,b,c,diedro,offset,nperfil,twist,h,Vinf,damping(i));
end

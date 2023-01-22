clear
clc

b = 4.572/2;
c = [0.813 0.325];
diedro = [0 0];
nperfil = [1 2];
offset = 0.122;
particoes = 1;
twist = -2.4;
h = 0;


aoa = [-5.9
-4.9
-3.8
-2.8
-1.7
-0.7
0.4
1.5
2.4
3.5
4.5
5.6
6.6
7.6
8.7
9.7
10.7
11.8
12.7
13.8
14.8
15.9
16.8
17.84
18.8
19.8
20.8];

vetor_CL = zeros(size(aoa,1),2);

rho = 2.83;
mi = 1.81e-5;
nu = 6.39e-6;
Vinf = 45.48;

% Inicialização dos perfis
[wing_airfoil_set,tail_airfoil_set] = llt_initialize_airfoils;
airfoil_set = wing_airfoil_set;
AC_check = 'n';
[asa_geometria] = AeroGeometry(particoes,b,c,diedro,offset,nperfil,twist,airfoil_set,AC_check,Vinf);

% Fatores da simulação
damping = 0.05;
iteracoes = 1000;
precisao = 1e-3;
N = asa_geometria{8,1};
G_inicial = ones(N,1)*0.1;

airfoil_data = asa_geometria{9,1};
linear_check = 'nlinear';
referencial = [0 0 0];

tic()
for i=1:size(aoa,1)
[G] = llt_main_solver(asa_geometria,airfoil_data,Vinf,h,aoa(i),damping,iteracoes,precisao,G_inicial,linear_check);
if isnan(G) == 0
        G_inicial = G;
end
[CF,CM,Cl] = llt_coeficientes(asa_geometria,Vinf,referencial,G,h,aoa(i));
vetor_CL(i,2) = CF(3);
vetor_CL(i,1) = aoa(i);
end
toc()

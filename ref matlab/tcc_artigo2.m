clear
clc

b = 4.572/2;
c = [0.726 0.290];
diedro = [0 3];
nperfil = [5 5];
offset = 0.109;
particoes = 1;
twist = -2;
h = 0;


aoa = [-3.0
-2.0
-0.8
0.3
1.3
2.4
3.4
4.4
5.5
6.5
7.6
8.8
9.7
10.8
11.3
11.9
12.4
12.9
13.5
14.0
14.4
15.1
15.4];

vetor_CL = zeros(size(aoa,1),2);

rho = 2.43;
mi = 1.81e-5;
nu = 7.44e-6;
Vinf = 57.85;

% Inicialização dos perfis
[wing_airfoil_set,tail_airfoil_set] = llt_initialize_airfoils;
airfoil_set = wing_airfoil_set;
AC_check = 'n';
[asa_geometria] = AeroGeometry(particoes,b,c,diedro,offset,nperfil,twist,airfoil_set,AC_check,Vinf,nu);

% Fatores da simulação
damping = 0.8;
iteracoes = 1000;
precisao = 1e-3;
N = asa_geometria{8,1};
G_inicial = ones(N,1);

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

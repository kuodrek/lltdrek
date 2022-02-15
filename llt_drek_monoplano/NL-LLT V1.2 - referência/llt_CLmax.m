function [CLmax,alfa_estol,Cl] = llt_CLmax(asa_geometria,airfoil_data,Vinf)

% Extração das variáveis

b = asa_geometria{1,1};
c = asa_geometria{2,1};
offset = asa_geometria{3,1};
diedro = asa_geometria{4,1};
twist = asa_geometria{5,1};
nperfil = asa_geometria{6,1};
particoes = asa_geometria{7,1};
h = asa_geometria{8,1};
N = asa_geometria{9,1};

collocation_points = asa_geometria{1,2};
vertices = asa_geometria{2,2};
csi = asa_geometria{3,2};
area = asa_geometria{4,2};
mac = asa_geometria{5,2};
u_n = asa_geometria{6,2};
u_a = asa_geometria{7,2};
reynolds_info = asa_geometria{8,2};
airfoil_info = asa_geometria{9,2};
MACw = asa_geometria{10,2};
referencial = zeros(N,3);

CL_conv = 0;
aoa_inicial = 10;
aoa_passo = [5 2.5 1 0.5];
cont = 0;
G_inicial = ones(N,1); %Primeiro chute
%Início da solução não linear
damping = 0.05;                         %Coeficiente de amortecimento entre iterações
iteracoes = 300;
precisao = 0.001;
aux_passo = 1;
while CL_conv == 0
    cont = cont + 1
    if cont == 1
        aoa = aoa_inicial;
        CL_anterior = 0;
    end
    
    aoa_rad = aoa*pi/180;           %Transformação para radianos
    Vinf_vetor = Vinf*[cos(aoa_rad) 0 sin(aoa_rad)];
    vinf = Vinf_vetor/Vinf;           %Vetor unitário da corrente livre
    
    [G] = llt_main_solver(asa_geometria,airfoil_data,Vinf,aoa,damping,iteracoes,precisao,G_inicial);
    if isnan(G(1)) == 0
        [CF,~,Cl_aux] = llt_coeficientes(vertices,collocation_points,csi,b,nperfil,vinf,area,mac,reynolds_info,airfoil_info,airfoil_data,MACw,u_n,u_a,referencial,G,h,N,aoa);
        aoa
        CL = CF(3)
        if abs(CL-CL_anterior) < 0.05
            CL_conv = 1;
            CLmax = CL;
            alfa_estol = aoa;
            Cl = Cl_aux;
        end
        if CL > CL_anterior
            aoa = aoa + aoa_passo(aux_passo);
            CL_anterior = CL;
        else
            aoa = aoa - aoa_passo(aux_passo);
            if aux_passo < size(aoa_passo,2)
                aux_passo = aux_passo + 1;
                aoa = aoa + aoa_passo(aux_passo);
            end
        end
    else
        aoa = aoa - aoa_passo(aux_passo);
            if aux_passo < size(aoa_passo,2)
                aux_passo = aux_passo + 1;
                aoa = aoa + aoa_passo(aux_passo);
            end
    end
    if cont > 25
        CL_conv = 1;
        CLmax = 0;
        alfa_estol = 0;
        Cl = zeros(N,1);
        fprintf('Resultado não convergiu');
    end
    
end
end
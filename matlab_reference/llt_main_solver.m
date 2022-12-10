function [G] = llt_main_solver(asa_geometria,airfoil_data,Vinf,h,aoa,damping,iteracoes,precisao,G_inicial,linear_check)

    b = asa_geometria{1,1};
    nperfil = asa_geometria{6,1};
    N = asa_geometria{8,1};
    
    collocation_points = asa_geometria{1,2};
    vertices = asa_geometria{2,2};
    csi = asa_geometria{3,2};
    mac = asa_geometria{5,2};
    u_n = asa_geometria{6,2};
    u_a = asa_geometria{7,2};
    reynolds_info = asa_geometria{8,2};
    airfoil_info = asa_geometria{9,2};
    
    aoa_rad = aoa*pi/180;           %Transformação para radianos
    Vinf_vetor = Vinf*[cos(aoa_rad) 0 sin(aoa_rad)];
    vinf = Vinf_vetor/Vinf;           %Vetor unitário da corrente livre
    
    convergencia = 1;                      %Variável lógica de convergência
    cont = 1;                        %Contador pra verificar limite máximo de tentativas
    G = G_inicial;                  %Chute inicial do conjunto de circulações adimensionalizadas
    while convergencia == 1
        [R,aoa_eff] = llt_solver1(vertices,collocation_points,csi,b,nperfil,vinf,u_n,u_a,mac,reynolds_info,airfoil_info,airfoil_data,G,h,N,linear_check);
        deltaG = llt_solver2(vertices,collocation_points,b,nperfil,csi,vinf,u_n,u_a,mac,reynolds_info,airfoil_info,airfoil_data,aoa_eff,G,-R,h,N,linear_check);
        if cont > iteracoes
            convergencia = 2;
            fprintf("aoa: %g\niteracao: %g\n", aoa,cont)
            fprintf('resultado não convergiu\n');
            G = ones(N,1)*NaN;
        end
        
        if abs(max(R)) < precisao
             fprintf("aoa: %g\niteracao: %g\n", aoa,cont)
            convergencia = 2;
        else
            cont = cont + 1;
            G = G + damping*deltaG;
        end
    end
    end
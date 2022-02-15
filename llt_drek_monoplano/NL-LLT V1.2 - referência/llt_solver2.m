function [deltaG] = llt_solver2(vertices,collocation_points,b,nperfil,csi,vinf,u_n,u_a,mac,reynolds_info,airfoil_info,airfoil_data,aoa_eff,G,R,h,N)
% V1.2

J = zeros(2*N, 2*N);
R_aux = zeros(2*N, 1);

u_n_oposto = zeros(N, 3);
u_a_oposto = zeros(N, 3);
csi_oposto = zeros(N, 3);

for i=1:N
    u_n_oposto(i, :) = [u_n(i, 1) -u_n(i, 2) u_n(i, 3)];
    u_a_oposto(i, :) = [u_a(i, 1) -u_a(i, 2) u_a(i, 3)];
    csi_oposto(i, :) = [-csi(i, 1) csi(i, 2) -csi(i, 3)];
end

%i: variação dos pontos de colocação
%j: variação dos vórtices

for i=1:N
    soma = zeros(1, 3);
    
    [vji,vji_sim] = llt_vind(vertices, collocation_points, vinf, mac, h, N, i);
    for j=1:N
        soma = soma + (vji(j, :)+vji_sim(j, :))*G(j);
    end
    
    w = cross(vinf+soma,csi(i,:));
    vn = dot(vinf+soma,u_n(i, :));
    va = dot(vinf+soma, u_a(i, :));
    
    w_oposto = [w(1, 1) -w(1, 2) w(1, 3)];
    
    cl_alfa = data_solver(reynolds_info,airfoil_info,airfoil_data,b,nperfil,aoa_eff(i),i,'cl_alfa');
    
    for j=1:N
        aux_J1 = cross(vji(j, :), csi(i, :));
        aux_J2 = dot(vji(j, :), u_n(i, :));
        aux_J3 = dot(vji(j, :), u_a(i, :));
        
        aux_J1_sim = cross(vji_sim(j, :), csi(i, :));
        aux_J2_sim = dot(vji_sim(j, :), u_n(i, :));
        aux_J3_sim = dot(vji_sim(j, :), u_a(i, :));
        
        aux_J1_op_sim = [aux_J1_sim(1, 1) -aux_J1_sim(1, 2) aux_J1_sim(1, 3)];
        aux_J1_op = [aux_J1(1, 1) -aux_J1(1, 2) aux_J1(1, 3)];
        
        if i == j
            %Painel original, vórtice original
            J(N+i, N+j) =   2*norm(w) + (1/norm(w))*dot(2*w, aux_J1)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2 - vn*aux_J3);
            %Painel oposto, vórtice oposto
            J(N+1-i, N+1-j) =   2*norm(w_oposto) + (1/norm(w_oposto))*dot(2*w_oposto, aux_J1_op)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2 - vn*aux_J3);
        else
            %Painel original, vórtice original
            J(N+i, N+j) =   (1/norm(w))*dot(2*w, aux_J1)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2 - vn*aux_J3) ;
            %Painel original, vórtice oposto
            J(N+i, N+1-j) =   (1/norm(w))*dot(2*w, aux_J1_sim)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2_sim - vn*aux_J3_sim);
            
            %Painel oposto, vórtice original
            J(N+1-i, N+j) =   (1/norm(w_oposto))*dot(2*w_oposto, aux_J1_op_sim)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2_sim - vn*aux_J3_sim);
            %Painel oposto, vórtice oposto
            J(N+1-i, N+1-j) =   (1/norm(w_oposto))*dot(2*w_oposto, aux_J1_op)*G(i) - cl_alfa/(vn^2+va^2)*(va*aux_J2 - vn*aux_J3);
        end
    end
end
for i=1:N
    R_aux(N+1-i, 1) = R(i, 1);
    R_aux(N+i, 1) = R(i, 1);
end
deltaG_aux = linsolve(J, R_aux);
deltaG = deltaG_aux((N+1):2*N);
end
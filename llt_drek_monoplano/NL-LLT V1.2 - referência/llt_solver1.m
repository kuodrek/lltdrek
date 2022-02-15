function [R,aoa_eff] = llt_solver1(vertices,collocation_points,csi,b,nperfil,vinf,u_n,u_a,mac,reynolds_info,airfoil_info,airfoil_data,G,h,N)
%V1.2
R = zeros(N, 1);

%i: variação dos pontos de colocação
%j: variação dos vórtices
aoa_eff = zeros(N,1);
for i=1:N
    soma = zeros(1, 3);
    
    [vji,vji_sim] = llt_vind(vertices, collocation_points, vinf, mac, h, N, i);
    for j=1:N
        soma = soma + (vji(j,:)+vji_sim(j,:))*G(j);
    end
    aoa_i = atan(dot(soma+vinf, u_n(i,:))/dot(soma+vinf,u_a(i,:)));
    aoa_eff(i) = aoa_i;
    [cl_eff] = data_solver(reynolds_info,airfoil_info,airfoil_data,b,nperfil,aoa_i,i,'cl');
    R(i, 1) = 2*norm(cross(soma+vinf, csi(i,:)))*G(i) - cl_eff;
end
end
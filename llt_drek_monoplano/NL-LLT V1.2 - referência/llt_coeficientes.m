function [F,M,Cl] = llt_coeficientes(vertices,collocation_points,csi,b,nperfil,vinf,area,mac,reynolds_info,airfoil_info,airfoil_data,MACw,u_n,u_a,referencial,G,h,N,aoa)
%Nl LLT - v1.2
u_n_oposto = zeros(N, 3);
u_a_oposto = zeros(N, 3);
csi_oposto = zeros(N, 3);
for i=1:N
    u_n_oposto(i, :) = [u_n(i, 1) -u_n(i, 2) u_n(i, 3)];
    u_a_oposto(i, :) = [u_a(i, 1) -u_a(i, 2) u_a(i, 3)];
    csi_oposto(i, :) = [-csi(i, 1) csi(i, 2) -csi(i, 3)];
end

u_s = cross(u_a(:, :), u_n(:, :)); %Vetor pra fazer o coeficiente de momentos
u_s_oposto = cross(u_a_oposto(:, :), u_n_oposto(:, :));
r = collocation_points - referencial;
r_oposto = [collocation_points(:, 1) -collocation_points(:, 2) collocation_points(:, 3)] - referencial;

soma_forcas = zeros(N, 3);
soma_forcas_oposto = zeros(N, 3);
soma_momentos = zeros(N, 3);
soma_momentos_oposto = zeros(N, 3);

Cl = zeros(N,1);
%i: variação dos pontos de colocação
%j: variação dos vórtices

%Ordem das velocidades induzidas: 
%Painel original, vórtice original
%Painel original, vórtice oposto
%Painel oposto, vórtice original
%Painel oposto, vórtice oposto
for i=1:N
    soma = zeros(1, 3);
    soma_oposto = zeros(1, 3);
    
    [vji,vji_sim] = llt_vind(vertices, collocation_points, vinf, mac, h, N, i);
    vji_op = [vji(:,1) -vji(:,2) vji(:,3)];
    vji_sim_op = [vji_sim(:,1) -vji_sim(:,2) vji_sim(:,3)];
    
    for j=1:N
        soma = soma + (vji(j, :)+vji_sim(j,:))*G(j);
        soma_oposto = soma_oposto + (vji_op(j, :)+vji_sim_op(j,:))*G(j);
    end
   
    Cl(i) = 2*G(i);
    
    soma = (soma + vinf)*G(i);
    soma_oposto = (soma_oposto + vinf)*G(i);
    
    cm = data_solver(reynolds_info,airfoil_info,airfoil_data,b,nperfil,aoa,i,'cm');
    
    soma_momentos(i, :) = ( cross(2*r(i, :), cross(soma, csi(i, :))) - cm*mac(i)*u_s(i, :))*area(i, 1)/(sum(2*area)*MACw);
    soma_momentos_oposto(i, :) = ( cross(2*r_oposto(i, :), cross(soma_oposto, csi_oposto(i, :))) - cm*mac(i)*u_s_oposto(i, :))*area(i, 1)/(sum(2*area)*MACw);
    soma_forcas(i, :) = 2*cross(soma, csi(i, :))*area(i, 1)/sum(2*area);
    soma_forcas_oposto(i, :) = 2*cross(soma_oposto, csi_oposto(i, :))*area(i, 1)/sum(2*area);
end
M = sum(soma_momentos)+sum(soma_momentos_oposto);
aoa_rad = aoa*pi/180;
F = [cos(aoa_rad) 0 sin(aoa_rad); 0 1 0; -sin(aoa_rad) 0 cos(aoa_rad)]*(sum(soma_forcas) + sum(soma_forcas_oposto))';
end
function [G] = llt_simplificado(vertices,collocation_points,b,nperfil,csi,vinf,u_n,mac,reynolds_info,airfoil_info,airfoil_data,h,N)
%V1.2

A = zeros(2*N, 2*N);
B = zeros(2*N, 1);

u_n_oposto = zeros(N, 3);
csi_oposto = zeros(N, 3);

for i=1:N
    u_n_oposto(i, :) = [u_n(i, 1) -u_n(i, 2) u_n(i, 3)];
    csi_oposto(i, :) = [-csi(i, 1) csi(i, 2) -csi(i, 3)];
end

%i: variação dos pontos de colocação
%j: variação dos vórtices
for i=1:N
    a0 = 0;
    a1 = 10;
    [cl0] = data_solver(reynolds_info,airfoil_info,airfoil_data,b,nperfil,a0*pi/180,i,'cl');
    [cl1] = data_solver(reynolds_info,airfoil_info,airfoil_data,b,nperfil,a1*pi/180,i,'cl');
    cl_alpha = (cl1-cl0)/(a1-a0)*pi/180;
    al_0 = -cl0/cl_alpha;
    A(N+i, N+i) = 2*norm(cross(vinf, csi(i, :)));
    A(N+1-i, N+1-i) = 2*norm(cross(vinf, csi_oposto(i, :)));
    B(N+i, 1) = cl_alpha*(dot(vinf, u_n(i, :))-al_0);
    B(N+1-i, 1) = cl_alpha*(dot(vinf, u_n_oposto(i, :))-al_0);
    
   [vji,vji_sim] = llt_vind(vertices, collocation_points, vinf, mac, h, N, i);
   vji_op = [vji(:,1) -vji(:,2) vji(:,3)];
   vji_sim_op = [vji_sim(:,1) -vji_sim(:,2) vji_sim(:,3)];
   for j=1:N
        A(N+i, N+j) = A(N+i, N+j) -cl_alpha*(dot(vji(j,:),u_n(i, :)));
        A(N+i, N+1-j) =  -cl_alpha*(dot(vji_sim(j,:),u_n(i, :)));
        
        A(N+1-i, N+j) =  -cl_alpha*(dot(vji_op(j,:),u_n_oposto(i, :)));
        A(N+1-i, N+1-j) = A(N+1-i, N+1-j) -cl_alpha*(dot(vji_sim_op(j,:),u_n_oposto(i, :)));
    end
end
X = linsolve(A, B);
G = X( (N+1) :2*N, 1);
end
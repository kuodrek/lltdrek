function [vertices, collocation_points, csi, area, MACw, mac, AR, u_n, u_a, reynolds_info, airfoil_info, cordas] = llt_geometria(b_particoes,c,offset,tw,diedro,sweep,nperfil,Vinf,nu,particoes,N)
%{
NL LLT
V1.2
%}
b_particoes_backup = b_particoes;
b_particoes = b_particoes*2;                                               %O input é de meia-envergadura porém as contas são realizadas como envergaduras inteiras
b_real = b_particoes;
b = sum(b_particoes);
b_particoes = b_particoes.*cosd(diedro);
part = zeros(1, particoes);                                                %Vetor que contém o número de painéis por partição
S = zeros(1, particoes);                                                   %Área de cada partição da asa
if particoes > 1
    for i=1:(particoes-1)
        %Distribuição do número de paineis por partição de asa.
        num_part = (b_particoes(i)/b)*N;                                   %Número de paineis em função do comprimento relativo da partição
        part(i) = floor(num_part);
    end
    part(end) = N - sum(part);                                             %A última partição recebe o número de painéis restantes
else
    part(end) = N;                                                         %Caso onde existe somente uma partição
end
collocation_points = zeros(N, 3);                                          %Pontos onde serão realizados os cálculos
vertices = zeros(N+1, 3);                                                  %Início e fim de cada painel

XYZ = zeros(particoes,3);
for i=1:particoes
    if i==1
        XYZ(i,1) = 0.25*c(i+1) + offset(i);
        XYZ(i,2) = b_particoes(i)/2;
        XYZ(i,3) = tan(diedro(i)*pi/180)*b_real(i)/2;
    else
        XYZ(i,1) = offset(i)+0.25*c(i+1);
        XYZ(i,2) = XYZ(i-1,2)+b_particoes(i)/2;
        XYZ(i,3) = XYZ(i-1,3)+tan(diedro(i)*pi/180)*b_real(i)/2;
    end
end
soma_part=0;
for i=1:particoes
    if i==1
        X12 = XYZ(i,:) - [0.25*c(i) 0 0];
        X1 = [0.25*c(i) 0 0];
    else
        X12 = XYZ(i,:)-XYZ(i-1,:);
        X1 = XYZ(i-1,:);
    end
    aspace1=0:pi/(part(i)):pi;
    aspace2=pi/(2*(part(i))):pi/(part(i)):pi-pi/(2*(part(i)));
    Lspace=(1-cos(aspace1))/2;
    Mspace=(1-cos(aspace2))/2;
    %     Lspace=(1-sin(aspace1))/2;
    %     Mspace=(1-sin(aspace2))/2;
    X = zeros(part(i)+1,3);
    Y = zeros(part(i),3);
    X(:,1)=X1(1)+Lspace.*X12(1);
    X(:,2)=X1(2)+Lspace.*X12(2);
    X(:,3)=X1(3)+Lspace.*X12(3);
    Y(:,1)=X1(1)+Mspace.*X12(1);
    Y(:,2)=X1(2)+Mspace.*X12(2);
    Y(:,3)=X1(3)+Mspace.*X12(3);
    comeco=soma_part+1;
    fim=soma_part + part(i);
    collocation_points(comeco:fim,1) = Y(:,1);
    collocation_points(comeco:fim,2) = Y(:,2);
    collocation_points(comeco:fim,3) = Y(:,3);
    vertices(comeco:fim+1,1) = X(:,1);
    vertices(comeco:fim+1,2) = X(:,2);
    vertices(comeco:fim+1,3) = X(:,3);
    
    soma_part = soma_part + part(i);
end
% for i=1:N
%     collocation_points(i,2) = (vertices(i+1,2) + vertices(i,2))/2;
% end




mac = zeros(N, 1);
u_n = zeros(N, 3);                                                         %Vetor unitário normal à corda de cada painel
u_a = zeros(N, 3);                                                         %Vetor colinear à corda de cada painel
area = zeros(N, 1);                                                        %Área de cada painel [m^2]
length = zeros(N, 3);                                                      %Vetor de comprimento de cada painel da asa [m]
csi = zeros(N, 3);                                                         %Coeficiente adimensional que relaciona a geometria e vorticidade causada por cada painel

twist_rad = zeros(N, 1);
sweep_rad = zeros(N, 1);
diedro_rad = zeros(N, 1);

j = 1;
soma_part = 0;
for i=1:N
    if i <= (part(j)+soma_part)
        sweep_rad(i) = sweep(j)*pi/180;
        diedro_rad(i) = -diedro(j)*pi/180;
        if j == particoes                                                  %Caso o for esteja na última partição, aplica-se o twist geométrico
            twist_rad(i) = -((pi/180)*(tw-0))/(0.5*b_particoes(j)) * (collocation_points(i, 2)-vertices(soma_part+1,2));
        end
    else
        soma_part = soma_part + part(j);
        j = j + 1;
        sweep_rad(i) = sweep(j)*pi/180;
        diedro_rad(i) = -diedro(j)*pi/180;
    end
end
j = 1;
soma_part = 0;
reynolds_info = zeros(N,1);
airfoil_info = zeros(N,1);
cordas = zeros(N,1);
for i = 1:N                                                                %For para várias geometrias dos paineis
    %Estas matrizes de rotação utilizam os ângulos de Euler
    %Convenção: Z Y X
    u_n(i, :) = [1 0 0; 0 cos(diedro_rad(i)) sin(diedro_rad(i)); 0 -sin(diedro_rad(i)) cos(diedro_rad(i))]...
        *[cos(twist_rad(i)) 0 -sin(twist_rad(i));  0 1 0; sin(twist_rad(i)) 0 cos(twist_rad(i))]...
        *[cos(sweep_rad(i)) sin(sweep_rad(i)) 0; -sin(sweep_rad(i)) cos(sweep_rad(i)) 0; 0 0 1]*[0; 0; 1];
    
    u_a(i, :) = [1 0 0; 0 cos(diedro_rad(i)) sin(diedro_rad(i)); 0 -sin(diedro_rad(i)) cos(diedro_rad(i))]...
        *[cos(twist_rad(i)) 0 -sin(twist_rad(i));  0 1 0; sin(twist_rad(i)) 0 cos(twist_rad(i))]...
        *[cos(sweep_rad(i)) sin(sweep_rad(i)) 0; -sin(sweep_rad(i)) cos(sweep_rad(i)) 0; 0 0 1]*[1; 0; 0];
    
    x = vertices(i+1, 1) - vertices(i, 1);                                 %x, y e z são coordenadas do vetor que representa o comprimento de linha de cada seção
    y = vertices(i+1, 2) - vertices(i, 2);
    z = vertices(i+1, 3) - vertices(i, 3);
    length(i, :) = [x y z];
    
    if i <= (part(j)+soma_part)
        % c1 (mais perto da raíz) e c2 (mais perto da ponta) são as cordas de cada seção
        c1 = ( c(j+1)-c(j) )/(0.5*b_particoes(j)) * ( vertices(i, 2)-vertices(soma_part+1,2) ) + c(j);
        c2 = ( c(j+1)-c(j) )/(0.5*b_particoes(j)) * ( vertices(i+1, 2)-vertices(soma_part+1,2) ) + c(j);
        mac(i) = (2/3)*(c1^2 + c1*c2 + c2^2)/(c1+c2);
        area(i, 1) = 0.5*(c1+c2)*sqrt(y^2+z^2);
        S(j) = S(j) + area(i,1);
    else
        soma_part = soma_part + part(j);
        j=j+1;
        c1 = ( c(j+1)-c(j) )/(0.5*b_particoes(j)) * ( vertices(i, 2)-vertices(soma_part+1,2) ) + c(j);
        c2 = ( c(j+1)-c(j) )/(0.5*b_particoes(j)) * ( vertices(i+1, 2)-vertices(soma_part+1,2) ) + c(j);
        mac(i) = (2/3)*(c1^2 + c1*c2 + c2^2)/(c1+c2);
        area(i, 1) = 0.5*(c1+c2)*sqrt(y^2+z^2);
        S(j) = S(j) + area(i,1);
    end
    cordas(i) = c1;
    reynolds_info(i) = c1*Vinf/nu;
    if nperfil(j) == nperfil(j+1)
        airfoil_info(i) = nperfil(j);
    else
        airfoil_info(i) = (collocation_points(i,2)-vertices(soma_part+1,2))/(0.5*b_particoes(j));
    end
    
    csi(i, 1) = (mac(i)/area(i))*length(i,1);
    csi(i, 2) = (mac(i)/area(i))*length(i,2);
    csi(i, 3) = (mac(i)/area(i))*length(i,3);
end
MAC=0;
for i=1:particoes
    lambda = c(i+1)/c(i);
    MAC = MAC + S(i)*(2/3)*c(i)*(1+lambda+lambda^2)/(1+lambda);
end
MACw = MAC/sum(S);
AR = (b)^2/(2*sum(area));
end
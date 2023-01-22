function [vji,vji_sim] = llt_vind(vertices, collocation_points, vinf, mac, h, N, i)
%V1.2
ri1j = zeros(1, 3);
ri2j = zeros(1, 3);
vinf_ge = [vinf(1) vinf(2) -vinf(3)];

vji = zeros(N,3);
vji_sim = zeros(N,3);
%i: ponto de colocação
%j: variação dos vórtices
    for j=1:N
        if j == i
            %Painel original, vórtice original
            ri2j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
            ri2j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2);
            ri2j(1, 3) = collocation_points(i, 3) - vertices(j+1, 3);
            
            ri1j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
            ri1j(1, 2) = collocation_points(i, 2) - vertices(j, 2);
            ri1j(1, 3) = collocation_points(i, 3) - vertices(j, 3);
            
            aux1_cross = cross(vinf, ri1j);
            aux2_cross = cross(vinf, ri2j);
            aux1_dot = dot(vinf, ri1j);
            aux2_dot = dot(vinf, ri2j);
            
            vji(j,:) = mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)));
            
            if h~=0
                %Adição do efeito solo
                ri1j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
                ri1j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2);
                ri1j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j+1, 3));
                
                ri2j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
                ri2j(1, 2) = collocation_points(i, 2) - vertices(j, 2);
                ri2j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j, 3));
                
                aux1_cross = cross(vinf_ge, ri1j);
                aux2_cross = cross(vinf_ge, ri2j);
                aux1_dot = dot(vinf_ge, ri1j);
                aux2_dot = dot(vinf_ge, ri2j);
                
                vji_groundeffect = mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)));
                vji(j,:) = vji(j,:) + vji_groundeffect;
            end
            
            %Painel original, vórtice oposto
            ri1j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
            ri1j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2)*-1;
            ri1j(1, 3) = collocation_points(i, 3) - vertices(j+1, 3);
            
            ri2j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
            ri2j(1, 2) = collocation_points(i, 2) - vertices(j, 2)*-1;
            ri2j(1, 3) = collocation_points(i, 3) - vertices(j, 3);
            
            aux1_cross = cross(vinf, ri1j);
            aux2_cross = cross(vinf, ri2j);
            aux1_dot = dot(vinf, ri1j);
            aux2_dot = dot(vinf, ri2j);
            aux_cross = cross(ri1j, ri2j);
            aux_dot = dot(ri1j, ri2j);
            
            vji_sim(j,:) =  mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot)));
            
            if h~=0
                %Efeito solo
                ri2j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
                ri2j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2)*-1;
                ri2j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j+1, 3));
                
                ri1j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
                ri1j(1, 2) = collocation_points(i, 2) - vertices(j, 2)*-1;
                ri1j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j, 3));
                
                aux1_cross = cross(vinf_ge, ri1j);
                aux2_cross = cross(vinf_ge, ri2j);
                aux1_dot = dot(vinf_ge, ri1j);
                aux2_dot = dot(vinf_ge, ri2j);
                aux_cross = cross(ri1j, ri2j);
                aux_dot = dot(ri1j, ri2j);
                
                vji_sim_groundeffect =  mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                    + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot)));
                vji_sim(j,:) = vji_sim(j,:) + vji_sim_groundeffect;
            end
        else
            %Painel original, vórtice original
            ri2j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
            ri2j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2);
            ri2j(1, 3) = collocation_points(i, 3) - vertices(j+1, 3);
            
            ri1j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
            ri1j(1, 2) = collocation_points(i, 2) - vertices(j, 2);
            ri1j(1, 3) = collocation_points(i, 3) - vertices(j, 3);
            
            aux1_cross = cross(vinf, ri1j);
            aux2_cross = cross(vinf, ri2j);
            aux1_dot = dot(vinf, ri1j);
            aux2_dot = dot(vinf, ri2j);
            aux_cross = cross(ri1j, ri2j);
            aux_dot = dot(ri1j, ri2j);
            
            vji(j,:) = mac(i)/(4*pi)*((aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot))));
            
            if h~=0
                %Efeito solo
                ri1j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
                ri1j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2);
                ri1j(1, 3) = collocation_points(i, 3) - (2*h+2*vertices(j+1, 3));
                
                ri2j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
                ri2j(1, 2) = collocation_points(i, 2) - vertices(j, 2);
                ri2j(1, 3) = collocation_points(i, 3) - (2*h+2*vertices(j, 3));
                
                aux1_cross = cross(vinf_ge, ri1j);
                aux2_cross = cross(vinf_ge, ri2j);
                aux1_dot = dot(vinf_ge, ri1j);
                aux2_dot = dot(vinf_ge, ri2j);
                aux_cross = cross(ri1j, ri2j);
                aux_dot = dot(ri1j, ri2j);
                
                vji_groundeffect = mac(i)/(4*pi)*((aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                    + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot))));
                
                vji(j,:) = vji(j,:) + vji_groundeffect;
            end
            
            %Painel original, vórtice oposto
            ri1j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
            ri1j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2)*-1;
            ri1j(1, 3) = collocation_points(i, 3) - vertices(j+1, 3);
            
            ri2j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
            ri2j(1, 2) = collocation_points(i, 2) - vertices(j, 2)*-1;
            ri2j(1, 3) = collocation_points(i, 3) - vertices(j, 3);
            
            aux1_cross = cross(vinf, ri1j);
            aux2_cross = cross(vinf, ri2j);
            aux1_dot = dot(vinf, ri1j);
            aux2_dot = dot(vinf, ri2j);
            aux_cross = cross(ri1j, ri2j);
            aux_dot = dot(ri1j, ri2j);
            
            vji_sim(j,:) =  mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot)));
            
            if h~=0
                %Efeito solo
                ri2j(1, 1) = collocation_points(i, 1) - vertices(j+1, 1);
                ri2j(1, 2) = collocation_points(i, 2) - vertices(j+1, 2)*-1;
                ri2j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j+1, 3));
                
                ri1j(1, 1) = collocation_points(i, 1) - vertices(j, 1);
                ri1j(1, 2) = collocation_points(i, 2) - vertices(j, 2)*-1;
                ri1j(1, 3) = collocation_points(i, 3) - (2*h + 2*vertices(j, 3));
                
                aux1_cross = cross(vinf_ge, ri1j);
                aux2_cross = cross(vinf_ge, ri2j);
                aux1_dot = dot(vinf_ge, ri1j);
                aux2_dot = dot(vinf_ge, ri2j);
                aux_cross = cross(ri1j, ri2j);
                aux_dot = dot(ri1j, ri2j);
                
                vji_sim_groundeffect = mac(i)/(4*pi)*(aux2_cross/(norm(ri2j)*(norm(ri2j)-aux2_dot)) - aux1_cross/(norm(ri1j)*(norm(ri1j)-aux1_dot)) ...
                    + (norm(ri1j) + norm(ri2j))*(aux_cross)/( norm(ri1j)*norm(ri2j) * (norm(ri1j)*norm(ri2j) + aux_dot)));
                
                vji_sim(j,:) = vji_sim(j,:) + vji_sim_groundeffect;
            end
        end
    end
end
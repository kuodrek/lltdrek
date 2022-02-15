function [particoes_novo,b_novo,c_novo,diedro_novo,offset_novo,nperfil_novo]=llt_autogap(particoes,b,c,diedro,offset,nperfil,gap)
offset = [0 offset];
aux_part = 0;
aux_index = 0;
for i=1:size(nperfil,2)-1
    if nperfil(i) ~= nperfil(i+1) && 0.7*b(i) > gap
        aux_index(i) = i + aux_part;
        aux_part = aux_part + 1;
    else
        aux_index(i) = -(i + aux_part);
    end
end
aux_b = zeros(particoes+aux_part,1);
aux_c = zeros(particoes+aux_part+1,1);
aux_c(1) = c(1);
aux_diedro = zeros(size(aux_b,1),1);
aux_nperfil = zeros(size(aux_b,1)+1,1);
aux_nperfil(1) = nperfil(1);
aux_offset = zeros(size(aux_b,1)+1,1);
for i=1:size(aux_index,2)
    if aux_index(i)>0
    aux_b(aux_index(i)) = b(i)-gap;
    aux_c(aux_index(i)+1) = ( c(i+1)-c(i) )/b(i)*( b(i) - gap) + c(i);
    aux_c(aux_index(i)+2) = c(i+1);
    aux_diedro(aux_index(i)) = diedro(i);
    aux_nperfil(aux_index(i)) = nperfil(i);
    aux_nperfil(aux_index(i)+1) = nperfil(i);
    aux_offset(aux_index(i)+1) = ( offset(i+1)-offset(i) )/b(i)*( b(i) - gap) + offset(i);
    aux_offset(aux_index(i)+2) = offset(i+1);
    else
        aux_b(abs(aux_index(i))) = b(i);
        aux_c(abs(aux_index(i))+1) = c(i+1);
        aux_diedro(abs(aux_index(i))) = diedro(i);
        aux_nperfil(abs(aux_index(i))) = nperfil(i);
        aux_offset(abs(aux_index(i))+1) = offset(i+1);
    end
end
aux_offset(end) = offset(end);
aux_nperfil(end) = nperfil(end);
for i=1:size(aux_b,1)
    if aux_b(i) == 0
        aux_b(i) = gap;
        aux_diedro(i) = aux_diedro(i-1);
    end
end
aux_c(end) = c(end);

particoes_novo = particoes + aux_part;
b_novo = aux_b';
c_novo = aux_c';
diedro_novo = aux_diedro';
nperfil_novo = aux_nperfil';
offset_novo = aux_offset(2:end)';
end
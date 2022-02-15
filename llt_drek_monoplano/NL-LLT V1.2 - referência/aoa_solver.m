function [data] = aoa_solver(cl_lista,aoa,tipo)
% V1.2

if strcmp(tipo,'cl') == 1
%     if aoa<cl_lista(1,1) || aoa>cl_lista(end,1)
%         fprintf("\nAlfa fora da faixa de valores.\nAlfa: %g", aoa);
%     end
    if aoa <= cl_lista(1,1)
        cl = ( cl_lista(2,2)-cl_lista(1,2) )/( cl_lista(2,1)-cl_lista(1,1) )*(aoa-cl_lista(1,1)) + cl_lista(1,2);
        data = cl;
    elseif aoa >= cl_lista(end,1)
        cl = ( cl_lista(end,2)-cl_lista(end-1,2) )/( cl_lista(end,1)-cl_lista(end-1,1) )*(aoa-cl_lista(end-1,1)) + cl_lista(end-1,2);
        data = cl;
    else
        for i=1:size(cl_lista,1)-1
            if aoa == cl_lista(i,1)
                cl = cl_lista(i,2);
                data = cl;
            elseif aoa > cl_lista(i,1) && aoa < cl_lista(i+1,1)
                aoa1 = cl_lista(i,1);
                cl1 = cl_lista(i,2);
                aoa2 = cl_lista(i+1,1);
                cl2 = cl_lista(i+1,2);
                cl = cl1 + (cl2-cl1)/(aoa2-aoa1)*(aoa-aoa1);
                data = cl;
            end
        end
        
    end
elseif strcmp(tipo,'cl_alfa') == 1
    if aoa <= cl_lista(1,1)
        cl = (cl_lista(2,2)-cl_lista(1,2))/(cl_lista(2,1)-cl_lista(1,1));
        data = cl;
    elseif aoa >= cl_lista(end,1)
        cl = (cl_lista(end,2)-cl_lista(end-1,2))/(cl_lista(end,1)-cl_lista(end-1,1));
        data = cl;
    else
        for i=1:size(cl_lista,1)-1
            if aoa == cl_lista(i,1)
                cl = ( cl_lista(i+1,2) - cl_lista(i,2) )/( cl_lista(i+1,1) - cl_lista(i,1) );
                data = cl;
            elseif aoa > cl_lista(i,1) && aoa < cl_lista(i+1,1)
                aoa1 = cl_lista(i,1);
                cl1 = cl_lista(i,2);
                aoa2 = cl_lista(i+1,1);
                cl2 = cl_lista(i+1,2);
                cl = (cl2-cl1)/(aoa2-aoa1);
                data = cl;
            end
        end
        
    end
else
    error('Erro em cl_alfa_solver: input tipo usado de maneira incorreta (esperava-se cl ou cl_alfa)')
end

end
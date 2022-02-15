function data = re_solver(re,v_coletor,aoa,tipo)
% V1.2

% Encontrar o Reynolds relevante
if re <= v_coletor{1,1}
    if strcmp(tipo,'cm') == 1
        cm = cell2mat(v_coletor(3,1));
        data = cm;
    else
        cl_lista = cell2mat(v_coletor(2,1));
        cl = aoa_solver(cl_lista,aoa,tipo);
        data = cl;
    end
elseif re >= v_coletor{1,end}
    if strcmp(tipo,'cm') == 1
        cm = cell2mat(v_coletor(3,end));
        data = cm;
    else
        cl_lista = cell2mat(v_coletor(2,end));
        cl = aoa_solver(cl_lista,aoa,tipo);
        data = cl;
    end
else
    for j=1:size(v_coletor,2)-1
        if re == v_coletor{1,j}
            if strcmp(tipo,'cm') == 1
                cm = cell2mat(v_coletor(3,j));
                data = cm;
            else
                cl_lista = cell2mat(v_coletor(2,j));
                cl = aoa_solver(cl_lista,aoa,tipo);
                data = cl;
            end
        elseif re > v_coletor{1,j} && re < v_coletor{1,j+1}
            if strcmp(tipo,'cm') == 1
                re1 = v_coletor{1,j};
                cm1 = cell2mat(v_coletor(3,j));
                re2 = v_coletor{1,j+1};
                cm2 = cell2mat(v_coletor(3,j+1));
                cm = (cm2-cm1)/(re2-re1)*(re-re1) + cm1;
                data = cm;
            else
                re1 = v_coletor{1,j};
                cl_lista1 = cell2mat(v_coletor(2,j));
                cl1 = aoa_solver(cl_lista1,aoa,tipo);
                re2 = v_coletor{1,j+1};
                cl_lista2 = cell2mat(v_coletor(2,j+1));
                cl2 = aoa_solver(cl_lista2,aoa,tipo);
                cl = (cl2-cl1)/(re2-re1)*(re-re1) + cl1;
                data = cl;
            end
        end
    end
end

end
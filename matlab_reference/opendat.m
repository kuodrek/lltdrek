function [airfoildat] = opendat(airfoil)

datfile = fopen(strcat(airfoil,'.dat'));

cont_dat = 0;
cont = 0;
airfoildat = cell(3,2);

airfoildat{1,2} = 'Perfil';
airfoildat{2,2} = 'Extradorso';
airfoildat{3,2} = 'Intradorso';
while cont_dat == 0
    if cont == 0
        nome_perfil = fscanf(datfile,'%s',1);
        cont = cont + 1;
        airfoildat{1,1} = nome_perfil;
    else
        dat_aux = fscanf(datfile,'%f %f',[1 2]);
        if feof(datfile) == 1
            cont_dat = 1;
        else
            dat(cont,:) = dat_aux;
            cont = cont + 1;
        end
    end
    
end
fclose(datfile);

extradorso = [1 0];
aux = dat(1,1);
for i=1:size(dat,1)
    if dat(i,1) < aux || i == 1
        extradorso(i,:) = dat(i,:);
        aux = dat(i,1);
    else
        cont = i;
        break
    end
end
intradorso(1,:) = [0 0];
j = 2;
for i=cont:size(dat,1)
    intradorso(j,:) = dat(i,:);
    j = j + 1;
end
aux = zeros(size(extradorso,1),2);
final = size(extradorso,1);
for i=1:size(extradorso,1)
    aux(i,:) = extradorso(final+1-i,:);
end
extradorso = aux;

airfoildat{2,1} = extradorso;
airfoildat{3,1} = intradorso;
end
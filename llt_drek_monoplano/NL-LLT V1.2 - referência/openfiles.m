function [airfoildata] = openfiles(airfoil)
% openfiles - V1.2

% INPUT: nome do perfil (string) - lembrar de colocar aspas simples ('optfoilB2')
% OUTPUT: Célula cuja primeira linha são os números de Reynolds e a segunda linha são os dados (Cl x alpha)
% Lembrando que o comando cell2mat converte a célula para matriz (útil para manipular dados)

% Conversão de tipo de variável - char (int) (tabela ASCII) (obs: os valores variam com o caps lock)
% R (82) E (69)
% f (102) i (105) m (109)

% strcat é uma função que concatena horizontalmente strings ('perfil' + '.txt')
% fopen vai abrir o .txt do perfil solicitado
airfoilfile = fopen(strcat(airfoil,'.txt'));

cont_reynolds = 0;  % Variável auxiliar que serve pra terminar o laço de repetição mestre
cont = 0;           % Contador que serve para separar as informações entre Reynolds
airfoildata = {};   % Célula de output (pré alocar ela no futuro!)

while cont_reynolds == 0
    % fscanf lê a linha do arquivo aberto; Aqui, a função está tentando
    % achar 'RE', onde a variável %f será o numero de reynolds
    string_reynolds = fscanf(airfoilfile,'%s %f',[1 2]);
    
    if isempty(string_reynolds) == 0
        if string_reynolds(1) == 82
            % Ângulo final do Cl x alpha. Assim, é possível saber quando os dados chegam no final
            angulo_aux = fscanf(airfoilfile,'%s %f',[1 2]);
            angulo_maximo = angulo_aux(end);
            Cm0_aux = fscanf(airfoilfile,'%s %f',[1 2]);
            Cm0 = Cm0_aux(end);
            cont = cont + 1;
            % Número de Reynolds dos dados
            airfoildata{1,cont} = string_reynolds(end);
            cont_data = 0;
            data_vetor = [];
            while cont_data == 0
                data_aux = fscanf(airfoilfile,'%f %f',[1 2]);
                % Alocação dos dados Cl x alpha em um vetor auxiliar (pré alocá-lo no futuro!)
                data_vetor(end+1,:) = [data_aux(1),data_aux(2)]; 
                if data_aux(1) == angulo_maximo
                    airfoildata{2,cont} = data_vetor;
                    airfoildata{3,cont} = Cm0;
                    cont_data = 1;
                end
            end
        end
    else
        cont_reynolds = 1;
    end
end
% Fechar a variável que representa o arquivo do perfil
fclose(airfoilfile);
end
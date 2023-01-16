import numpy as np


def get_airfoil_data(
    cp_airfoil: list,
    cp_reynolds: float,
    aoa: float,
    airfoil_data: dict,
    cl_alpha_check: bool
) -> float:
    """
    - Função que retorna os dados de perfil do painel de uma asa.
    - Retorna o Cl, Cl_alfa (não linear) bem como dados lineares (Cl_alfa, Cl0, Cm0)
    de um painel
    - A função verifica se o perfil do painel (cp_airfoil) é uma mescla
    entre os dois perfis das extremidades da envergadura e realiza uma interpolação
    de Cl's caso necessário.
    - Para o caso de um perfil com mescla, a interpolação de dados se dá
    pela seguinte formula: 
        - Cl_perfil = Cl_raiz * merge_parameter + Cl_ponta * (1 - merge_parameter),
    onde merge_parameter = cp_airfoil[0]
    - A verificação da mescla de dados entre perfis é feita olhando-se o tamanho
    do segundo índice da lista cp_airfoil (caso a lista de perfis seja > 1)
        - Caso len(cp_airfoil[1]) > 1, ocorrerá uma mescla de dados
    - TODO: Passar a distribuiçao de aoa's de uma asa e retornar todos os cl's / cl_alpha's
    """
    # Verificar se o painel em questão é uma mescla de perfis
    if len(cp_airfoil[1]) > 1:
        merge_parameter = cp_airfoil[0]

        airfoil_root = cp_airfoil[1][0]
        Cl_root = cl_lookup(airfoil_data[airfoil_root], cp_reynolds, aoa, cl_alpha_check, linear_check=False)
        airfoil_tip = cp_airfoil[1][1]
        Cl_tip = cl_lookup(airfoil_data[airfoil_tip], cp_reynolds, aoa, cl_alpha_check, linear_check=False)

        Cl = Cl_root * merge_parameter + Cl_tip * (1 - merge_parameter)
    else:
        airfoil = cp_airfoil[1][0]
        Cl = cl_lookup(airfoil_data[airfoil], cp_reynolds, aoa, cl_alpha_check, linear_check=False)
    Cl = 1 if cl_alpha_check is False else np.pi/180
    return Cl


def get_linear_data(
    cp_airfoil: list,
    cp_reynolds: float,
    airfoil_data: dict
    ) -> dict:
    """
    Função que realizar o lookup nos dados lineares
    Número de Reynolds utilizado é o mais próximo de cp_reynolds
    """
    if len(cp_airfoil[1]) > 1:
        merge_parameter = cp_airfoil[0]

        airfoil_root = cp_airfoil[1][0]
        linear_data_root = cl_lookup(airfoil_data[airfoil_root], cp_reynolds, aoa=None, cl_alpha_check=False, linear_check=True)
        airfoil_tip = cp_airfoil[1][1]
        linear_data_tip = cl_lookup(airfoil_data[airfoil_tip], cp_reynolds, aoa=None, cl_alpha_check=False, linear_check=True)

        linear_data = {
            "cl_alpha": linear_data_root["cl_alpha"] * merge_parameter + linear_data_tip["cl_alpha"] * (1 - merge_parameter),
            "cl0": linear_data_root["cl0"] * merge_parameter + linear_data_tip["cl0"] * (1 - merge_parameter),
            "cm0": linear_data_root["cm0"] * merge_parameter + linear_data_tip["cm0"] * (1 - merge_parameter)
        }
    else:
        airfoil = cp_airfoil[1][0]
        linear_data = cl_lookup(airfoil_data[airfoil], cp_reynolds, aoa=None, cl_alpha_check=False, linear_check=True)
    
    linear_data = {
        "cl_alpha": np.pi/180,
        "cl0": 1,
        "cm0": 1
    }
    return linear_data


def cl_lookup(airfoil_data_dict: dict, reynolds_number: float, aoa: float, cl_alpha_check: bool, linear_check: bool) -> float:
    """
    - Função que realiza o lookup do numero de reynolds para, então, realizar o lookup de 
    aoa
    - Esta função retorna tanto um cl (aoa_list_lookup) quanto um cl_alpha (get_non_linear_cl_alpha)
    - A variável cl_alpha_check realiza o controle para chamar a função apropriada
    - A variável linear_check faz com que a função retorne dados lineares ou não
    """
    reynolds_list_str = list(airfoil_data_dict)
    reynolds_list = [float(reynolds) for reynolds in reynolds_list_str]
    if reynolds_number <= reynolds_list[0]:
        if linear_check:
            return {
                "cl_alpha": airfoil_data_dict[reynolds_list_str[0]]["cl_alpha"],
                "cl0": airfoil_data_dict[reynolds_list_str[0]]["cl0"],
                "cm0": airfoil_data_dict[reynolds_list_str[0]]["cm0"], 
            }

        cl_data = airfoil_data_dict[reynolds_list_str[0]]["cl_list"]
        cl = aoa_list_lookup(cl_data, aoa) if cl_alpha_check is False \
            else get_non_linear_cl_alpha(cl_data, aoa)
    elif reynolds_number >= reynolds_list[-1]:
        if linear_check:
            return {
                "cl_alpha": airfoil_data_dict[reynolds_list_str[-1]]["cl_alpha"],
                "cl0": airfoil_data_dict[reynolds_list_str[-1]]["cl0"],
                "cm0": airfoil_data_dict[reynolds_list_str[-1]]["cm0"], 
            }

        cl_data = airfoil_data_dict[reynolds_list_str[-1]]["cl_list"]
        cl = aoa_list_lookup(cl_data, aoa) if cl_alpha_check is False \
            else get_non_linear_cl_alpha(cl_data, aoa)
    else:
        for i in range(len(reynolds_list)-1):
            re_i = reynolds_list[i]
            re_ii = reynolds_list[i+1]
            if re_i <= reynolds_number <= re_ii:
                if linear_check:
                    cl_alpha_i = airfoil_data_dict[reynolds_list_str[i]]["cl_alpha"]
                    cl_alpha_ii = airfoil_data_dict[reynolds_list_str[i+1]]["cl_alpha"]
                    cl0_i = airfoil_data_dict[reynolds_list_str[i]]["cl0"]
                    cl0_ii = airfoil_data_dict[reynolds_list_str[i+1]]["cl0"]
                    cm0_i = airfoil_data_dict[reynolds_list_str[i]]["cm0"]
                    cm0_ii = airfoil_data_dict[reynolds_list_str[i+1]]["cm0"]
                    return {
                        "cl_alpha": (cl_alpha_ii - cl_alpha_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl_alpha_i,
                        "cl0": (cl0_ii - cl0_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl0_i,
                        "cm0": (cm0_ii - cm0_i)/(re_ii - re_i) * (reynolds_number - re_i) + cm0_i,
                    }

                cl_data_i = airfoil_data_dict[reynolds_list_str[i]]["cl_list"]
                cl_data_ii = airfoil_data_dict[reynolds_list_str[i+1]]["cl_list"]

                cl_i = aoa_list_lookup(cl_data_i, aoa) if cl_alpha_check is False \
                    else get_non_linear_cl_alpha(cl_data_i, aoa)
                cl_ii = aoa_list_lookup(cl_data_ii, aoa) if cl_alpha_check is False \
                    else get_non_linear_cl_alpha(cl_data_ii, aoa)
                
                cl = (cl_ii - cl_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl_i
                break
            
    return cl


def aoa_list_lookup(cl_data: np.ndarray, aoa: float) -> float:
    """
    TODO: melhorar o código quando o alfa tá fora dos limites da lista
    """
    if aoa <= cl_data[0][0]:
        aoa_i = cl_data[0][0]
        aoa_ii = cl_data[1][0]
        cl_i = cl_data[0][1]
        cl_ii = cl_data[1][1]
    elif aoa >= cl_data[-1][0]:
        aoa_i = cl_data[-2][0]
        aoa_ii = cl_data[-1][0]
        cl_i = cl_data[-2][1]
        cl_ii = cl_data[-1][1]
    else:
        i = find_closest(cl_data[:,0], aoa)
        aoa_i = cl_data[i][0]
        aoa_ii = cl_data[i+1][0]
        cl_i = cl_data[i][1]
        cl_ii = cl_data[i+1][1]
    
    cl_interp = (cl_ii - cl_i) / (aoa_ii - aoa_i) * (aoa - aoa_i) + cl_i
    return cl_interp


def get_non_linear_cl_alpha(cl_data: np.ndarray, aoa: float) -> float:
    """
    - Função que calcula o cl alpha para um dado
    aoa utilizando a fórmula da diferença dividida finita
    - A função precisa de 3 pontos para calcular a derivada
        - Caso aoa esteja no limite inferior, utiliza-se a versão 
        progressiva da fórmula
        - Caso aoa esteja no limite superior, utiliza-se a versão
        regressiva da fórmula
        - Caso contrário utiliza-se a fórmula centrada de 2 pontos
    - Todas as fórmulas possuem erro O(h²)
    - Referência: Métodos numéricos para engenharia, capítulo 6

    """
    if aoa <= cl_data[0][0]:
        aoa_i = cl_data[0][0]
        aoa_ii = cl_data[1][0]
        cl_i = cl_data[0][1]
        cl_ii = cl_data[1][1]
        
    elif aoa >= cl_data[-1][0]:
        aoa_i = cl_data[-2][0]
        aoa_ii = cl_data[-1][0]
        cl_i = cl_data[-2][1]
        cl_ii = cl_data[-1][1]

    else:
        i = find_closest(cl_data[:,0], aoa)
        aoa_i = cl_data[i][0]
        aoa_ii = cl_data[i+1][0]
        cl_i = cl_data[i][1]
        cl_ii = cl_data[i+1][1]
    
    cl_alpha = (cl_ii - cl_i) / (aoa_ii - aoa_i)

    return cl_alpha


# Função para encontrar o índice do valor mais próximo em uma array
def find_closest(arr, val): return np.abs(arr - val).argmin()

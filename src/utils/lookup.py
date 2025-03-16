import numpy as np


def get_airfoil_data(
    cp_airfoil: list,
    cp_reynolds: float,
    aoa: float,
    airfoil_data: dict,
    cl_alpha_check: bool,
    show_logs: bool = True
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
        - Cl_perfil = Cl_raiz * (1 - merge_parameter) + Cl_ponta * (1 - merge_parameter),
    onde merge_parameter = cp_airfoil[0]
    - A verificação da mescla de dados entre perfis é feita olhando-se o tamanho
    do segundo índice da lista cp_airfoil (caso a lista de perfis seja > 1)
        - Caso len(cp_airfoil[1]) > 1, ocorrerá uma mescla de dados
    - TODO: Passar a distribuiçao de aoa's de uma asa e retornar todos os cl's / cl_alpha's
    """
    # Verificar se o painel em questão é uma mescla de perfis
    linear_check = False
    if len(cp_airfoil[1]) > 1:
        merge_parameter = cp_airfoil[0]

        airfoil_root = cp_airfoil[1][0]
        Cl_root = cl_lookup(airfoil_data[airfoil_root], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)
        airfoil_tip = cp_airfoil[1][1]
        Cl_tip = cl_lookup(airfoil_data[airfoil_tip], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)

        Cl = Cl_root * (1 - merge_parameter) + Cl_tip * merge_parameter
    else:
        airfoil = cp_airfoil[1][0]
        Cl = cl_lookup(airfoil_data[airfoil], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)
    return Cl


def get_linear_data_and_clmax(
    cp_airfoil: list,
    cp_reynolds: float,
    airfoil_data: dict,
    show_logs: bool = True
    ) -> dict:
    """
    Função que realizar o lookup nos dados lineares e clmax
    Número de Reynolds utilizado é o mais próximo de cp_reynolds
    """
    aoa = None
    cl_alpha_check = False
    linear_check = True
    if len(cp_airfoil[1]) > 1:
        merge_parameter = cp_airfoil[0]

        airfoil_root = cp_airfoil[1][0]
        lookup_data_root = cl_lookup(airfoil_data[airfoil_root], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)
        airfoil_tip = cp_airfoil[1][1]
        lookup_data_tip = cl_lookup(airfoil_data[airfoil_tip], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)

        lookup_data = {
            "cl_alpha": lookup_data_root["cl_alpha"] * (1 - merge_parameter) + lookup_data_tip["cl_alpha"] * merge_parameter,
            "cl0": lookup_data_root["cl0"] * (1 - merge_parameter) + lookup_data_tip["cl0"] * merge_parameter,
            "cm0": lookup_data_root["cm0"] * (1 - merge_parameter) + lookup_data_tip["cm0"] * merge_parameter,
            "clmax": lookup_data_root["clmax"] * (1 - merge_parameter) + lookup_data_tip["clmax"] * merge_parameter,
        }
    else:
        airfoil = cp_airfoil[1][0]
        lookup_data = cl_lookup(airfoil_data[airfoil], cp_reynolds, aoa, cl_alpha_check, linear_check, show_logs)
    
    return lookup_data


def cl_lookup(airfoil_data_dict: dict, reynolds_number: float, aoa: float, cl_alpha_check: bool, linear_check: bool, show_logs: bool = True) -> float:
    """
    - Função que realiza o lookup do numero de reynolds para, então, realizar o lookup de 
    aoa
    - Esta função retorna tanto um cl (aoa_list_lookup) quanto um cl_alpha (get_non_linear_cl_alpha)
    - A variável cl_alpha_check realiza o controle para chamar a função apropriada
    - A variável linear_check faz com que a função retorne dados lineares ou não
    """
    reynolds_list_str = list(airfoil_data_dict)
    reynolds_list = [float(reynolds) for reynolds in reynolds_list_str]
    if show_logs is True and (aoa < reynolds_list[0] or aoa > reynolds_list[-1]):
        print(f"Warning: Reynolds {aoa} out of bounds")
    if reynolds_number <= reynolds_list[0]:
        if linear_check:
            return {
                "cl_alpha": airfoil_data_dict[reynolds_list_str[0]]["cl_alpha"],
                "cl0": airfoil_data_dict[reynolds_list_str[0]]["cl0"],
                "cm0": airfoil_data_dict[reynolds_list_str[0]]["cm0"], 
                "clmax": airfoil_data_dict[reynolds_list_str[0]]["clmax"], 
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
                "clmax": airfoil_data_dict[reynolds_list_str[-1]]["clmax"],
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
                    clmax_i = airfoil_data_dict[reynolds_list_str[i]]["clmax"]
                    clmax_ii = airfoil_data_dict[reynolds_list_str[i+1]]["clmax"]
                    return {
                        "cl_alpha": (cl_alpha_ii - cl_alpha_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl_alpha_i,
                        "cl0": (cl0_ii - cl0_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl0_i,
                        "cm0": (cm0_ii - cm0_i)/(re_ii - re_i) * (reynolds_number - re_i) + cm0_i,
                        "clmax": (clmax_ii - clmax_i)/(re_ii - re_i) * (reynolds_number - re_i) + clmax_i,
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


def aoa_list_lookup(cl_data: np.ndarray, aoa: float, show_logs: bool = True) -> float:
    """
    TODO: melhorar o código quando o alfa tá fora dos limites da lista
    """
    if show_logs is True and (aoa < cl_data[0][0] or aoa > cl_data[-1][0]):
        print(f"Warning: Alpha {aoa} out of bounds")
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
        for i, _ in enumerate(cl_data[0:-1,0]):
            if cl_data[i][0] < aoa < cl_data[i+1][0]:
                aoa_i = cl_data[i][0]
                aoa_ii = cl_data[i+1][0]
                cl_i = cl_data[i][1]
                cl_ii = cl_data[i+1][1]
    
    cl_interp = (cl_ii - cl_i) / (aoa_ii - aoa_i) * (aoa - aoa_i) + cl_i
    return cl_interp


def get_non_linear_cl_alpha(cl_data: np.ndarray, aoa: float, show_logs: bool = True) -> float:
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
    if show_logs is True and (aoa < cl_data[0][0] or aoa > cl_data[-1][0]):
        print(f"Warning: Alpha {aoa} out of bounds")
    
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
        for i, _ in enumerate(cl_data[0:-1,0]):
            if cl_data[i][0] < aoa < cl_data[i+1][0]:
                aoa_i = cl_data[i][0]
                aoa_ii = cl_data[i+1][0]
                cl_i = cl_data[i][1]
                cl_ii = cl_data[i+1][1]
    
    cl_alpha = (cl_ii - cl_i) / (aoa_ii - aoa_i)
    return cl_alpha

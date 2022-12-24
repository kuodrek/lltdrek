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
        Cl_root = cl_lookup(airfoil_data[airfoil_root], cp_reynolds, aoa, cl_alpha_check)
        airfoil_tip = cp_airfoil[1][1]
        Cl_tip = cl_lookup(airfoil_data[airfoil_tip], cp_reynolds, aoa, cl_alpha_check)

        Cl = Cl_root * merge_parameter + Cl_tip * (1 - merge_parameter)
    else:
        airfoil = cp_airfoil[1][0]
        Cl = cl_lookup(airfoil_data[airfoil], cp_reynolds, aoa, cl_alpha_check)
    return Cl


def cl_lookup(airfoil_data_dict: dict, reynolds_number: float, aoa: float, cl_alpha_check: bool) -> float:
    """
    - Função que realiza o lookup do numero de reynolds para, então, realizar o lookup de 
    aoa
    - Esta função retorna tanto um cl (aoa_list_lookup) quanto um cl_alpha (get_non_linear_cl_alpha)
    - A variável cl_alpha_check realiza o controle para chamar a função apropriada
    """
    reynolds_list_str = list(airfoil_data_dict)
    reynolds_list = [float(reynolds) for reynolds in reynolds_list_str]
    if reynolds_number <= reynolds_list[0]:
        cl_data = airfoil_data_dict[reynolds_list_str[0]]["cl_list"]
        cl = aoa_list_lookup(cl_data, aoa) if cl_alpha_check is False \
            else get_non_linear_cl_alpha(cl_data, aoa)
    elif reynolds_number >= reynolds_list[-1]:
        cl_data = airfoil_data_dict[reynolds_list_str[-1]]["cl_list"]
        cl = aoa_list_lookup(cl_data, aoa) if cl_alpha_check is False \
            else get_non_linear_cl_alpha(cl_data, aoa)
    else:
        for i in range(len(reynolds_list)-1):
            re_i = reynolds_list[i]
            re_ii = reynolds_list[i+1]
            if re_i <= reynolds_number <= re_ii:
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
        for i in range(len(cl_data)-1):
            aoa_i = cl_data[i][0]
            aoa_ii = cl_data[i+1][0]
            if aoa_i <= aoa <= aoa_ii:
                cl_i = cl_data[i][1]
                cl_ii = cl_data[i+1][1]
                break
    
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


def get_linear_data(
    cp_airfoil: list,
    cp_reynolds: float,
    airfoil_data: dict
    ) -> dict:
    """
    Função que realizar o lookup nos dados lineares
    Número de Reynolds utilizado é o mais próximo de cp_reynolds
    """
    airfoil_root = cp_airfoil[1][0]
    # dados = airfoil_data["airfoil_name"]["reynolds_number"]["cl_alpha"]
    reynolds_list_str = list(airfoil_data[airfoil_root])
    reynolds_list = [float(reynolds) for reynolds in reynolds_list_str]

    reynolds_dict_translation = {}
    for i, reynolds in enumerate(reynolds_list):
        reynolds_dict_translation[reynolds] = reynolds_list_str[i]
    reynolds_number = reynolds_list[min(range(len(reynolds_list)), key = lambda i: abs(reynolds_list[i]-cp_reynolds))]
    reynolds_key = reynolds_dict_translation[reynolds_number]

    if len(cp_airfoil[1]) > 1:
        merge_parameter = cp_airfoil[0]
        airfoil_tip = cp_airfoil[1][1]

        cl_alpha_root = airfoil_data[airfoil_root][reynolds_key]["cl_alpha"]
        cl0_root = airfoil_data[airfoil_root][reynolds_key]["cl0"]
        cm0_root = airfoil_data[airfoil_root][reynolds_key]["cm0"]

        cl_alpha_tip = airfoil_data[airfoil_tip][reynolds_key]["cl_alpha"]
        cl0_tip = airfoil_data[airfoil_tip][reynolds_key]["cl0"]
        cm0_tip = airfoil_data[airfoil_tip][reynolds_key]["cm0"]
        

        cl_alpha = cl_alpha_root * merge_parameter + cl_alpha_tip * (1 - merge_parameter)
        cl0 = cl0_root * merge_parameter + cl0_tip * (1 - merge_parameter)
        cm0 = cm0_root * merge_parameter + cm0_tip * (1 - merge_parameter)
    else:
        airfoil = airfoil_root
        cl_alpha = airfoil_data[airfoil][reynolds_key]["cl_alpha"]
        cl0 = airfoil_data[airfoil][reynolds_key]["cl0"]
        cm0 = airfoil_data[airfoil][reynolds_key]["cm0"]

    linear_data = {
        "cl_alpha": cl_alpha,
        "cl0": cl0,
        "cm0": cm0,
    }

    return linear_data
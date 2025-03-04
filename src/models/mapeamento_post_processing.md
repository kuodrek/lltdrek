1. verificar se resposta convergiu
    1.1 se nao, retorna tudo nan
2. pegar v_inf, angulo de ataque
3. para cada asa na wing pool
    3.1 Calcular distancia de cada painel em relacao a referencia
    3.2 Para cada cp
        3.2.1 Pegar dado cm0
        3.2.2 Para cada asa na wing pool
            3.2.2.1 Pegar distribuicao de velocidades  
            3.2.2.2 Para cada cp da asa j
                3.2.2.2.1 Calcular produto v_ij * G_j
        3.2.2 Calcular velocidade + v_inf * G_i
        3.2.3 Calcular CF
        3.2.4 Calcular CM
    3.3 Criar matriz de rotacao em relacao a alfa
    3.4 Rebater CF
    Retorna CF, CM, Cl_distr


-> Criar uma classe de retorno que retorna tudo para um alfa
- coeficientes globais do sistema
- coeficientes de cada superficie
- Distribuicoes de cada superficie
- Precisa ser extensivel para retornar mais informacoes no futuro
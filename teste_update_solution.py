complete_wing_pool = [3, 4, 10]
total_value = sum(complete_wing_pool)
G_solution = [i+1 for i in range(total_value)]

def update_solution(G_solution):
    G_list_output = []
    global_counter = 0
    for wing in complete_wing_pool:
        counter = 0
        G_list = []
        while counter < wing:
            G_list.append(G_solution[global_counter])
            counter += 1
            global_counter += 1
        G_list_output.append(G_list)
    return G_list_output

G_list_output = update_solution(G_solution)
print(G_list_output)
airfoils_data_dict = {
    "naca0012": {
        "1.2e5": {
            "cl_list": [[1, 1.1], [2, 1.3], [3, 1.35], [4, 1.4], [5, 1.46]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
        "1e6": {
            "cl_list": [[1, 1.14], [2, 1.34], [3, 1.33], [4, 1.5], [5, 1.6]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
        "1e7": {
            "cl_list": [[1, 1.2], [2, 1.5], [3, 1.7], [4, 1.8], [5, 1.9]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
    },
    "s1223": {
        "1.2e5": {
            "cl_list": [[1, 1.1], [2, 1.3], [3, 1.35], [4, 1.4], [5, 1.46]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
        "1e6": {
            "cl_list": [[1, 1.14], [2, 1.34], [3, 1.33], [4, 1.5], [5, 1.6]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
        "1e7": {
            "cl_list": [[1, 1.2], [2, 1.5], [3, 1.7], [4, 1.8], [5, 1.9]],
            "cl_alfa": 1,
            "cl0": 1,
            "cm0": 1,
        },
    },
}


def get_cl(airfoil_data_dict, reynolds_number, aoa):
    reynolds_list = list(airfoil_data_dict)
    for reynolds_i, cl_data in airfoil_data_dict.items():
        reynolds_i = int(reynolds_i)


def aoa_list_lookup(cl_list, aoa):
    if aoa <= cl_list[0][0]:
        # cl = ( cl_lista(2,2)-cl_lista(1,2) )/( cl_lista(2,1)-cl_lista(1,1) )*(aoa-cl_lista(1,1)) + cl_lista(1,2);
        aoa_i = cl_list[0][0]
        aoa_ii = cl_list[1][0]
        cl_i = cl_list[0][1]
        cl_ii = cl_list[1][1]
    elif aoa >= cl_list[-1][0]:
        aoa_i = cl_list[-2][0]
        aoa_ii = cl_list[-1][0]
        cl_i = cl_list[-2][1]
        cl_ii = cl_list[-1][1]
    else:
        for i in range(len(cl_list)-1):
            aoa_i = cl_list[i][0]
            aoa_ii = cl_list[i+1][0]
            if aoa_i <= aoa <= aoa_ii:
                cl_i = cl_list[i][1]
                cl_ii = cl_list[i+1][1]
                break
    
    cl_interp = (cl_ii - cl_i) / (aoa_ii - aoa_i) * (aoa - aoa_i) + aoa_i
    return cl_interp

# cl = get_cl(airfoils_data_dict)
for index, items in airfoils_data_dict["s1223"].items():
    print(f"item: {items}")
    print(f"index: {index}")

print(list(airfoils_data_dict["s1223"]))

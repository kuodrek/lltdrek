def get_cl(airfoil_data_dict, reynolds_number, aoa):
    reynolds_list_str = list(airfoil_data_dict)
    reynolds_list = [float(reynolds) for reynolds in reynolds_list_str]
    if reynolds_number <= reynolds_list[0]:
        cl_data = airfoil_data_dict[reynolds_list_str[0]]
        cl = aoa_list_lookup(cl_data, aoa)
    elif reynolds_number >= reynolds_list[-1]:
        cl_data = airfoil_data_dict[reynolds_list_str[-1]]
        cl = aoa_list_lookup(cl_data, aoa)
    else:
        for i in range(len(reynolds_list)-1):
            re_i = reynolds_list[i]
            re_ii = reynolds_list[i+1]
            if re_i <= reynolds_number <= re_ii:
                cl_data_i = airfoil_data_dict[reynolds_list_str[i]]
                cl_data_ii = airfoil_data_dict[reynolds_list_str[i+1]]

                cl_i = aoa_list_lookup(cl_data_i, aoa)
                cl_ii = aoa_list_lookup(cl_data_ii, aoa)
                
                cl = (cl_ii - cl_i)/(re_ii - re_i) * (reynolds_number - re_i) + cl_i
                break
    return cl


def aoa_list_lookup(cl_data, aoa):
    cl_list = cl_data["cl_list"]
    if aoa <= cl_list[0][0]:
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
    
    cl_interp = (cl_ii - cl_i) / (aoa_ii - aoa_i) * (aoa - aoa_i) + cl_i
    return cl_interp

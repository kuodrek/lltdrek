import os
import numpy as np
from typing import Union


AOA_POLYFIT_MIN = 0
AOA_POLYFIT_MAX = 8


def load_folder(folder_name: str, aoa_polyfit_min: Union[float, int] = None, aoa_polyfit_max: Union[float, int] = None):
    if aoa_polyfit_min is None or aoa_polyfit_max is None:
        aoa_polyfit_min = AOA_POLYFIT_MIN
        aoa_polyfit_max = AOA_POLYFIT_MAX
    # Get full directory to lookup files
    current_dir = os.path.abspath(os.getcwd())
    target_dir = os.path.join(current_dir, folder_name)
    files_list = []
    for (dir_path, _, filenames) in os.walk(target_dir):
        files_full_path = [(dir_path+"/"+file) for file in filenames]
        files_list.append(files_full_path)
    
    # Transform list of lists into a flat list
    files_flat_list = [item for sublist in files_list for item in sublist]
    txt_list = []
    dat_list = []
    for file in files_flat_list:
        file_extension = file.split(".")[-1]
        if file_extension == "txt": txt_list.append(file)
        elif file_extension == "dat": dat_list.append(file)
    
    airfoils_data_dict = {}
    for file in txt_list:
        airfoil_name = file.split("/")[-1].split(".")[0]
        airfoil_data = cl_file_to_dict(file, aoa_polyfit_min, aoa_polyfit_max)
        airfoils_data_dict[airfoil_name] = airfoil_data

    airfoils_dat_dict = {}
    for file in dat_list:
        airfoil_name = file.split("/")[-1].split(".")[0]
        dat_list = dat_file_to_dict(file)
        airfoils_dat_dict[airfoil_name] = dat_list
    
    return airfoils_data_dict, airfoils_dat_dict



def cl_file_to_dict(file_path: str, aoa_polyfit_min, aoa_polyfit_max):
    airfoil_data_dict = {}
    with open(file_path, "r") as f:
        while True:
            line = f.readline()
            if line == "": break
            reynolds = line.rstrip().split(",")[-1]
            cm0_line = f.readline()
            cm0 = float(cm0_line.rstrip().split(",")[-1])
            cl_list_temp = []
            while True:
                cl_line = f.readline()
                if cl_line == "\n" or cl_line == "": break
                else:
                    cl_line = cl_line.rstrip().split(",")
                    cl_list_temp.append([float(value) for value in cl_line])
            cl_linear_coefs = get_linear_coefs(cl_list_temp, aoa_polyfit_min, aoa_polyfit_max)
            cl_list = np.array(cl_list_temp)
            clmax = cl_list[:,1].max()
            airfoil_data_dict[reynolds] = {
                "cl_list": cl_list,
                "cl_alpha": cl_linear_coefs["cl_alpha"],
                "cl0": cl_linear_coefs["cl0"],
                "cm0": cm0,
                "clmax": clmax
                }
    return airfoil_data_dict


def dat_file_to_dict(file_path: str):
    dat_list = []
    dat_file = np.loadtxt(file_path,skiprows=1)
    for row in dat_file:
        dat_list.append([coordinate for coordinate in row])
    return dat_list


def get_linear_coefs(cl_list: list, aoa_polyfit_min, aoa_polyfit_max):
    aoa_interp_list = []
    cl_interp_list = []
    for values in cl_list:
        if aoa_polyfit_min <= values[0] <= aoa_polyfit_max:
            aoa_interp_list.append(values[0])
            cl_interp_list.append(values[1])
    linear_coefs = np.polyfit(aoa_interp_list, cl_interp_list,1)
    cl_linear_coefs_dict = {
        "cl_alpha": linear_coefs[0],
        "cl0": linear_coefs[1]
    }
    return cl_linear_coefs_dict
        

def xfoil_file_to_cl_file(file_path: str):
    # TODO: implementar
    pass

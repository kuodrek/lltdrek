import os


def load_folder(folder_name: str):
    # Get full directory to lookup files
    current_dir = os.path.abspath(os.getcwd())
    target_dir = os.path.join(current_dir, folder_name)
    files_list = []
    for (dir_path, _, filenames) in os.walk(target_dir):
        files_full_path = [(dir_path+file) for file in filenames]
        files_list.append(files_full_path)
    
    # Transform list of lists into a flat list
    files_flat_list = [item for sublist in files_list for item in sublist]

    for file in files_flat_list:
        if file.endswith('.txt') is True:
            airfoil_cl_list = cl_to_list(file)
        else:
            airfoil_dat_list = dat_to_list(file)



def cl_to_list(file: str):
    pass


def dat_to_list(file: str):
    pass


def xfoil_file_to_cl_file(file: str):
    pass
import os


def load_folder(folder_name: str, data_dict: dict = None):
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
            cl_file_to_dict(file)
        else:
            dat_file_to_dict(file)



def cl_file_to_dict(file: str):
    pass


def dat_file_to_dict(file: str):
    pass


def xfoil_file_to_cl_file(file: str):
    pass
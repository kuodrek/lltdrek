from utils.lookup import get_cl
from airfoil_data_test import airfoils_data_dict


s1223_data = airfoils_data_dict["s1223"]
reynolds_number = 5e6
aoa = 2
cl_1 = get_cl(s1223_data, reynolds_number, aoa)
print(cl_1)
a=1

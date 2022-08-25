cl_dict = {
    "naca0012": {
        "1.2e5": [[1, 1.1], [2, 1.3], [3, 1.35], [4, 1.4], [5, 1.46]],
        "1e6": [[1, 1.14], [2, 1.34], [3, 1.33], [4, 1.5], [5, 1.6]],
        "1e7": [[1, 1.2], [2, 1.5], [3, 1.7], [4, 1.8], [5, 1.9]],
    },
    "s1223": {
        "1.2e5": [[1, 1.1], [2, 1.3], [3, 1.35], [4, 1.4], [5, 1.46]],
        "1e6": [[1, 1.14], [2, 1.34], [3, 1.33], [4, 1.5], [5, 1.6]],
        "1e7": [[1, 1.2], [2, 1.5], [3, 1.7], [4, 1.8], [5, 1.9]],
    },
}


for Re_i, Cl_i in cl_dict.items():
    print(float(Re_i))
    print(Cl_i)

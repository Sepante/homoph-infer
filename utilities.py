import numpy as np
def growth_network_param_extractor( degrees, modes, phi ):
    s_f, s_m = degrees[ modes==0 ].sum(), degrees[ modes==1 ].sum()

    s_f = int(s_f)
    s_m = int(s_m)
    phi = phi.astype('int')

    s_min = np.min([s_f, s_m])

    if s_f % 2 == 1:
        fixing_value = (np.random.randint(2) * 2) - 1
        # print(s_f, s_m, phi)
        s_f += fixing_value
    if s_m % 2 == 1:
        fixing_value = (np.random.randint(2) * 2) - 1
        s_m += fixing_value
        # print(s_f, s_m, phi)
    if phi.sum() % 2 == 1:
        fixing_value = (np.random.randint(2) * 2) - 1
        phi[0] += fixing_value
        # print(s_f, s_m, phi)

    return s_f, s_m, phi

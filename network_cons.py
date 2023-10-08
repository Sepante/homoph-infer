import numpy as np

def init( step_num, mode_num, proportion, h, late_comer_mode, t_emerge, symmetric = True, h2 = 0 ):
    degrees = np.zeros( step_num, dtype = 'int' )

    #init state
    degrees[:2] = 1 #a dyad
    modes = np.random.choice(mode_num, size = step_num, p = [proportion,  1 - proportion])

    late_entrance( modes, late_comer_mode, t_emerge )

    if symmetric:
        homophily_mat = symmetric_homophily_construct( h, step_num, mode_num, modes )
        # print('symmetric')
    else:
        homophily_mat = asymmetric_homophily_construct( h, h2, step_num, mode_num, modes )

    # print(homophily_mat)
    return degrees, modes, homophily_mat

def fitness_init( step_num, distribution ):
    if distribution == 'identical':
        return np.ones( step_num )
    elif distribution == 'uniform':
        return np.random.uniform( size = step_num )

def late_entrance( modes, late_comer_mode, t_emerge ): #currently for only two modes

    modes[ :t_emerge ] = 1 - late_comer_mode #currently for only two modes
    #modes[ len(modes) - t_emerge: ] = late_comer_mode
    #currently only replaces the late_comers in the begining, with the other group


def symmetric_homophily_construct( h, step_num, mode_num, modes ):
    homophily_mat = np.zeros( ( mode_num , step_num) )
    for mode in range( mode_num ):
        homophily_mat[mode, modes == mode] = h
        homophily_mat[mode, modes != mode] = 1-h
    return homophily_mat

def asymmetric_homophily_construct( h1, h2, step_num, mode_num, modes ):
    homophily_mat = np.zeros( ( mode_num , step_num) )
    h_array = np.array([ h1, h2 ])
    for mode in range( mode_num ):
        homophily_mat[mode, modes == mode] = h_array[mode]
        homophily_mat[mode, modes != mode] = 1 - h_array[mode]
    return homophily_mat


def net_grow( step_num, mode_num, proportion, h, k_init, late_comer_mode\
 , t_emerge, t_unbias, distribution, history_step_size, eps\
 , inter_group_pair_num_output = False, inter_group_per_node = False, symmetric = True, h2 = None\
 , degree_symmetry = True, k_init_2 = None, P_A = True):
#    print('symmetric=', symmetric)
    

    degrees, modes, homophily_mat = init( step_num, mode_num, proportion, h, late_comer_mode, t_emerge, symmetric, h2)

    fitnesses = fitness_init( step_num, distribution )

    importances = degrees.astype( 'float' )
    attachment_probabilities = degrees.astype( 'float' )



    if history_step_size:
        history_steps = np.arange( history_step_size, step_num, history_step_size )
        degrees_history = np.zeros( shape = [ len( history_steps )  ] + list( degrees.shape ) )

    if not degree_symmetry:
        k_init_seq = k_init, k_init_2

    #value = 1
    #eps_inv = 1 / eps
#     print(importances)
    if inter_group_pair_num_output:
        # inter_group_pair_num = 0
        inter_group_pair_num = np.zeros(2)
    if inter_group_per_node:
        external_links = np.zeros_like( degrees )

    for step in range(2, step_num):
        # print(  (importances == degrees).sum() )


#         print( 'step: ', step )

#         print( degrees[ degrees > 0 ] )


        mode = modes[step]

        if not degree_symmetry: k_init = k_init_seq[mode]


        degrees[step] += k_init #newcomer getting edges
        importances[step] += k_init #* value



        #old nodes getting edges
        #attachment_probabilities = ( degrees[:step] * homophily_mat[mode, :step] )
        #attachment_probabilities /= attachment_probabilities.sum()
        last_attachment_probabilities = np.copy( attachment_probabilities )


        # attachment_probabilities = ( importances[:step] * homophily_mat[mode, :step] * fitnesses[:step] )
        # attachment_probabilities = ( importances[:step] * homophily_mat[mode, :step] * fitnesses[:step] * degrees[:step] )
        # print( np.sum( degrees == importances ) )
        #super temp:
        # attachment_probabilities = ( homophily_mat[mode, :step] * degrees[:step] ) * 1
        attachment_probabilities = ( homophily_mat[mode, :step] ) * 1

        if attachment_probabilities[:step].sum() == 0:
            #If there is zero preference to connect to anyone, just connect based on a random uniform distribution
            attachment_probabilities[:step] = 1
            pass

        if P_A == True:
            attachment_probabilities *= degrees[:step] #P.A

        attachment_probabilities[:step] /= attachment_probabilities[:step].sum()

        old_attached = np.random.choice(step, size = k_init, p = attachment_probabilities, replace = True) # may output repeated index

        for old_index in old_attached:
            if inter_group_pair_num_output:
                if modes[old_index] != mode:
                    inter_group_pair_num[mode] += 1

                    if inter_group_per_node:

                        external_links[old_index] += 1
                        external_links[step] += 1

#                     print(inter_group_pair_num)
            degrees[old_index] += 1
            importances[old_index] += 1

#         importances[:step] /= importances[:step].sum()
        importances[:step] *= eps


        #value *= eps_inv

        #if (step % 5000) == 0:
            #print(step)
        #print( degrees[ degrees > 0 ] )

        if history_step_size:
            if (step % history_step_size) == 0:
                degrees_history[ step // history_step_size - 1 ] = degrees


        if step == t_unbias:
            homophily_mat = symmetric_homophily_construct( 0.5, step_num, mode_num, modes )

    if inter_group_per_node:
        return degrees, modes, external_links, inter_group_pair_num

    if inter_group_pair_num_output:
        return degrees, modes, inter_group_pair_num

    if history_step_size:
        return degrees_history, history_steps, degrees, modes, fitnesses
    else:
        return degrees, modes, fitnesses

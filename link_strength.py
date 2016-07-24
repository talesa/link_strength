
# coding: utf-8

# In[160]:

# define variables

number_of_conditioned_variables = 2
number_of_states_in_variable = 3


# In[161]:

# 

number_of_states_in_conditioned_variable = []
for i in range(number_of_conditioned_variables):
    number_of_states_in_conditioned_variable.append(0)
    
link_strength = []
for i in range(number_of_conditioned_variables):
    link_strength.append(0)


# In[172]:

# define variables

number_of_states_in_conditioned_variable[0] = 3
number_of_states_in_conditioned_variable[1] = 3

link_strength[0] = 0
link_strength[1] = 0


# In[173]:

number_of_states_in_conditioned_variable


# In[174]:

# create placeholders for P(u'|u)

import numpy as np

P = list()

for i in range(number_of_conditioned_variables):
    P.append(np.matrix(np.zeros((number_of_states_in_conditioned_variable[i], number_of_states_in_conditioned_variable[i]))))


# In[175]:

for state_of_variable in range(number_of_states_in_variable):
    for conditioned_variable_no in range(number_of_conditioned_variables):
        for state_of_conditioned_variable in range(number_of_states_in_conditioned_variable[conditioned_variable_no]):
            r = state_of_variable
            c = state_of_conditioned_variable
            m_i = float(number_of_states_in_conditioned_variable[conditioned_variable_no])
            abs_eta_i = np.abs(link_strength[conditioned_variable_no])
            K = 1 - 1/m_i
            if r == c:
                P[conditioned_variable_no][r, c] = 1/m_i + abs_eta_i * K
            else:
                sum_term = 0
                for j in range(int(m_i)):
                    if j != r:
                        sum_term += 1/(j-r)**2
                P[conditioned_variable_no][r, c] = (abs_eta_i/(c-r)**2/sum_term + (1-abs_eta_i)/(m_i-1))*(1-1/m_i-abs_eta_i*K)


# In[176]:

def F(u_prime):
#     print enumerate(zip(link_strength, u_prime))
    res = sum([np.abs(eta_i)*H(u_prime_i, i) for i, (eta_i, u_prime_i) in enumerate(zip(link_strength, u_prime))])
    if sum(link_strength) != 0:
        res /= sum(link_strength)
#     print('Result of F(): {}, u_prime: {}'.format(res, u_prime))
    return res


# In[177]:

def H(u_prime_i, i):
    if link_strength[i] >= 0:
        return u_prime_i
    else:
        return -u_prime_i + 1 + number_of_states_in_conditioned_variable[i]


# In[178]:

import itertools

def generate_all_possible_u():
    groups = list()
    
    for i in range(number_of_conditioned_variables):
        group = [j for j in range(number_of_states_in_conditioned_variable[i])]
        groups.append(group)
        
    return list(itertools.product(*groups))


# In[186]:

P


# In[195]:

b = [1.5, 0.6, 1.5]


# In[210]:

# P(x|u)
# = sum over all u' such that F(u') = x

from functools import reduce
import operator

def Pr(x, u):
    sum_term = 0
    for u_prime in generate_all_possible_u():
        # do the rounding here
        
        F_value = F(u_prime)
        print('{} '.format(F_value))
        if F_value == x:
            weight = 1
        elif np.ceil(F_value) == x:
            weight = x - F_value
            F_value = x
            assert weight > 0
        elif np.floor(F_value) == x:
            weight = F_value - x
            F_value = x
            assert weight > 0
        else:
            continue
            
        prod_term = 1
        for i in range(number_of_conditioned_variables):
#             prod_term *= P[i][u_prime[i], u[i]]
            prod_term *= P[i][u[i], u_prime[i]]
        
        sum_term += prod_term + correction(weight, x)
        
    return sum_term


# In[211]:

def correction(weight2, x):
    c1 = np.max(np.abs(link_strength)) * weight2 + (1 - np.max(np.abs(link_strength))) * b[x] * weight2
    if np.max(np.abs(link_strength)) == 0:
        min_max = 1
    else:
        min_max = np.min(np.abs(link_strength))/np.max(np.abs(link_strength))
    c2 = min_max*c1 + weight2 * (1 - min_max)
    return c2


# In[212]:

for us in generate_all_possible_u():
    l = [Pr(i, us) for i in range(number_of_states_in_variable)]
    suml = sum(l)
    if suml != 0:
        l = [i/suml for i in l]
    print (l)


# In[ ]:




# In[ ]:




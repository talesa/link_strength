
# coding: utf-8

# In[114]:

number_of_conditioned_variables = 1


# In[115]:

number_of_states_in_conditioned_variable = []
for i in range(number_of_conditioned_variables):
    number_of_states_in_conditioned_variable.append(0)
    
link_strength = []
for i in range(number_of_conditioned_variables):
    link_strength.append(0)


# In[116]:

number_of_states_in_variable = 3


# In[117]:

number_of_states_in_conditioned_variable[0] = 3
# number_of_states_in_conditioned_variable[1] = 3


# In[118]:

link_strength[0] = 0.5
# link_strength[1] = 0.5


# In[119]:

P = np.matrix(np.zeros((reduce(lambda x, y: x * y, number_of_states_in_conditioned_variable), number_of_states_in_variable)))


# In[120]:

P


# In[ ]:

number_of_conditioned_variables 
states = number_of_conditioned_variables


# In[125]:

states = []
last_state = tuple([0 for i in range(number_of_conditioned_variables)])
for conditioned_variable_no in range(number_of_conditioned_variables):
    for state_of_conditioned_variable in range(number_of_states_in_conditioned_variable[conditioned_variable_no]):
        states.append(())


# In[126]:

last_state


# In[123]:

for state_of_variable in range(number_of_states_in_variable):
    for conditioned_variable_no in range(number_of_conditioned_variables):
        for state_of_conditioned_variable in range(number_of_states_in_conditioned_variable[conditioned_variable_no]):
            r = state_of_variable
            c = state_of_conditioned_variable
            m_i = float(number_of_states_in_conditioned_variable[conditioned_variable_no])
            abs_eta_i = np.abs(link_strength[conditioned_variable_no])
            K = 1 - 1/m_i
            if r == c:
                P[r, c] = 1/m_i + abs_eta_i * K
            else:
                sum_term = 0
                for j in range(int(m_i)):
                    if j != r:
                        sum_term += 1/(j-r)**2
                P[r, c] = (abs_eta_i/(c-r)**2/sum_term + (1-abs_eta_i)/(m_i-1))*(1-1/m_i-abs_eta_i*K)


# In[124]:

P


# In[ ]:




# In[ ]:





def del_h(a, b):
    if (a=='Fe' and b=='Co'):
        return(-2)
    if (a=='Co' and b=='Ni'):
        return(-1)
    if (a=='Ni' and b=='Fe'):
        return(-6)
    if (a=='Co' and b=='Fe'):
        return(-2)
    if (a=='Ni' and b=='Co'):
        return(-1)
    if (a=='Fe' and b=='Ni'):
        return(-6)



def vol(a):
    if (a=='Fe'):
        return(7.09)
    if (a=='Co'):
        return(6.7)
    if (a=='Ni'):
        return(6.6)




def C_S_ele(C_a,C_b ,a,b):

    C_S_elem = (C_a * (vol(a)** 0.667)) / ((C_a *( vol(a)**0.667)) + (C_b *( vol(b)**0.667)))

    return(C_S_elem)

def func(C_a,C_b, gamma , a ,b):


    f_ab = C_S_ele(C_b,C_a,a,b)*(1 + gamma*(np.power(C_S_ele(C_a,C_b,a,b) * C_S_ele(C_b,C_a,a,b), 2)))
    f_ba = C_S_ele(C_a,C_b,a,b)*(1 + gamma*(np.power(C_S_ele(C_a,C_b,a,b) * C_S_ele(C_b,C_a,a,b), 2)))

    del_h__ = C_a * C_b*(f_ab * (del_h(a, b)) + f_ba * (del_h(b, a)))

    return(del_h__)

import numpy as np
import pandas as p
a='Fe'
b='Co'
c='Ni'

X_a=(float)(input('concentration of Fe '))
X_b=(float)(input('concentration of Co '))
X_c=1-X_a-X_b

import sympy as sp
p1,p2,p3 = sp.symbols("p1,p2,p3", real=True)
C_a,C_b,C_c = sp.symbols("C_a,C_b,C_c", real=True)

del_h_mm=[]
del_h_mm.append(func(C_a,1-C_a,0 ,a,b))
del_h_mm.append(func(C_b,1-C_b,0,b,c))
del_h_mm.append(func(C_c,1-C_c,0,c,a))

eq=p1*del_h_mm[0]+p2*del_h_mm[1]+p3*del_h_mm[2]
print('minimizing the formation enthalpy')
#minimizing the equations using the derivative wrt p1 , p2 , p3
solutions=sp.solve((del_h_mm[0],del_h_mm[1],del_h_mm[2]),(C_a,C_b,C_c))
C_a=solutions[0]
C_b=solutions[1]
C_c=solutions[2]

solution_final = sp.solve((p1*C_a+p3*(1-C_c)-X_a,p1+p2+p3-1,p2*C_b+p1*(1-C_a)-X_b,p3*C_c+p2*(1-C_b)-X_c),(p1,p2,p3))

sp.simplify(solutions)


print('the minimized formation enthalpy is ' ,solution_final)






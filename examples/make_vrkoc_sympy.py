#Pyhon script to obtain an use symbolic order conditions for VRK
#similar to the notebook files.

from nviepy import vrk
from sympy import solve


method='BVRK' #'PVRK'  'BVRK'  'VRK'

matrixA='explicit' # 'explicit' 'diagimp' 'implicit'
p_max=3; s=3; 


arraylist=vrk.create_arrays(method,matrixA,s)

if method=='PVRK':
    b,c,A=arraylist
elif method=='BVRK':
    b,c,A,d=arraylist 
else:
    b,c,A,e,D=arraylist


Alleq=vrk.gen_order_conditions(p_max,arraylist)

solve(Alleq)

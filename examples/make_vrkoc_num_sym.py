#Make pyhon script to check if elementary_weight_str produces the same
# mathematical output as his symbolic conterpart.


from nviepy import volterra_rooted_3s as vrt

#xa='a'; method='PVRK'
xa='xa'; method='BVRK'
#xa='xa'; method='VRK'

p_max=6; s=2;
#Oss. if explicit must be p<s+2
#     if diagimp  or implicit: OK


f=open('oc_vrkSY.py','w')
 
f.write("from __future__ import division\n")
f.write("from nviepy import vrk\n")
f.write("from sympy import Eq\n")
f.write("from numpy import dot, ones\n\n")

f.write("s="+str(s)+";xa='"+xa+"';p_max="+str(p_max)+";method='"+method+"'\n\n")


f.write("arraylist=vrk.create_arrays(method,'diagimp',s)\n") 

if method=='PVRK':
    f.write("b,c,A=arraylist\n\n") 
elif method=='BVRK':
    f.write("b,c,A,d=arraylist\n\n") 
else:
    f.write("b,c,A,e,D=arraylist\n\n") 


f.write("Alleq=vrk.gen_order_conditions(p_max,arraylist)\n") 

f.write("\ncheck=[]\n") 
ioc=0
p=0
many_forest=vrt.list_trees(p_max,'upto',xa=xa)
for forest in many_forest:
    p+=1
    f.write("\n    # order "+str(p)+" conditions:\n")
    for tree in forest:
        vrk_oc=vrt.elementary_weight_str(tree,method,style='python')
        rhs =vrt.right_hand_side_str(tree,style='python')
        f.write("check.append(Alleq["+str(ioc)+"].expand()==Eq(dot(b.T,"+vrk_oc+")[0,0],"+rhs+").expand())\n")
        #f.write("Alleq["+str(ioc)+"].expand()==Eq(dot(b.T,"+vrk_oc+")[0,0],"+rhs+").expand()\n")
        ioc+=1

f.write("\nprint check")

f.close()

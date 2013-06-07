#Make Octave script to evaluate order conditions for VRK, with a given formula

from nviepy import volterra_rooted_3s as vrt

p_max=5

#s=3; xa='xa'          s=2; xa='a'  
s=3; xa='xa' 

f=open('oc_vrk.m','w')
 
 
f.write("s="+str(s)+" # number of stages \n")

#PVRK s=2
f.write("A=[1/4, (3-2*sqrt(3))/12;(3+2*sqrt(3))/12, 1/4];\n") 
f.write("b=[1/2;1/2];    c=[(3-sqrt(3))/6;(3+sqrt(3))/6]; \n")

#BVRK s=3
f.write("A=[0,0,0;1,0,0;1/9,2/9,0];d=[1/2;1;1]; \n") 
f.write("b=[0;1/4;3/4];  c=[0;1;1/3];    \n")



for p in range(1,p_max):
    ioc=1
    f.write("\n    # order "+str(p)+" conditions:\n")
    forest = vrt.list_trees(p,xa=xa)
    for tree in forest:
        vrk_oc=vrt.elementary_weight_str(tree,method='BVRK',style='octave')
        rhs =vrt.right_hand_side_str(tree,style='python')
        f.write("coneq("+str(ioc)+")=dot(b',"+vrk_oc+")-"+rhs+"\n")
        ioc+=1

f.write("\n")
f.close()

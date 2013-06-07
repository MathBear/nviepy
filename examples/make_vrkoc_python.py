#Make pyhon script to evaluate order conditions for VRK, with a given formula

from nviepy import volterra_rooted_3s as vrt

p_max=5


#xa='a'; method='PVRK'
xa='xa'; method='BVRK'

f=open('oc_vrk.py','w')
 
f.write("from __future__ import division\n")
f.write("from numpy import sqrt, array, ones, dot, zeros\n\n")


if xa=='xa':
    #BVRK s=p=3
    s=3
    f.write("s="+str(s)+" # number of stages \n")
    f.write("A=array([[0,0,0],[1,0,0],[1/9,2/9,0] ]);d=array([[1/2],[1],[1]]) \n") 
    f.write("b=array([[0],[1/4],[3/4]]);  c=array([[0],[1],[1/3]])    \n")
    
    #BVRK s=5 p=4
    s=5
    f.write("s="+str(s)+" # number of stages \n")
    f.write("c2=(3-sqrt(3))/8; a52=(sqrt(3)-2)/(12*c2); a42=(2544-807*sqrt(3))/(13754*c2); a43=(1245*sqrt(3)-90)/6877 \n")
    f.write("A=array([ [0,0,0,0,0],[c2,0,0,0,0],[(3-sqrt(3))/6+a52 ,-a52 ,0,0,0],")
    f.write(" [(2781-647*sqrt(3))/6877-a42,a42,a43,0,0],[(-3+2*sqrt(3))/9-a52,a52,1/5,(57-5*sqrt(3))/90,0]  ])\n")
    f.write("d=array([[(3-sqrt(3))/4],[(3-sqrt(3))/4-c2],[1],[(57+5*sqrt(3))/92],[1]]) \n")
    f.write("c=array([[0],[c2],[(3-sqrt(3))/6],[(9+2*sqrt(3))/23],[(3+sqrt(3))/6]]) \n")
    f.write("b=array([[0],[0],[1/2],[0],[1/2]]) \n")
    
else:
    #PVRK s=2
    s=2
    f.write("s="+str(s)+" # number of stages \n")
    f.write("A=array([[1/4, (3-2*sqrt(3))/12],[(3+2*sqrt(3))/12, 1/4]])\n")
    f.write("b=array([[1/2],[1/2]]);    c=array([[(3-sqrt(3))/6],[(3+sqrt(3))/6]]) \n")

p=0
list_of_forests=vrt.list_trees(p_max,'upto',xa=xa)
for forest in list_of_forests:
    ioc=0
    p+=1
    f.write("\n    # order "+str(p)+" conditions:\n")
    f.write("coneq = zeros(("+str(len(forest))+"))\n")
    for tree in forest:
        vrk_oc=vrt.elementary_weight_str(tree,method='BVRK',style='python')
        rhs =vrt.right_hand_side_str(tree,style='python')
        f.write("coneq["+str(ioc)+"]=dot(b.T,"+vrk_oc+")-"+rhs+"\n")
        ioc+=1
    f.write("print 'Estimate for order "+str(p)+"'\n")
    f.write("print coneq\n")


f.write("\n")
f.close()

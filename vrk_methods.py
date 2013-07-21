from nviepy import volterra_rooted_3s as vrt
from sympy import Symbol, Eq, Matrix
from numpy import array, dot, size, sum


#=====================================================
# Functions for creating Matrix A
#=====================================================
# c for character, s for string, and d for a integer decimal value
def _exp(M,i,j):
     if i > j:
             return Symbol('%c_%d%d'%(M,i,j)) 
     else:
             return 0

def _imp(M,i,j):
    return Symbol('%c_%d%d'%(M,i,j)) 

def _diagimp(M,i,j):
     if i >= j:
             return Symbol('%c_%d%d'%(M,i,j)) 
     else:
             return 0


#=====================================================
# Functions for symbolic order conditions
#=====================================================
def create_arrays(vrktype,Amatrix,stage):
    """
        Constructs a list () of matrix and vector containing the symbols of a VRK method.
        
        INPUT:
            
        - vrktype -- select which class of VRK formulas (VRK, BVRK or PVRK) to use, 
                    it determines the number of output arrays.
        - Amatrix -- select the structure of matrix A, they are:
                    'explicit', 'diagonally_implicit' or 'implicit'.
        - stage   -- number of stages, determines the length of the output arrays.
        
        OUTPUT: Is a list of 3,4 or 5 arrays, whether the vrktype is VRK, BVRK or PVRK .
                The output is a Tuple [] containing several arrays of type:
                <type 'numpy.ndarray'>.

        **Examples**:

            >>> from nviepy import vrk
            >>> vrk.create_arrays('PVRK','implicit',2)
                (array([[b_1],
                   [b_2]], dtype=object), array([[c_1],
                   [c_2]], dtype=object), array([[a_11, a_12],
                   [a_21, a_22]], dtype=object))
            >>> vrk.create_arrays('BVRK','explicit',2)
                (array([[b_1],
                        [b_2]], dtype=object), array([[c_1],
                        [c_2]], dtype=object), array([[0, 0],
                        [a_21, 0]], dtype=object), array([[d_1],
                        [d_2]], dtype=object))
            
    """
    
    Range=range(1,stage+1)
    
    b = array([[Symbol('b_%d' % i)] for i in Range])
    
    c = array([[Symbol('c_%d' % i)] for i in Range])
    
    if Amatrix == 'explicit':
        A=array([[_exp('a',i,j) for j in Range] for i in Range]) 
    elif (Amatrix == 'diagonally_implicit' or Amatrix=='diagimp'):
        A=array([[_diagimp('a',i,j) for j in Range] for i in Range])
    else:
        A=array([[_imp('a',i,j) for j in Range] for i in Range])
    
    if (vrktype=='pouzet' or vrktype=='PVRK'):
        return b, c, A
    elif (vrktype=='beltyukov' or vrktype=='BVRK'):
        d=array([[Symbol('d_%d'% i)] for i in Range])
        return b, c, A, d
    else:
        e=array([[Symbol('e_%d'% i)] for i in Range])
        D=array([[Symbol('d_%d%d'%(i,j)) for j in Range] for i in Range])
        return b, c, A, e, D
    

        
def gen_order_conditions(order,arraylist,RSC=True):
    """ Generates a list [] of order conditions of a given 'order', using 'arraylist'.
        
        INPUT:
            
        - order     -- order of the formula.
        - arraylist -- input that comes from the output of create_arrays.
        - RSC       -- choose whether to add to the list of 
                        equations the Row Sum Conditions.
        
        OUTPUT: list of equations.

        **Examples**:

            >>> from nviepy import vrk
            
            >>> arraylist=vrk.create_arrays('PVRK','explicit',2)
            >>> vrk.gen_order_conditions(2,arraylist)
            [b_1 + b_2 == 1, b_1*c_1 + b_2*c_2 == 1/2, c_1 == 0, c_2 == a_21]
            >>> vrk.gen_order_conditions(3,arraylist)
            [b_1 + b_2 == 1, b_1*c_1 + b_2*c_2 == 1/2, a_21*b_2*c_1 == 1/6, 
             b_1*c_1**2 + b_2*c_2**2 == 1/3, c_1 == 0, c_2 == a_21]
             
            >>> arraylist=vrk.create_arrays('BVRK','explicit',2)
            >>> vrk.gen_order_conditions(2,arraylist,'n')
            [b_1 + b_2 == 1, b_1*d_1 + b_2*d_2 == 1, b_1*c_1 + b_2*c_2 == 1/2]
            >>> vrk.gen_order_conditions(3,arraylist,'n')
            [b_1 + b_2 == 1, b_1*d_1 + b_2*d_2 == 1, b_1*c_1 + b_2*c_2 == 1/2,
             b_1*d_1**2 + b_2*d_2**2 == 1, a_21*b_2*d_1 == 1/3, a_21*b_2*c_1 == 1/6,
             b_1*c_1**2 + b_2*c_2**2 == 1/3, b_1*c_1*d_1 + b_2*c_2*d_2 == 1/2]
    """
    
    numarray=size(arraylist,0)
    if not 3 <= numarray <= 5:
        raise ValueError('arraylist must be between 3 and 5')
    
    if not size(arraylist[0],0)==size(arraylist[1],0)==size(arraylist[2],0):
        raise ValueError('First dimension does not match')
    
    if not size(arraylist[0],1)==size(arraylist[1],1)==1:
        raise ValueError('First two elements are not column vector')
        
    if not size(arraylist[2],0)==size(arraylist[2],1):
        raise ValueError('Third element is not a square matrix')
    
    s=size(arraylist[0],0)
    b = array(arraylist[0])
    if numarray==3:
        xa='a'
        method='PVRK'
    elif numarray==4:
        xa='xa'
        method='BVRK'
    elif numarray==5:
        xa='xa'
        method='VRK'
    else: raise TypeError('Unexpected Error!!! :P')
    
    
    Alleq=[]
    many_forest=vrt.list_trees(order,'upto',xa)
    for forest in many_forest:
        for tree in forest:
            vrk_oc=vrt.elementary_weight(tree,s,arraylist[1:],method)
            rhs =vrt.right_hand_side(tree)
            Alleq.append(Eq(dot(b.T,vrk_oc)[0,0],rhs).expand())
    
    if RSC:
        #Finally add the Row Sum Conditions: c=Au. There are s of them.
        c=arraylist[1]
        Au=sum(arraylist[2], axis=1)
        for i in range(0,s):
            Alleq.append(Eq(c.T[0][i],Au[i])) #Alleq.append(Eq((c.T[0]-Au)[i]))
        
    
    return Alleq
    


  
def substitute(arraylist,listsubs, executexpand=True):
    """ Function useful to substitute many values in many arrays.
    
        Substitutes all the conditions listed in 'listsubs' in all the 
        arrays contained in 'arraylist'.
        
        INPUT:
            
        - arraylist -- list containing multiple arrays.
        - listsubs  -- dictionary containing variables in arrraylist to substitue.
        - executexpand -- choose wether to expand the output or not (False or True);
                        if expanded the output can be less readable, but unquely defined.
        
        OUTPUT: whatever array or matrix in input, gives a Matrix in output.

        **Examples**:

            >>> from nviepy import vrk
            >>> arraylist=vrk.create_arrays('PVRK','explicit',2)
            >>> b,c,A=arraylist
            >>> vrk.substitute(arraylist,{b[0][0]:b[1][0], A[1][0]: c[0][0]})
            ([b_2]
            [b_2], [c_1]
            [c_2], [  0, 0]
            [c_1, 0])
    """
    
    matrixlist=()
    if executexpand==True:
        for arrayelement in arraylist:
            matrixlist+=(Matrix(arrayelement).subs(listsubs).expand(),)
    else:
        for arrayelement in arraylist:
            matrixlist+=(Matrix(arrayelement).subs(listsubs),)
    
    return matrixlist
    
    
    
def SymPowenest(arraylist):
    """ Function useful to apply powdenest to a list of arrays.
    
        uses powdenest from sympy on all the elementrs of the
        arrays contained in 'arraylist'.
        
        INPUT:
            
        - arraylist -- list containing multiple arrays.
        
        OUTPUT: whatever array or matrix in input, gives a Matrix in output.

    """
    from sympy import powdenest
    
    matrixlist=()
    for arrayelement in arraylist:
        #simplify(element,ratio=10)
        newarrayelement=[[ powdenest(element,force=True) for element in vector] for vector in array(arrayelement)]
        matrixlist+=(Matrix(newarrayelement),)
    
    return matrixlist


    
def subEqns(eqns,listsubs, executexpand=True):
    """ Function usefoul to substitute many values in many arrays.
    
        Substitutes all the conditions listed in 'listsubs' in all the 
        arrays contained in 'arraylist'.
        
        INPUT:
            
        - arraylist -- list containing multiple arrays.
        - listsubs  -- dictionary containing variables in arrraylist to substitue.
        - executexpand -- choose wether to expand the output or not (False or True);
                        if expanded the output can be less readable, but unquely defined.
        
        OUTPUT: whatever array or matrix in input, gives a Matrix in output.

        **Examples**:

            >>> from nviepy import vrk
            >>> arraylist=vrk.create_arrays('PVRK','explicit',2)
            >>> b,c,A=arraylist
            >>> vrk.substitute(arraylist,{b[0][0]:b[1][0], A[1][0]: c[0][0]})
            ([b_2]
            [b_2], [c_1]
            [c_2], [  0, 0]
            [c_1, 0])
    """
    
    print size(eqns)
    neweqns=[]
    for eq in eqns:
        print eq
        try:
            neweq=eq.subs(listsubs).expand()
            print neweq
            print (not neweq)
            if not neweq:    #ater expanding neweq may become True
                print '\n entrato'
                print neweq.expand()
                neweqns.append(neweq)
        
        except AttributeError: #AttributeError: 'bool' object has no attribute 'expand' 
            print '\n ErrorAtribute \n|'
            pass
        except: 
            print "Unexpected error!"
    
    print size(neweqns)
    print size(list(set(neweqns)))
    
    return list(set(neweqns)) # removes possible duplicate
    
    
    
    
    
def gen_Sharp_order_conditions(order,arraylist):
    """ Generates order condition using Sharp theory.
        
        INPUT:
            
        - order     -- order of the formula.
        - arraylist -- input that comes from the output of create_arrays.
        
        OUTPUT: list of equations.

        **Examples**:

            >>> from nviepy import vrk
            
            
    """
    from sympy import eye, ones
    
    numarray=size(arraylist,0)
    if numarray != 4:
        raise ValueError('arraylist must be 4')
    
    if not size(arraylist[0],0)==size(arraylist[1],0)==size(arraylist[2],0):
        raise ValueError('First dimension does not match')
    
    if not size(arraylist[0],1)==size(arraylist[1],1)==1:
        raise ValueError('First two elements are not column vector')
        
    if not size(arraylist[2],0)==size(arraylist[2],1):
        raise ValueError('Third element is not a square matrix')
    
    s=size(arraylist[0],0)
    u = array(ones((s,1)))
    I = array(eye(s))
    bT= array(arraylist[0]).T
    #b = array(arraylist[0])
    C = array(arraylist[1])*I
    #c = array(arraylist[1])
    A = array(arraylist[2])
    D = array(arraylist[3])*I
    #d = array(arraylist[3])
    
    L=[bT,dot(bT,(2*C-I)),dot(bT,(2*A-I)),dot(bT,(D-I))]
    #L=[bT,(b*(2*c-u)).T,dot(bT,(2*A-I)),(b*(d-u)).T]
    
    #q_20=dot((dot(A,C)-C**2/2),u)
    #q_20=dot(A*c,u)-c**2/2
    
    q_20=array([[0,-C[1][1]**2/2,0,0,0]]).T# Sharp Hypotesys
    #q_20=array([[0,-c[1]**2/2,0,0,0]]).T
    
    q_21=dot((dot(A,D)-C**2),u)
    #q_21=dot(A*d,u)-c**2
    
    R=[dot((D-I)*(2*D-I),u), dot((D-I)*(2*C-D),u), q_20, q_21, dot(6*C**2-6*C-I,u)]
    #R=[(d-u)*(2*d-u), (d-u)*(2*c-d), q_20, q_21, 6*c**2-6*c-u]
    
    Lu=dot(L,u)
    
    
    Alleq=[]
    
    #
    #Finally add the Row Sum Conditions: c=Au. There are s of them.
    c=arraylist[1]
    Au= sum(A, axis=1)
    for i in range(0,s):
        Alleq.append(Eq(c.T[0][i],Au[i])) #Alleq.append(Eq((c.T[0]-Au)[i]))
        
    
    Alleq.append(Eq(Lu[0,0,0],1).expand())
    for vector in Lu[1:]:
        Alleq.append(Eq(vector[0,0],0).expand())
    
    
    #Orthogoanlity
    for row_vector in L:
        for col_vector in R:
            Alleq.append(Eq(dot(row_vector,col_vector)[0,0],0).expand())
    
    q_20=dot((dot(A,C)-C**2/2),u)
    # Sharp Hypotesys
    for i in range(2,s):
        Alleq.append(Eq(q_20[i,0],0).expand())    
    
    
    return Alleq
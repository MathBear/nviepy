import numpy as np
import matplotlib.pyplot as pl
#from itertools import combinations_with_replacement, product
#from utils import  gen_partitions_ms
#from sympy import Rational, eye, ones
#from strmanip import get_substring, open_to_close, countlett, getint


#=====================================================
class RootedTree(str):
#=====================================================
    r"""
        A rooted tree is a directed acyclic graph with one node, which
        has no incoming edges, designated as the root.
        Rooted trees are useful for analyzing the order conditions of 
        multi-stage numerical VIEs solvers, such as Volterra Runge-Kutta methods.

        The trees are represented as strings, using one of the notations
        introduced by Butcher.  The characters 'a' and 'x' are used 
        to represent a vertex, and braces '{ }' are used instead of brackets 
        '[ ]' to indicate that everything inside the braces
        is joined to a single parent node.  Thus the trees up to order 3 are:

        'a';
        
        '{x}', ' {a}';
        
        '{x^2}', '{{x}}', '{{a}}', '{aa}', '{xa}'.
        

        These can be generated using the function list_trees(), which
        returns a list of all trees of a given order::

            >>> from nviepy import *
            >>> for p in range(4): print vrt.list_trees(p)
            ['']
            ['a']
            ['{x}', '{a}']
            ['{x^2}', '{{x}}', '{{a}}', '{aa}', '{xa}']
            
        Note that the tree of order 0 is indicated by an empty string.
        
        If the tree contains an edge from vertex A to vertex B, vertex
        B is said to be a child of vertex A.
        A vertex with no children is referred to as a leaf.

        .. warning::

            One important convention is assumed in the code; namely, that at each 
            level, 'x' leaves are listed first, then 'a' leaves, finally any 
            other subtrees, and only if there are $n$ 'x' leaves, we write 'x^n'.
        
        .. note::

            Powers cannot be used for subtrees; thus
            '{{a}{a}}'
            is valid, while
            '{{a}^2}' 
            is not.
             

        **Examples**::

            >>> from nviepy import volterra_rooted_3s as vrt
            >>> tree=vrt.RootedTree('{x^2{a{a}}{a}}')
            >>> tree.order()
            9
            >>> tree.density()
            112

        We can generate Python code to evaluate the elementary weight
        corresponding to a given tree for a given class of methods::

            >>> vrt.elementary_weight_str(tree)
            'd**2*dot(A,c*dot(A,c))*dot(A,c)'

        **References**:  
            #. [BrHaNo1982]_
            #. [Harier81]_
            #. [brunner1986numerical]_
    """
    def __init__(self,strg):
        if any([strg[i] not in '{}ax^1234567890' for i in range(len(strg))]):
            raise Exception('Not a valid rooted tree string (illegal character)')
        op,cl=strg.count('{'),strg.count('}')
        if op!=cl or (op+cl>0 and (strg[0]!='{' or strg[-1]!='}')):
            raise Exception('Not a valid rooted tree string')

    def order(self):
        """
        The order of a rooted tree, denoted $\\rho(t)$, is the number
        of vertexes in the tree.

        **Examples**::

            >>> from nviepy import volterra_rooted_3s as vrt
            >>> tree=vrt.RootedTree('{a^2{a{x}}}')
            >>> tree.order()
            7
        """
        from strmanip import countlett
        if self=='a': return 1
        if self=='x': return 1 #useful when called from 'density'
        if self=='':  return 0
        
        # Count brackets
        r=self.count('{')
        
        # Count x
        r+=countlett('x',self)

        # Count a
        r+=countlett('a',self)
        
        return r

    def density(self):
        r"""
        The density of a rooted tree, denoted by $\\gamma(t)$,
        is the product of the orders of the subtrees.

        **Examples**::

            >>> from nviepy import volterra_rooted_3s as vrt
            >>> tree=vrt.RootedTree('{a^2{x{a}}}')
            >>> tree.density()
            42

        **Reference**: 

            - [BrHaNo1982]_ p. 156
        """
        nxleaves,naleaves,subtrees=self._parse_subtrees()
        gamma=self.order()-nxleaves
        for tree in subtrees:
            gamma*=tree.density()
        return gamma

    def plot(self,nrows=1,ncols=1,iplot=1,ttitle='',order=1):
        """
            Plots the rooted tree.

            *INPUT*: (optional)
                * nrows, ncols -- number of rows and columns of subplots
                  in the figure
                * iplot        -- index of the subplot in which to plot
                  this tree

            These are only necessary if plotting more than one tree
            in a single figure using subplot.

            *OUTPUT*: <matplotlib.figure.Figure>.

            The plot is created recursively by
            plotting the root, parsing the subtrees, plotting the 
            subtrees' roots, and calling _plot_subtree on each child
        """
        if iplot==1: pl.clf()
        splot=pl.subplot(nrows,ncols,iplot)
        pl.hold(True)
        if self!='':
            pl.scatter([0],[0],s=120, c='w')
            pl.scatter([0],[0],s=20, c='k')
            if self!='a': self._plot_subtree(0,0,0)
        
        ypos=(4-nrows)*0.2  #1 2 3  -4-   5 6 7
        pl.text(0,ypos-0.8, ttitle, ha="center", va="bottom",size='small')#size=fs  size='small'
        pl.hold(False)
        
        pl.xticks([])
        pl.yticks([])
        pl.ylim((-0.8, order-0.6))
        semiXlength=order/2.+0.1   # the output is a square
        pl.xlim((-semiXlength, semiXlength))
        pl.axis('off')
        return splot
        
        
    def _plot_subtree(self,xroot,yroot,xwidth):
        """
            Recursively plots subtrees.  Should only be called from plot().

            INPUT:
                xroot, yroot -- coordinates at which root of this subtree 
                                is plotted
                xwidth -- width in which this subtree must fit, in order
                          to avoid possibly overlapping with others
                          If it is 0, the tree has been growing straight above.
        """
        ychild=yroot+1
        nxleaves,naleaves,subtrees=self._parse_subtrees()
        nleaves=nxleaves+naleaves
        nchildren=nleaves+len(subtrees)
        
        assert nchildren>=1, 'Unexpected Error!'
        if nchildren==2: 
            dist=0.5  # prettier output
        else:
            dist=(nchildren-1)/4.
        
        CONST=1.0   # for dense trees, e.g. tree.density()/tree.order()>1000  
                    # this number should be smaller then one.
        if dist > 0.2: # at least dist > 0
            if xwidth==0:
                xwidth=CONST*dist
            else:
                dist*=CONST*xwidth
                xwidth-=dist  #xwidth-dist
        
        xchild=np.linspace(xroot-dist,xroot+dist,nchildren)
        for i in range(nchildren):
            pl.plot([xroot,xchild[i]],[yroot,ychild],'-k')
            if i>nleaves-1:
                subtrees[i-nleaves]._plot_subtree(xchild[i],ychild,xwidth)
        #pl.scatter(xchild,ychild*np.ones(nchildren),s=20,c=['w']*nxleaves+['k']*(nchildren-nxleaves))
        pl.plot(xchild[:nxleaves],ychild*np.ones(nxleaves),'wo')
        pl.plot(xchild[nxleaves:],ychild*np.ones(nchildren-nxleaves),'ko')
        
    def _parse_subtrees(self):
        """ 
            Returns the number of 'x' and 'a' leaves and a list of the subtrees,
            for a given rooted tree.

            OUTPUT:
                nXleaves  -- number of 'x' leaves attached directly to the root
                nAleaves  -- number of 'a' leaves attached directly to the root
                subtrees -- list of subtrees attached to the root

            The method can be thought of as returning what remains if the
            root of the tree is removed.  For efficiency, instead of
            returning possibly many copies of 'a' and 'x', the leaves are just
            returned as a number.
        """
        from strmanip import get_substring, open_to_close, countlett
        if str(self)=='a' or str(self)=='': return 0,0,[]
        
        # current level        
        cl=self[1:self.find('{',1)+1]+self[self.rfind('}',0,-1)+1:-1]
        
        # Count 'x' leaves at current level
        nXleaves=countlett('x',cl)
        
        # Count 'a' leaves at current level
        nAleaves=countlett('a',cl)
        
        subtrees=[]
        pos=self.find('{',1)
        while pos!=-1:
            subtrees.append(RootedTree(get_substring(self,pos)))
            pos=open_to_close(self,pos)
            pos=self.find('{',pos+1)
        
        return nXleaves,nAleaves,subtrees

    def __mul__(self,tree2):
        """ 
            Returns Butcher's product: t*u is the tree obtained by
            attaching the root of u as a child to the root of t. 
        """
        from strmanip import getint
        if self=='a': return RootedTree('{'+tree2+'}')
        if self=='x': 
            print "Error!" 
            return False
        if tree2=='x':  # We're just adding an x leaf to self
            nxleaves,naleaves,subtrees=self._parse_subtrees()
            if nxleaves==0: return RootedTree(self[0]+'x'+self[1:])
            if nxleaves==1: return RootedTree(self[0]+'x^2'+self[2:])
            if nxleaves>1:
                n = getint(self[3:])
                return RootedTree(self[0:3]+str(n+1)+self[(3+len(str(n))):])
        else: return RootedTree(self[:-1]+tree2+'}') # tree2 wasn't just 'x'
#=====================================================
#End of RootedTree class
#=====================================================

#=====================================================
def plot_all_trees(p,subtitle=True,TypeOfTree='xa'):
#=====================================================
    """ Plots all rooted trees of order p """
    MaxSubplot=49 #maximum number of trees in one figure= 48+1
    
    forest=list_trees(p,xa=TypeOfTree)
    nplots=len(forest)
    nfigure=nplots//MaxSubplot+1
    nplot4fig=nplots//nfigure
    nrows=int(np.ceil(np.sqrt(float(nplot4fig))))
    ncols=int(np.round(np.sqrt(float(nplot4fig))))
    
    if nfigure>1:
        fig=[]
        fig.append(pl.figure(1,facecolor='white'))
    else:
        fig=pl.figure(1,facecolor='white')
    ithfig=0
    ith3=0
    for tree in forest:
        if ith3==nplot4fig and ithfig<nfigure-1:
            pl.suptitle('Rooted Trees of Order '+str(p),fontsize='large')
            pl.tight_layout()
            #close the previous figure and open a new one
            ith3=0
            ithfig+=1
            fig.append(pl.figure(ithfig+1,facecolor='white'))
        if subtitle: ttitle=tree
        else: ttitle=''
        ith3+=1
        #could use fig[ithfig].add_subplot(111)
        tree.plot(nrows,ncols,ith3,ttitle,p) #ith3= forest.index(tree)+1
        
    pl.suptitle('Rooted Trees of Order '+str(p),fontsize='large')
    pl.tight_layout()
    return fig

#=====================================================
def list_trees(p,ind='all',xa='xa'):
#=====================================================
    r""" 
    Returns rooted trees of order p.
    This algorithm generates all rooted trees for a given order p,
    
    We use the root-oriented recursive approach that was chosen by Sofroniou 
    in his Mathematica package Butcher.m.
    He generates all trees t of orde p by first, listing all integer partitions 
    $p-1=p_1 +...+p_n, n=1,...,p-1$, and next, 
    setting $t=[t_1,...,t_n]$ for all trees $t_1,...,t_n$ of order
    $#t_1=p_1<p,...,,#t_n=p_n<p.These trees have already been generated by the recursion.
    
    INPUT: 
    
        - p   -- order of trees desired
        - ind -- if given, as a number returns a single tree corresponding
                           to that index. Not very useful since 
                           the ordering isn't obvious.
                 If is given 'upto', in output is given an array which contains p arrays,
                 each i-th array has trees of i-th order.
        - xa  -- chose whether to have trees containing both flags 'x' and 'a',
                 or only the ones with 'a', this is useful for Pouzet's formulas or
                 trees of Runge-Kutta methods for ODEs.
                
    OUTPUT: list of all trees of order p. 
            If 'ind' is provided than we can have just one tree, 
            or all of them up to order p.
    
    **Examples**:

    Produce number of order condition of a BVRK formula ::
        
        >>> from nviepy import volterra_rooted_3s as vrt
        >>> for i in range(1,11): 
        ...     forest=vrt.list_trees(i)
        ...     print len(forest)
        1
        2
        5
        13
        37
        108
        332
        1042
        3360
        11019


    Produce column of Butcher's Table 302(I) or number of order condition
    of a PVRK formula::
        
        >>> from nviepy import volterra_rooted_3s as vrt
        >>> for i in range(1,11): 
        ...     forest=vrt.list_trees(i,xa='a')
        ...     print len(forest)
        1
        1
        2
        4
        9
        20
        48
        115
        286
        719
        
    **Reference**:
            #. [BrHaNo1982]_
            #. [butcher2008]_
            #. [Sofroniou1994]_
    """

    if p==0: return [RootedTree('')]
    
    W=[[],[]]
    W[0].append(RootedTree("a"))
    # Now add the order two Volterra trees:
    if xa=='xa': W[1].append(RootedTree("{x}"))
    W[1].append(RootedTree("{a}"))
    
    if xa=='xa':
        for i in range(2,p):
            # Construct W[i]
            ps=_powerString("x",i,powchar="^") # Construct the x 'bloom'
            W.append([RootedTree("{"+ps+"}")])
            for k in range(0,i):
                # Start from 'little' x bloom
                ps=_powerString("x",k,powchar="^")
                W=_genSubTrees(i,k,W,ps)
    else:
        for i in range(2,p):
            W.append([])
            W=_genSubTrees(i,0,W,"")
        
    
    if ind=='all': return W[p-1]
    elif ind=='upto': return W
    else: return W[p-1][ind]
    

def _genSubTrees(i,k,W,ps):
    """Generates stings of trees of order i, given the ones before
       with the root-oriented technique."""
       
    from itertools import combinations_with_replacement, product 
    from utils import  gen_partitions_ms
    
    for r in gen_partitions_ms(i-k):
        Rp=[] # R products
        for l in r:
            Rt=[] # R temporary
            cwr=combinations_with_replacement(W[l-1],r[l])
            for elem in cwr:
                Rt.append(reduce(lambda x, y: x+y, elem))
            Rp.append(Rt)
        Rp=product(*Rp) # Cartesian Product
        for newsubtree in Rp:
            W[i].append(RootedTree("{"+ps+reduce(lambda x, y: x+y, newsubtree)+"}"))
    return W


def _powerString(s,npow,powchar="**",trailchar=''):
    """Raise string s to power npow with additional formatting."""
    if npow==0:
        return ""
    else:
        if npow==1:
            return s+trailchar
        else:
            return s+powchar+str(npow)+trailchar



#=====================================================
# Functions on trees
#=====================================================
 
def right_hand_side(tree):
    """ 
    Gives the $1/\\gamma(t)$ value as a Rational number.
    The output belongs to:
        <class 'sympy.core.numbers.Rational'>.

    **Examples**::

        >>> from nviepy import volterra_rooted_3s as vrt
        >>> tree=vrt.RootedTree('{a^2{a{a}}}')
        >>> vrt.right_hand_side(tree)
        1/56

    **Reference**: 

        [BrHaNo1982]_
    """
    from sympy import Rational
    return Rational(1,(tree.density()))

def right_hand_side_str(tree,style='python'):
    """ 
    Gives the $1/\\gamma(t)$ value as a string, it can be latex or python style.

    **Examples**::

        >>> from nviepy import volterra_rooted_3s as vrt
        >>> tree=vrt.RootedTree('{a^2{a{a}}}')
        >>> vrt.right_hand_side_str(tree)
        '1/56'

    **Reference**: 

        [BrHaNo1982]_
    """
    from sympy import Rational
    dens=tree.density()
    if style=='latex': 
        if dens==1:return "1"
        else: return r"\frac{1}{"+str(dens)+"}" 
    else: 
        return str(Rational(1,(dens)))



def elementary_weight(tree,s,arrays,method):
    """
        Constructs elementary weights for a Volterra Runge-Kutta method,
        supposing the row sum condition.
        The output needs to be multiplied by b^T and equated to LHS
        to obtain the order condition.
        
        It is used by gen_order_conditions in vrk_methods
        
        INPUT:
            
        - tree   -- input tree, must be a RootedTree.
        - s      -- number of stages
        - arrays -- it depends on the method input, is a list containing two
                    or three arrays, the order should be: c,A,e,D (if VRK)
                    c,A,d (if BVRK) or c,A if (if PVRK)
        - method -- select which type of method to use: VRK, BVRK or if a PVRK
                    method is wanted, the three must be created with the xa='a' flag.
        
        OUTPUT: it is a column vector belonging to 
                <type 'numpy.ndarray'>
                e.g. like array([[],...,[]], dtype=object)

    """
    from sympy import eye, ones
    if tree=='': return ''    
    
    u=np.array(ones((s,1)))   #np.ones((s, 1), dtype=np.int)
    if tree=='a': return u # Matrix(u)

    I=np.array(eye(s))    # np.eye(s, dtype=np.int)
    ew=u.copy()
    
    c=arrays[0]
    A=arrays[1]
    if method=='VRK':
        d=arrays[2]
        D=arrays[3]
    elif method=='BVRK':
        d=arrays[2]
    else: #method=='PVRK'
        d=arrays[0]
    
    
    
    nx,na,subtrees=tree._parse_subtrees()
    ew*=c**na #na and nx can also be zero, then we have a vector of ones.
    ew*=d**nx
    
    # Two curly bracket contain at least a symbol, 'a' or 'x', i.e {} is not a leaf
    if len(subtrees)>0:
        for subtree in subtrees:
            if method=='VRK':
                ew=ew*_elem_weight_sub3(subtree,method,c,D,A,I,u)
            else:
                ew=ew*_elem_weight_sub3(subtree,method,c,d,A,I,u)
    
    return ew #returns a column


def _elem_weight_sub3(tree,method,c,D,A,I,u):
    
    ews=I.copy() #ews is a matrix
    
    nx,na,subtrees=tree._parse_subtrees()
    lens3s=len(subtrees)
    
    ews*=I*c**na  # ews=np.dot(ews,I*c**na)
    
    #we do not multiply by A prior entering _elem_weight_sub3, so we multiply now.
    if method!='VRK':  # D is a column
        ews*=I*D**nx
        ews=np.dot(A,ews)
    else:
        ews=np.dot(A*D**nx,ews)
    # ews is still a matrix
    
    
    # Two curly bracket contain al least a symbol, 'a' or 'x', i.e {} is not a leaf
    if lens3s>0:
        ewsi=I.copy() 
        for sub3 in subtrees:     #ewsi is always a matrix
            ewsi=np.dot(ewsi,I*_elem_weight_sub3(sub3,method,c,D,A,I,u))
        ews=np.dot(ews,ewsi) # ews is still a matrix
    
    return np.dot(ews,u) #returns a column



def elementary_weight_str(tree,method='BVRK',style='python'):
    """
        Constructs elementary weights for a Volterra Runge-Kutta method
        as strings, supposing the row sum condition.
        The output needs to be multiplied by b^T and equated to LHS_str
        in order to obtain the order condition.
        
        INPUT:
            
        - tree   -- input tree, must be a RootedTree.
        - method -- select which type of method to use: VRK, BVRK or if a PVRK
                    method is wanted, the three must be crated with the xa='a' flag.
        - style  -- output style: 'octave' or 'matlab', 'python' and 'latex'.
        
        OUTPUT: a string, mathematically is a column vector.

        **Examples**:

            >>> from nviepy import vrt
            >>> tree = vrt.list_trees(5)[6]
            >>> vrt.elementary_weight_str(tree)
            'dot(A,dot(A,c*d))'
            >>> vrt.elementary_weight_str(tree,style='matlab')
            '(A*(A*c.*d))'
            >>> vrt.elementary_weight_str(tree,style='latex')
            'AACD\\mathbf{u}'
            >>> vrt.elementary_weight_str(vrt.RootedTree('{a^10}'))
            'c**10*ones((s, 1),dtype=int)'
            >>> vrt.elementary_weight_str(vrt.RootedTree('{{a^11}x}'))
            'd*dot(A,c**11)'
    """
    
    if tree=='': return ''
    
    
    if tree=='a': 
        if style=='python':
            return 'ones((s, 1),dtype=int)'
        elif style=='latex':
            return r'\mathbf{u}'
        else: 
            return 'ones(s, 1)'
        
    if style=='python':
        c=C='*c'
        A="A,"
        opn="*dot("
        cls=")"
        if method=='VRK':
            d='*e'
            D='D'
        else :
            D=d='*d'
    elif style=='latex':
        c=C='C'
        A='A'
        opn="("
        cls=")"
        if method=='VRK':
            d='E'
            D=r'\!\circ\!D'
        else :
            d=D='D'
    else: 
        c=C='.*c'
        A='A*'
        opn=".*("
        cls=")"
        if method=='VRK':
            d='.*e'
            D="D"
        else :
            D=d='.*d'
            
    ewstr=''
    nx,na,subtrees=tree._parse_subtrees()
    lenst=len(subtrees)
    
    if na>0:
        ewstr+=c
        if na>1: ewstr+='^'+str(na)
    if nx>0:
        ewstr+=d
        if nx>1: ewstr+='^'+str(nx)
    
    
    # Two curly bracket contain at least a symbol, 'a' or 'x', i.e {} is not a leaf
    if lenst>0:
        for subtree in subtrees:
            if (lenst>1) or (style!='latex'): ewstr+=opn
            ewstr+=A+_elem_weight_sub3_str(subtree,method,style,C,D,A,opn,cls)
            if (lenst>1) or (style!='latex'): ewstr+=cls
    
    
    
    
    if style=='python':
        ewstr=ewstr.replace('*', '', 1) #must be first command
        ewstr=ewstr.replace('A,D', 'A*D')
        ewstr=ewstr.replace(',)', ',ones((s, 1),dtype=int))')
        ewstr=ewstr.replace(',*', ',')
        ewstr=ewstr.replace('^', '**')
        if ewstr[-1]!=')': ewstr+='*ones((s, 1),dtype=int)'  #no subtrees
    elif style=='latex':
        ewstr=ewstr.replace(')(', r')\!\circ\!(')
        ewstr=ewstr.replace(')', r'\mathbf{u})')
        if ewstr[-1]!=')': ewstr+=r'\mathbf{u}'
    else:
        ewstr=ewstr.replace('.*', '', 1) #must be first command
        ewstr=ewstr.replace('A*D', 'A.*D')
        ewstr=ewstr.replace(',)', '*ones(s, 1))')
        ewstr=ewstr.replace('*.*', '*')
        ewstr=ewstr.replace('^', '.^')
        if ewstr[-1]!=')': ewstr+='.*ones(s, 1)'  #no subtrees
    
    
    return ewstr
        


def _elem_weight_sub3_str(tree,method,style,C,D,A,opn,cls):
    
    if tree!='':
        
        
        ewsstr=''
        nx,na,subtrees=tree._parse_subtrees()
        lens3s=len(subtrees)
        
        if method=='VRK':
            if nx>0: 
                ewsstr+=D
                if nx>1: ewsstr+='^'+str(nx)
                if style!='latex': ewsstr+=','
            if na>0: 
                ewsstr+=C
                if na>1: ewsstr+='^'+str(na)
        else:
            if na>0:
                ewsstr+=C
                if na>1: ewsstr+='^'+str(na)
            if nx>0:
                ewsstr+=D
                if nx>1: ewsstr+='^'+str(nx)
            
        
        # Two curly bracket contain al least a symbol, 'a' or 'x', i.e {} is not a leaf
        if lens3s>0:
            for sub3 in subtrees:
                if (lens3s>1) or (style!='latex'): ewsstr+=opn
                ewsstr+=A+_elem_weight_sub3_str(sub3,method,style,C,D,A,opn,cls)
                if (lens3s>1) or (style!='latex'): ewsstr+=cls
    return ewsstr

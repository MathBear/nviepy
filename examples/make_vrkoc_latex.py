#Generates a latex file in order to tabulate order conditions for VRK and BVRK methods

from nviepy import volterra_rooted_3s as vrt


p_max=7  # choose maximum order to write

Order_colum=False # True or False; decide if you want the extra order column

LongTable=True # True or False,
               #if false, it will be a normal Table (useful for small p_max)
               #otherwise it will be a Long table spanning on multiple pages.



bexp="\\mathbf{b}^T "  # you can also choose the math render of the transpose
#bexp="\\mathbf{b}^\intercal "




# ------------   Code starts here -------------

f=open('oc_vrk_table.tex','w')

f.write("\documentclass{article}\n")
f.write("\usepackage[utf8]{inputenc} \n")
f.write("\usepackage{booktabs} \n")
f.write("\usepackage{tabularx} \n")
f.write("\usepackage{multirow} \n")
f.write("\usepackage{amsmath} \n")
if LongTable==True: f.write("\usepackage{longtable} \n")

f.write("\usepackage[a4paper,margin=0.8in]{geometry}  %  needed for big orders \n")

f.write("\\newcolumntype{L}{>{$}l<{$}} \n")
f.write("\\newcolumntype{C}{>{$}c<{$}} \n\n")

f.write("\\begin{document} \n\n")

f.write("\\renewcommand{\\arraystretch}{1.5} %  creates vertical space, and serves to show well frac{}{} \n\n")


if LongTable==False:
    f.write("\\begin{table}\n")
else:
    f.write("\\begin{center}\n")

f.write("\\footnotesize   %\\tiny  %\\scriptsize %\\footnotesize  %\\small %\\normalsize (default) %\\large \n")

if LongTable==True:
    f.write("\\begin{longtable}")
    if Order_colum==True:
        f.write("{CLr@{ }Lr@{ }L}\n")
    else:
        f.write("{Lr@{ }Lr@{ }L}\n")

f.write("\\caption{Order conditions from 1 to "+str(p_max)+"}\n")
f.write("\\label{tab:}")

if LongTable==False:
    f.write(" \n")
    f.write("\\centering\n")
    if Order_colum==True:
        f.write("\\begin{tabular}{CLr@{ }Lr@{ }L}\n")
    else:
        f.write("\\begin{tabular}{Lr@{ }Lr@{ }L}\n")
else:
    f.write(" \\\\ \n")

f.write("\\toprule\n")
if Order_colum==True: f.write("\\text{Order} & ")
f.write(" \\text{VRK} & & \\text{PVRK}")
f.write(" &  & \\text{BVRK} \\\\ \n")
if Order_colum==True: f.write("\\mathbf{p} &")
f.write("       &  \multicolumn{2}{l}{$(e_i=1,d_{ij}=c_i)$}")
f.write(" & \multicolumn{2}{l}{$(e_i=d_i, d_{ij}=d_j)$} \\\\ \n")
f.write("\\midrule\n")

if LongTable==True:
    f.write("\\endfirsthead\n")
    if Order_colum==True:
        f.write("\\multicolumn{6}{l}")
    else:
        f.write("\\multicolumn{5}{l}")
    f.write("{\\footnotesize\\itshape Continua dalla pagina precedente} \\\\ \n")
    #     {{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}}
    #     {\\tablename\\ \\thetable\\ -- \\textit{Continued from previous page}}
    f.write("\\toprule\n")
    if Order_colum==True: f.write("\\text{Order} & ")
    f.write(" \\text{VRK} & & \\text{PVRK}")
    f.write(" &  & \\text{BVRK} \\\\ \n")
    if Order_colum==True: f.write("\\mathbf{p} &  ")
    f.write("        &  \multicolumn{2}{l}{$(e_i=1,d_{ij}=c_i)$}")
    f.write(" & \multicolumn{2}{l}{$(e_i=d_i, d_{ij}=d_j)$} \\\\ \n")
    f.write("\\midrule\n")
    f.write("\\endhead\n")
    
    f.write("\\midrule\n")
    if Order_colum==True:
        f.write("\\multicolumn{6}{r}")
    else:
        f.write("\\multicolumn{5}{r}")
    f.write("{\\footnotesize\\itshape Continua nella prossima pagina} \\\\ \n")
    #        {\\textit{Continued on next page}}
    f.write("\\endfoot\n")
    f.write("\\bottomrule\n")
    f.write("\\endlastfoot\n")
    
ioc=1
iocP=1
p=0
many_Bforest=vrt.list_trees(p_max,'upto',xa='xa')
many_Pforest=vrt.list_trees(p_max,'upto',xa='a')
for Bforest in many_Bforest:
    Pforest=many_Pforest[p]
    iP=0
    p+=1
    if Order_colum==True:
        if LongTable==True:
            f.write("\\textbf{"+str(p)+"} ")
        else:
            f.write("\\multirow{"+str(len(Bforest))+"}*{\\textbf{"+str(p)+"}} ")
            # Alternative way for order column, when the table is on one page!
    for tree in Bforest:
        vrk_oc=vrt.elementary_weight_str(tree,method='VRK',style='latex')
        bvrk_oc=vrt.elementary_weight_str(tree,method='BVRK',style='latex')
        rhs =vrt.right_hand_side_str(tree,style='latex')
        if Order_colum==True: f.write("  &")
        f.write(" "+bexp+vrk_oc+"="+rhs)
        if tree==Pforest[iP]:
            if iP<len(Pforest)-1: iP+=1
            pvrk_oc=vrt.elementary_weight_str(tree,method='PVRK',style='latex')
            f.write("  & "+str(iocP)+". & "+bexp+pvrk_oc+"="+rhs)
            iocP+=1
        else:
            f.write("  &    &   ")
        f.write("  & "+str(ioc)+". & "+bexp+bvrk_oc+"="+rhs+" \\\\ \n")
        ioc+=1
    if p_max!=p: f.write("\\midrule % ---------- Here starts order "+str(p+1)+" ----------\n")

if LongTable==False:
    f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
else:
    f.write("\\end{longtable}\n")
    f.write("\\end{center}\n\n")

f.write("\\end{document} \n")
f.close()

PDBFILE=$1  
                   # bfactor |     
                   # cutoff  V
awk '{if ($1=="ATOM" && $11>.75) print $0}' ${PDBFILE} > B75_${PDBFILE}
ipython structure.py ${PDBFILE} B75_${PDBFILE} 
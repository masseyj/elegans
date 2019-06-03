from pyfasta import Fasta
import numpy as np

def complement(s): 
    newArray = np.copy(s)
    """Return the complementary sequence string.""" 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    for k, v in basecomplement.iteritems(): newArray[s==k] = v
    return newArray
#    letters = list(s) 
#    letters = [basecomplement[base] for base in letters] 
#    return ''.join(letters) 




#def reversecomplement(s):
    #"""Return the reverse complement of the dna string."""
    #s = s.reverse()
    #s = complement(s)
    #return s


    
#Replaced full chromosome names with just first word in command line
#sed 's/ .*//' GCF_000224195.1_Dele_2.0_genomic.fna > GCF_000224195.1_Dele_2.0_genomic.desc_stripped.fa
#sed 's/ .*//' GCF_000224195.1_Dele_2.0_genomic.gun_update.fa.updated.fasta  > GCF_000224195.1_Dele_2.0_genomic.gun_updated.desc_stripped.fa

#Toggle between first and second lines respectively to assemble the original and updated genome


#####
# Enter original genome name here
#####
#f = Fasta('GCF_000224195.1_Dele_2.0_genomic.ele_updated.desc_stripped.fa')
f = Fasta('GCF_000224195.1_Dele_2.0_genomic.gun_updated.desc_stripped.fa')

#####
# Enter assembled genome name here
#####
#faout = open('GCF_000224195.1_Dele_2.0_genomic.desc_stripped.ele_updated.chr_assembly_DLS_28x16.fa','w')
faout = open('GCF_000224195.1_Dele_2.0_genomic.gun_updated.desc_stripped.fa.chr_assembly_DLS_28x16.fa','w')

"""
relevant elegans scaffols


"""

"""
Chromosome arms - Not Muller chromosome arms!!!
X, B, C, D, E, F
"""


spacer=np.empty(1000,dtype='string')
spacer[:] = 'N'

#####
# List each contig that you want assembled for a particular chromosome here. 
# If the contig should be rev comped, then you need to add the two extra lines (like for X5)
# Then, in the np.concatenate line list the fragments (X1, X2, etc.) in the order they should be concatenated.
# If you want a spacer between contigs, then include "spacer"
#####


#Assemble chr X
print "Assembling chrA"
X1 = np.array(f['NW_016083076.1'])
X2 = np.array(f['NW_016083504.1'])
X3 = np.array(f['NW_016083526.1'])
X4 = np.array(f['NW_016083521.1'])
X5 = np.array(f['NW_016083551.1'])
revX5=X5[::-1]
rev_comp_X5=complement(revX5)
X6 = np.array(f['NW_016083509.1'])
revX6=X6[::-1]
rev_comp_X6=complement(revX6)
#X7 skipped intentionally
X8 = np.array(f['NW_016083864.1'])
X9 = np.array(f['NW_016082652.1'])
X10 = np.array(f['NW_016083701.1'])
lengthX10=len(X10)
X10_left = X10[0:699504]
X11 = np.array(f['NW_016083514.1'])
X12 = np.array(f['NW_016083021.1'])
X13 = np.array(f['NW_016083519.1'])
revX13=X13[::-1]
rev_comp_X13=complement(revX13)
X14 = np.array(f['NW_016083595.1'])
revX14=X14[::-1]
rev_comp_X14=complement(revX14)
X15 = np.array(f['NW_016083056.1'])
revX15=X15[::-1]
rev_comp_X15=complement(revX15)
X16 = np.array(f['NW_016083531.1'])
revX16=X16[::-1]
rev_comp_X16=complement(revX16)
X17 = np.array(f['NW_016082725.1'])
revX17=X17[::-1]
rev_comp_X17=complement(revX17)
X18 = np.array(f['NW_016083552.1'])
X19 = np.array(f['NW_016083780.1'])
lengthX19=len(X19)
X19_right = X19[100996:lengthX19]
X20 = np.array(f['NW_016083490.1'])
revX20=X20[::-1]
rev_comp_X20=complement(revX20)
X21 = np.array(f['NW_016083261.1'])
X22 = np.array(f['NW_016083690.1'])
X23 = np.array(f['NW_016083785.1'])
X24 = np.array(f['NW_016083750.1'])
revX24=X24[::-1]
rev_comp_X24=complement(revX24)
X25 = np.array(f['NW_016083621.1'])
revX25=X25[::-1]
rev_comp_X25=complement(revX25)

A = np.concatenate((X1,spacer,X2,spacer,X3,spacer,X4,spacer,rev_comp_X5,spacer,rev_comp_X6,spacer,X8,spacer,X9,spacer,X10_left,spacer,X11,spacer,X12,spacer,rev_comp_X13,spacer,rev_comp_X14,spacer,rev_comp_X15,spacer,rev_comp_X16,spacer,rev_comp_X17,spacer,X18,spacer,X19_right,spacer,rev_comp_X20,spacer,X21,spacer,X22,spacer,X23,spacer,rev_comp_X24,spacer,rev_comp_X25))

#Assemble chr 2
print "Assembling chrD"
B1 = np.array(f['NW_016083757.1'])
B2 = np.array(f['NW_016083763.1'])
B3 = np.array(f['NW_016083501.1'])
B4 = np.array(f['NW_016083074.1'])
revB4=B4[::-1]
rev_comp_B4=complement(revB4)
B5 = np.array(f['NW_016083776.1'])
B6 = np.array(f['NW_016083701.1'])
lengthB6=len(B6)
B6_right = B6[699505:lengthB6]
B7 = np.array(f['NW_016083508.1'])
B8 = np.array(f['NW_016083707.1'])
B9 = np.array(f['NW_016083541.1'])
revB9=B9[::-1]
rev_comp_B9=complement(revB9)
B10 = np.array(f['NW_016083775.1'])
B11 = np.array(f['NW_016083495.1'])
B12 = np.array(f['NW_016082605.1'])
revB12=B12[::-1]
rev_comp_B12=complement(revB12)
B13 = np.array(f['NW_016083533.1'])
B14 = np.array(f['NW_016082931.1'])

D = np.concatenate((B1, spacer, B2, spacer, B3, spacer, rev_comp_B4, spacer, B5, spacer, B6_right, spacer, B7, spacer, B8, spacer, rev_comp_B9, spacer, B10, spacer, B11, spacer, rev_comp_B12, spacer, B13, spacer, B14))

#Assemble chr C
print "Assembling chrC"
C1 = np.array(f['NW_016083615.1'])
C2 = np.array(f['NW_016083722.1'])
revC2=C2[::-1]
rev_comp_C2=complement(revC2)
C3 = np.array(f['NW_016083709.1'])
C4 = np.array(f['NW_016083723.1'])
revC4=C4[::-1]
rev_comp_C4=complement(revC4)
C5 = np.array(f['NW_016082211.1'])
revC5=C5[::-1]
rev_comp_C5=complement(revC5)
C6 = np.array(f['NW_016083535.1'])
revC6=C6[::-1]
rev_comp_C6=complement(revC6)
C7 = np.array(f['NW_016083857.1'])
C8 = np.array(f['NW_016083766.1'])
C9 = np.array(f['NW_016083740.1'])
C10 = np.array(f['NW_016083748.1'])
C11 = np.array(f['NW_016082993.1'])
revC11=C11[::-1]
rev_comp_C11=complement(revC11)
C12 = np.array(f['NW_016083507.1'])
C13 = np.array(f['NW_016083858.1'])
revC13=C13[::-1]
rev_comp_C13=complement(revC13)
C14 = np.array(f['NW_016083790.1'])
revC14=C14[::-1]
rev_comp_C14=complement(revC14)
C15 = np.array(f['NW_016083223.1'])
C16 = np.array(f['NW_016083530.1'])
C17 = np.array(f['NW_016083694.1'])
C17_left = C7[0:289487]
C18 = np.array(f['NW_016083791.1'])
C19 = np.array(f['NW_016083485.1'])
C20 = np.array(f['NW_016083719.1'])
C21 = np.array(f['NW_016083660.1'])
revC21=C21[::-1]
rev_comp_C21=complement(revC21)
C22 = np.array(f['NW_016083736.1'])
C23 = np.array(f['NW_016083725.1'])
C24 = np.array(f['NW_016083773.1'])
revC24=C24[::-1]
rev_comp_C24=complement(revC24)
C25 = np.array(f['NW_016083855.1'])
revC25=C25[::-1]
rev_comp_C25=complement(revC25)
C26 = np.array(f['NW_016083779.1'])
revC26=C26[::-1]
rev_comp_C26=complement(revC26)
C27 = np.array(f['NW_016083036.1'])
C28 = np.array(f['NW_016080073.1'])
C29 = np.array(f['NW_016083619.1'])
revC29=C29[::-1]
rev_comp_C29=complement(revC29)

C = np.concatenate((C1, spacer, rev_comp_C2, spacer, C3, spacer, rev_comp_C4, spacer,rev_comp_C5, spacer, rev_comp_C6, spacer, C7, spacer, C8, spacer, C9, spacer, C10, spacer, rev_comp_C11, spacer, C12, spacer, rev_comp_C13, spacer, rev_comp_C14, spacer, C15, spacer, C16, spacer, C17_left, spacer, C18, spacer, C19, spacer, C20, spacer, rev_comp_C21, spacer, C22, spacer, C23, spacer, rev_comp_C24, spacer, rev_comp_C25, spacer, rev_comp_C26, spacer, C27, spacer, C28, spacer, rev_comp_C29))

#Assemble chr D
print "Assembling chrD"

D1 = np.array(f['NW_016083746.1'])
D2 = np.array(f['NW_016083275.1'])

F = np.concatenate((D1, spacer, D2))

#Assemble chr E
print "Assembling chrE"

E1 = np.array(f['NW_016083612.1'])
E2 = np.array(f['NW_016083529.1'])
E3 = np.array(f['NW_016079081.1'])
revE3=E3[::-1]
rev_comp_E3=complement(revE3)
E4 = np.array(f['NW_016083502.1'])
revE4=E4[::-1]
rev_comp_E4=complement(revE4)
E5 = np.array(f['NW_016083476.1'])
E6 = np.array(f['NW_016083555.1'])
E7 = np.array(f['NW_016083769.1'])
E8 = np.array(f['NW_016083516.1'])
revE8=E8[::-1]
rev_comp_E8=complement(revE8)
E9 = np.array(f['NW_016083503.1'])
revE9=E9[::-1]
rev_comp_E9=complement(revE9)
E10 = np.array(f['NW_016083588.1'])
E11 = np.array(f['NW_016083525.1'])
E12 = np.array(f['NW_016082963.1'])
E13 = np.array(f['NW_016083788.1'])
E14 = np.array(f['NW_016083720.1'])
revE14=E14[::-1]
rev_comp_E14=complement(revE14)
E15 = np.array(f['NW_016083702.1'])
E16 = np.array(f['NW_016083704.1'])
revE16=E16[::-1]
rev_comp_E16=complement(revE16)

#E = np.concatenate((E1, spacer, rev_comp_E2, spacer, E3, spacer, E4, spacer, E5))

#early runs suggest that E4 and E5 are switched.
E = np.concatenate((E1, spacer, E2, spacer, rev_comp_E3, spacer, rev_comp_E4, spacer, E5, spacer, E6, spacer, E7, spacer, rev_comp_E8, spacer, rev_comp_E9, spacer, E10, spacer, E11, spacer, E12, spacer, E13, spacer, rev_comp_E14, spacer, E15, spacer, rev_comp_E16))


#Assemble chr F
print "Assembling chrF"

F1 = np.array(f['NW_016083198.1'])
F2 = np.array(f['NW_016083532.1'])
F3 = np.array(f['NW_016082909.1'])
F4 = np.array(f['NW_016083150.1'])
F5 = np.array(f['NW_016083479.1'])
F6 = np.array(f['NW_016083185.1'])
F7 = np.array(f['NW_016083694.1'])
lengthF7=len(F7)
F7_right = F7[289488:lengthF7]
F8 = np.array(f['NW_016083154.1'])
revF8=F8[::-1]
rev_comp_F8=complement(revF8)
F9 = np.array(f['NW_016083480.1'])
F9_left = F9[0:336703]
F10 = np.array(f['NW_016083536.1'])
revF10=F10[::-1]
rev_comp_F10=complement(revF10)
F11 = np.array(f['NW_016083089.1'])
revF11=F11[::-1]
rev_comp_F11=complement(revF11)
F12 = np.array(f['NW_016083554.1'])
F13 = np.array(f['NW_016083846.1'])
revF13=F13[::-1]
rev_comp_F13=complement(revF13)
F14 = np.array(f['NW_016083781.1'])
F15 = np.array(f['NW_016083277.1'])
F16 = np.array(f['NW_016083238.1'])
F17 = np.array(f['NW_016082968.1'])
F18 = np.array(f['NW_016083534.1'])
F19 = np.array(f['NW_016082727.1'])
F20 = np.array(f['NW_016082964.1'])

B = np.concatenate((F1, spacer, F2, spacer, F3, spacer, F4, spacer, F5, spacer, F6, spacer, F7_right, spacer, rev_comp_F8, spacer, F9_left, spacer, rev_comp_F10, spacer, rev_comp_F11, spacer, F12, spacer, rev_comp_F13, spacer, F14, spacer, F15, spacer, F16, spacer, F17, spacer, F18, spacer, F19, spacer, F20))


#get all other chrom keys
chroms = f.keys()
#####
# remove scaffolds that were used in new assembly
#####

chroms.remove("NW_016083076.1")
chroms.remove("NW_016083504.1")
chroms.remove("NW_016083526.1")
chroms.remove("NW_016083521.1")
chroms.remove("NW_016083551.1")
chroms.remove("NW_016083509.1")
chroms.remove("NW_016082652.1")
chroms.remove("NW_016083864.1")
chroms.remove("NW_016083701.1")
chroms.remove("NW_016083514.1")
chroms.remove("NW_016083021.1")
chroms.remove("NW_016083519.1")
chroms.remove("NW_016083595.1")
chroms.remove("NW_016083056.1")
chroms.remove("NW_016083531.1")
chroms.remove("NW_016082725.1")
chroms.remove("NW_016083552.1")
chroms.remove("NW_016083780.1")
chroms.remove("NW_016083490.1")
chroms.remove("NW_016083261.1")
chroms.remove("NW_016083690.1")
chroms.remove("NW_016083785.1")
chroms.remove("NW_016083750.1")
chroms.remove("NW_016083621.1")
chroms.remove("NW_016083757.1")
chroms.remove("NW_016083763.1")
chroms.remove("NW_016083501.1")
chroms.remove("NW_016083074.1")
chroms.remove("NW_016083776.1")
chroms.remove("NW_016083508.1")
chroms.remove("NW_016083707.1")
chroms.remove("NW_016083541.1")
chroms.remove("NW_016083775.1")
chroms.remove("NW_016083495.1")
chroms.remove("NW_016082605.1")
chroms.remove("NW_016083533.1")
chroms.remove("NW_016082931.1")
chroms.remove("NW_016083615.1")
chroms.remove("NW_016083722.1")
chroms.remove("NW_016083709.1")
chroms.remove("NW_016083723.1")
chroms.remove("NW_016082211.1")
chroms.remove("NW_016083535.1")
chroms.remove("NW_016083857.1")
chroms.remove("NW_016083766.1")
chroms.remove("NW_016083740.1")
chroms.remove("NW_016083748.1")
chroms.remove("NW_016082993.1")
chroms.remove("NW_016083507.1")
chroms.remove("NW_016083858.1")
chroms.remove("NW_016083790.1")
chroms.remove("NW_016083223.1")
chroms.remove("NW_016083530.1")
chroms.remove("NW_016083694.1")
chroms.remove("NW_016083791.1")
chroms.remove("NW_016083485.1")
chroms.remove("NW_016083719.1")
chroms.remove("NW_016083660.1")
chroms.remove("NW_016083736.1")
chroms.remove("NW_016083725.1")
chroms.remove("NW_016083773.1")
chroms.remove("NW_016083855.1")
chroms.remove("NW_016083779.1")
chroms.remove("NW_016083036.1")
chroms.remove("NW_016080073.1")
chroms.remove("NW_016083619.1")
chroms.remove("NW_016083746.1")
chroms.remove("NW_016083275.1")
chroms.remove("NW_016083612.1")
chroms.remove("NW_016083529.1")
chroms.remove("NW_016079081.1")
chroms.remove("NW_016083502.1")
chroms.remove("NW_016083476.1")
chroms.remove("NW_016083555.1")
chroms.remove("NW_016083769.1")
chroms.remove("NW_016083516.1")
chroms.remove("NW_016083503.1")
chroms.remove("NW_016083588.1")
chroms.remove("NW_016083525.1")
chroms.remove("NW_016082963.1")
chroms.remove("NW_016083788.1")
chroms.remove("NW_016083720.1")
chroms.remove("NW_016083702.1")
chroms.remove("NW_016083704.1")
chroms.remove("NW_016083198.1")
chroms.remove("NW_016083532.1")
chroms.remove("NW_016082909.1")
chroms.remove("NW_016083150.1")
chroms.remove("NW_016083479.1")
chroms.remove("NW_016083185.1")
chroms.remove("NW_016083154.1")
chroms.remove("NW_016083480.1")
chroms.remove("NW_016083536.1")
chroms.remove("NW_016083089.1")
chroms.remove("NW_016083554.1")
chroms.remove("NW_016083846.1")
chroms.remove("NW_016083781.1")
chroms.remove("NW_016083277.1")
chroms.remove("NW_016083238.1")
chroms.remove("NW_016082968.1")
chroms.remove("NW_016083534.1")
chroms.remove("NW_016082727.1")
chroms.remove("NW_016082964.1")

###########
# write new fasta file
###########
print "Writing file"
faout.write(">A\n")
faout.write(A.tostring())
faout.write("\n")
faout.write(">B\n")
faout.write(B.tostring())
faout.write("\n")
faout.write(">C\n")
faout.write(C.tostring())
faout.write("\n")
faout.write(">D\n")
faout.write(D.tostring())
faout.write("\n")
faout.write(">E\n")
faout.write(E.tostring())
faout.write("\n")
faout.write(">F\n")
faout.write(F.tostring())
faout.write("\n")

for contig in chroms:
    faout.write(">" + contig + "\n")
    faout.write(f[contig][:])
    faout.write("\n")

faout.close()




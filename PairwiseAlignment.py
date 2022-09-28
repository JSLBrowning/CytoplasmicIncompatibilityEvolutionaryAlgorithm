
# Import libraries
from Bio import pairwise2
from Bio.Seq import Seq

# "HWVTLVI#########YY#DSL########I###L#####D#########QQ#DG###CG####EN"  # ref for this consensus of nuc2 site is Gillespie et al., Tangle Web
# "DL#LL#R##########PIIIELK#####################DLVL##########PIGLELK"  # ref for this consensus of nuc2 site is Gillespie et al., Tangle Web
# "KRAR"
# "R-------R-R-R
# Creating sample sequences
seq1 = Seq("KRAR")
seq2 = Seq("MHWVTLVILIKRARAGRYLLTYYSDSLYVPSGHKFILTSLLTFLTDSPVSGKISIQQLDGRDFCGGEMDEN")

# Finding similarities
print(pairwise2.align.globalxs(seq1, seq2, -1, -.1, penalize_end_gaps=False)) #score_only=True
# gap opening, gap extension, gap opening, gap extension

"""
MELSIRVQPHKVLLRNQRLIGWPCDTDMCAVRKLSKQVCTLRARPDASFTGTERDTLKYADLVLLAQILFEHIVTEDVPALSARKLRDSRRIYILSVVRSSLVSTWPRPSTNRLKTCRSRVKSRVRIYKGDRGDYNADESMRWSTEAGIRRNCPEEYPEQERCYVEPHPPSLTERPPPEGLVMGLDVYKCVDPTLVIIGEIVLREILCKGMFCCKHTVNLLQTYCTVRLG*
TGHRHLKTAQNLAPKQRAVYHEPHRAAIATALAIPFIAAFSRHRYPRSAVIGDIKLAPCTREPNKPSADDLATYRAGSRAPSDCGNGRTVSGKSSLSGSEMT*
DUB fitness: 0.6613043478260882
NUC fitness: 0.7491304347826098
NLS fitness: 0.965
TYPEIV fitness: 0.9450000000000001
Binding fitness: 0.8181818181818182
combined fitness: 4.138616600790516
NLS Site: 1
TYPEIV Site: 1"""
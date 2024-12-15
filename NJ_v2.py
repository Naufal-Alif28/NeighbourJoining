'''
Simple Model Program for Constructing DNA Sequence-base Phylogenetic Tree using Neighbour-Joining Algorithm
by Naufal Muhammad Alif - 10422039
supplementary of Tugas Makalah Matematika Diskrit IF1202 2024-2025

*not equipped with input validation nor error-handling capability
*only accept multiple sequence alignment (MSA) file in FASTA format
*FASTA entry identifier must follow the standard NCBI nomenclature
'''

# IMPORT
import math

# OBJECT DEFINITION
# Defining object <Tree>
class MainTree:
    def __init__(self, taxa=[], children=[], height=0, n_intnodes=0):
        self.taxa = taxa
        self.children = children
        self.height = height
        self.n_intnodes = n_intnodes

# Defining object <Taxon>
class Taxon:
    def __init__(self, seq="", id="", branch_len=0):
        self.seq = seq
        self.id = id
        self.branch_len = branch_len

# Defining object <IntNodes>
class IntNodes:
    def __init__(self, children=[], child_taxon=[], branch_len=0):
        self.children = children
        self.child_taxon = child_taxon
        self.branch_len = branch_len

# FUNCTION DEFINITION FOR NEIGHBOUR-JOINING
# Transition & transversion ratio calculator for Kimura Two-Parameter model <TransRatio>
def TransRatio(taxon_i, taxon_j, len_msa):
    # Define transitions (p) and transversions (q)
    transitions = {
        "A": "G", 
        "G": "A", 
        "T": "C", 
        "C": "T"
    }
    transversions = {
        "A": {"T", "C"},
        "G": {"T", "C"},
        "T": {"A", "G"},
        "C": {"A", "G"}
    }
    p, q = 0, 0
    seq_i = taxon_i.seq
    seq_j = taxon_j.seq

    for idx in range(len_msa):
        base_i = seq_i[idx]
        base_j = seq_j[idx]
        if base_i != "-" and base_j != "-":  # Valid bases
            if transitions.get(base_i) == base_j:
                p += 1  # count transitions
            elif base_j in transversions.get(base_i, set()):
                q += 1  # count transversions

    p /= len_msa
    q /= len_msa

    return {"p": p, "q": q}

# Hamming distance calculator for Jukes-Cantor model
def HammingDist(taxon_i, taxon_j, len_msa):
    p = 0
    seq_i = taxon_i.seq
    seq_j = taxon_j.seq
    for idx in range(len_msa):
        base_i = seq_i[idx]
        base_j = seq_j[idx]
        if base_i != base_j:
            p += 1 # count difference
        
    p /= len_msa

    return p

# Distance matrix calculator (D(i,j)) <DistMatrix>
def DistMatrix(tree, len_msa, mode):
    degree = len(tree.children)
    dist_matrix = [[0 for i in range(degree)] for j in range(degree)]
    for i in range(degree):
        for j in range(degree):
            if i > j:   # filling bottom diagonal
                # flatten & dividing child taxa of the two nodes
                if isinstance(tree.children[i],Taxon):
                    bipartite_1 = [tree.children[i]]
                elif isinstance(tree.children[i],IntNodes):
                    bipartite_1 = [taxon for taxon in tree.children[i].child_taxon]
                if isinstance(tree.children[j],Taxon):
                    bipartite_2 = [tree.children[j]]
                elif isinstance(tree.children[j],IntNodes):
                    bipartite_2 = [taxon for taxon in tree.children[j].child_taxon]
                # calculate sum of distance between child taxa of two nodes
                dist = 0
                for k in bipartite_1:
                    for l in bipartite_2:
                        if mode == "K": # when using Kimura two-parameter model
                            p = TransRatio(k,l,len_msa)["p"]
                            q = TransRatio(k,l,len_msa)["q"]
                            dist +=  -0.5*math.log(1-2*p-q)-0.25*math.log(1-2*q)
                        elif mode == "J": # when using Jukes-Cantor model
                            p = HammingDist(k,l,len_msa)
                            dist += -0.75*math.log(1-(4/3)*p)
                dist /= len(bipartite_1)*len(bipartite_2)
                dist_matrix[i][j] = dist    # value assignment into distance matrix
    
    return dist_matrix

# Total divergence calculator (R(i)) <QMatrix>
def QMatrix(dist_matrix):
    q_matrix = []
    degree = len(dist_matrix)
    for k in range(degree):
        r_sum = 0
        for i in range(degree):
            for j in range(degree):
                if (i == k) ^ (j == k):
                    r_sum += dist_matrix[i][j]
        q_matrix.append(r_sum)     
    
    return q_matrix

# Adjusted pairwise distance calculator (M(i,j)) <AdjDist>
def AdjDist(dist_matrix, q_matrix, tree):
    D = dist_matrix
    R = q_matrix
    N = len(tree.taxa)
    degree = len(tree.children)
    dist_matrix_adj = [[0 for i in range(degree)] for j in range(degree)]
    for i in range(degree):
        for j in range(degree):
            dist_matrix_adj[i][j] = D[i][j]*(N-2)-((R[i]+R[j]))
    
    return dist_matrix_adj

# Branch length calculator <BranchLen>
def BranchLen(tree, imin, jmin, dist_matrix, q_matrix):
    N = len(tree.taxa)
    brlen_i = 0.5*(dist_matrix[imin][jmin]) + (1/(2*(N-2))*(q_matrix[imin] - q_matrix[jmin]))
    brlen_j = dist_matrix[imin][jmin] - brlen_i

    return {"i": brlen_i, "j": brlen_j}

# New node constructor & tree modificator <ModifTree>
def ModifTree(tree, dist_matrix_adj):
    # finding pair with minimum adjusted distance (M)
    degree = len(dist_matrix_adj)
    min_M = 1
    imin = -1
    jmin = -1
    for i in range(degree):
        for j in range(degree):
            if i > j:
                if min_M > dist_matrix_adj[i][j]:
                    min_M = dist_matrix_adj[i][j]
                    imin, jmin = i, j
    new_children = [tree.children[imin],tree.children[jmin]]
    new_taxon = []
    # compute branch length
    brlen = BranchLen(
        tree,imin,jmin,dist_matrix,q_matrix
    )
    tree.children[imin].branch_len = brlen["i"]
    tree.children[jmin].branch_len = brlen["j"]
    # create a new node
    for k in (imin,jmin):
        if isinstance(tree.children[k],IntNodes):
            for i in range(len((tree.children[k].child_taxon))):
                new_taxon.append(tree.children[k].child_taxon[i])
        elif isinstance(tree.children[k],Taxon):
            new_taxon.append(tree.children[k])
    globals()[f"IntNode_{tree.n_intnodes}"] = IntNodes(new_children,new_taxon)
    # replace merged OTU in tree with a new internal node
    tree.children.remove(tree.children[imin])
    tree.children.remove(tree.children[jmin])
    tree.children.append(globals()[f"IntNode_{tree.n_intnodes}"])
    tree.height += 1
    tree.n_intnodes += 1

# NEIGHBOUR-JOINING ALGORITHM
# MSA FASTA input 
print("\n" + "-"*20 + " NEIGHBOUR-JOINING " + "-"*20 + "\n")
dir_msa = input(">> INPUT\ninsert multiple sequence alignment FASTA file (directory)\n")
# choose distance model
print("\n" + "-"*50)
mod = input("\n>> MODEL\nselect distance model\nKimura[K]  Jukes-Cantor[J]\n")
print("\n" + "-"*20 + "end" + "-"*20 + "\n")

# reading MSA FASTA file
msa_fasta = open(dir_msa,"r")
msa_readlines = [line.replace("\n","") for line in msa_fasta.readlines()] #with \n removal
n_taxa = int(len(msa_readlines)/2)
len_msa = len(msa_readlines[1])

# Taxon assignment
taxa = []
idx = 0
for i in range(n_taxa):
    id = msa_readlines[idx].split(" ")
    id = f"{id[0]}{id[1]}{id[2][0]}.{id[3]}"
    globals()[f"taxon_{i}"] = Taxon(
        id = id,
        seq = msa_readlines[idx+1],
    )
    idx += 2
    taxa.append(globals()[f"taxon_{i}"])

# Tree assignment
main_tree = MainTree(
    taxa=taxa,
    children=taxa
)

# Neighbour-Joining
while len(main_tree.children) > 2:
    # distance matrix construction
    dist_matrix = DistMatrix(tree=main_tree, len_msa=len_msa, mode=mod)
    # Q-matrix construction
    q_matrix = QMatrix(dist_matrix=dist_matrix)
    # adjusted distance matrix construction
    dist_matrix_adj = AdjDist(dist_matrix=dist_matrix, q_matrix=q_matrix, tree=main_tree)
    # neighbour joining, creating a new node
    ModifTree(tree=main_tree, dist_matrix_adj=dist_matrix_adj)

# CONVERSION INTO NEWICK TREE
# converter from object <Tree> into a newick tree string
def MakeNewick(node):
    if isinstance(node, Taxon):     # leaf node
        return f"{node.id}:{node.branch_len:.6f}"
    elif isinstance(node, IntNodes):    # internal node with children
        child_newicks = [MakeNewick(child) for child in node.children]
        return f"({','.join(child_newicks)}):{node.branch_len:.6f}"
    elif isinstance(node, MainTree):    # main tree with children
        if len(node.children) == 2:  # root node with two children
            child_newicks = [MakeNewick(child) for child in node.children]
            return f"({','.join(child_newicks)});"

# newick tree construction
newick_tree = MakeNewick(main_tree)

# TREE VISUALISATION
# import
from Bio import Phylo # Phylo package from Biopython
from io import StringIO
import matplotlib.pyplot as plt # pyplot for visualization
# output
plt.rc('font', size=6) # font size adjustment
tree = Phylo.read(StringIO(newick_tree), "newick") 
Phylo.draw(tree, branch_labels=lambda c: c.branch_length) # draw tree with branch length (weight) label
plt.show()
import sys
import argparse
from argparse import RawTextHelpFormatter
from ete3 import Tree, NodeStyle, TextFace

parser = argparse.ArgumentParser(
    description = '''exnn_ete3.py: EXtract Neighbor Names with ETE3
    The ETE3-based script helping to extract names of neighbors
    from a phylogenetic tree.  The script draws the phylogenetic
    tree with highlighted pivotal points (target leaves and MRCA
    of targets). The user can inspect the tree and remove
    non-neighbor leaves interactively. Use the 'Cut partition'
    command from the context menu to remove clades.  After the
    ETE Tree Browser window is closed, the tree is saved and
    labels of the remained leaves are written.''',
    
    epilog = '''Author: Ivan Tsers
E-mail: tsers@evolbio.mpg.de
ETE reference: Jaime Huerta-Cepas, François Serra and Peer
               Bork. "ETE 3: Reconstruction, analysis and
               visualization of phylogenomic data." Mol Biol
               Evol (2016) doi: 10.1093/molbev/msw046
ETE GitHub repo: https://github.com/etetoolkit/ete''',
                formatter_class = RawTextHelpFormatter)
                
parser.add_argument('tree', type = str, default = "tree.nwk", help = 'Tree in the Newick format')
parser.add_argument('targets', type = str, default = "targets.txt", help = 'Newline-separated target names')
parser.add_argument('neighbors', type = str, default = "neighbors.txt", help = 'Newline-separated neighbor names (the output file)')

args = parser.parse_args()

# read the tree
t = Tree(args.tree, format = 1)
t_nleaves = len(t)

# read the targets
targ = []
targ_f = open(args.targets)
targ = targ_f.read().splitlines()
targ_f.close()
print("\nThere are", len(targ), "targets (marked in red).")
print(" The MRCA of all targets is shown with red circle.")
print(" The MRCA of the majority of targets is shown with blue circle.")
print(" Cut the clades you wish to exclude from the neighbor set.")

# find the MRCA of all targets
mrca_all = t.get_common_ancestor(targ)

# Find the MRCA of the majority of targets
cRoot = t.get_tree_root()
losVotes = 0
losLeaves = None

while losVotes != losLeaves:
  rChild = cRoot.get_children()[0]
  lChild = cRoot.get_children()[1]
  
  # Collect leaf labels
  rLabels = []
  for leaf in rChild:
    rLabels.append(leaf.name)
    
  lLabels = []
  for leaf in lChild:
    lLabels.append(leaf.name)

  # count votes
  rVotes = sum(label in targ for label in rLabels)
  lVotes = sum(label in targ for label in lLabels)

  # choose winner
  if rVotes > lVotes:
    cRoot = rChild
    losVotes = lVotes
    losLeaves = len(lLabels)
  else:
    cRoot = lChild
    losVotes = rVotes
    losLeaves = len(rLabels)

mrca_majority = cRoot.up

### Drawing ###

# initialize styles
#    NodeStyles
target_style = NodeStyle()
target_style["bgcolor"] = "#FF5555"

mrca_all_style = NodeStyle()
mrca_all_style["fgcolor"] = "red"
mrca_all_style["size"] = 15

mrca_majority_style = NodeStyle()
mrca_majority_style["fgcolor"] = "#2A7FFF"
mrca_majority_style["size"] = 13

#    TextFaces
mrca_all_label = TextFace("MRCA of all targets", fgcolor = "#FF5555")
mrca_majority_label = TextFace("MRCA of the majority of targets", fgcolor = "#2A7FFF")

# mark the target leaves
for node in t.traverse():
    if node.name in targ:
        node.set_style(target_style)
        
# mark the MRCA
mrca_all.set_style(mrca_all_style)
mrca_all.add_face(mrca_all_label, column = 1, position = "branch-bottom")

mrca_majority.set_style(mrca_majority_style)
mrca_majority.add_face(mrca_majority_label, column = 1, position = "branch-bottom")

# now show the tree
t.show()

# the cutting summary
cut_nleaves = t_nleaves - len(t)
print("You have removed", cut_nleaves, "of", t_nleaves, "leaves.")

# Return the labels of
print("\nLabels of the extracted leaves (neighbors)\n")

out_leaves = []
for leaf in t.iter_leaves():
    out_leaves.append(leaf.name)

print(*out_leaves, sep = ",")
print("\nThe list is written to", args.neighbors)

with open(args.neighbors, "w") as o:
    o.write("\n".join(out_leaves))
    o.close()

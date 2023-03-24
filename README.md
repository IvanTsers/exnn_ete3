# exnn_ete3.py: EXtract Neighbor Names with ETE3

The ETE3-based script helping to extract names of neighbors
from a phylogenetic tree.  The script draws the phylogenetic
tree with highlighted pivotal points (target leaves and MRCA
of targets). The user can inspect the tree and remove
non-neighbor leaves interactively. Use the 'Cut partition'
command from the context menu to remove clades.  After the
ETE Tree Browser window is closed, the tree is saved and
labels of the remained leaves are written.

## Install dependencies
`pip install ete3`

## Usage
`exnn_ete3.py [-h] tree targets neighbors`

**positional arguments**:
  `tree`        Tree in the Newick format
  `targets`     Newline-separated target names
  `neighbors`   Newline-separated neighbor names (the output file)

## References
Jaime Huerta-Cepas, François Serra and Peer Bork. 
"ETE 3: Reconstruction, analysis and visualization of phylogenomic data." 
*Mol Biol Evol* **(2016)** doi: 10.1093/molbev/msw046


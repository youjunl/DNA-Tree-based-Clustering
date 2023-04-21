import sys
sys.path.append('build')
import tree

tr = tree.new_tree(8)
print(tr)
inputs = []
inputs.append('ATTGCATA')
inputs.append('ATTGCATA')
inputs.append('ATTGCACA')
inputs.append('ATTGTACC')
inputs.append('CCTGAACT')
inputs.append('CCTGAACA')

clust_ind = 1
for inp in inputs:
    out = tree.search(tr, inp, 2)
    print('%s %d %d'%(inp, out.label, out.distance))
    if out.label < 0:
        tree.insert(tr, inp, clust_ind)
        clust_ind += 1
    
# print(out)
# print(out.label)
# print(out.distance)

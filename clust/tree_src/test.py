import sys
sys.path.append('build')
import tree
import random
from tqdm import tqdm

def editDistance(str1, str2, m, n):
 
    # If first string is empty, the only option is to
    # insert all characters of second string into first
    if m == 0:
        return n
 
    # If second string is empty, the only option is to
    # remove all characters of first string
    if n == 0:
        return m
 
    # If last characters of two strings are same, nothing
    # much to do. Ignore last characters and get count for
    # remaining strings.
    if str1[m-1] == str2[n-1]:
        return editDistance(str1, str2, m-1, n-1)
 
    # If last characters are not same, consider all three
    # operations on last character of first string, recursively
    # compute minimum cost for all three operations and take
    # minimum of three values.
    return 1 + min(editDistance(str1, str2, m, n-1),    # Insert
                   editDistance(str1, str2, m-1, n),    # Remove
                   editDistance(str1, str2, m-1, n-1)    # Replace
                   )

if __name__=='__main__':

    tr = tree.new_tree(8)
    inputs = []
    inputs.append('ATTGCATA')
    inputs.append('ATTGCATT')
    clust_ind = 1
    for inp in inputs:
        out = tree.search(tr, inp, 4)
        print('%s %d %d'%(inp, out.label, out.distance))
        if out.label < 0:
            tree.insert(tr, inp, clust_ind)
            clust_ind += 1
        
    n = 8
    sim = 1000
    tr = tree.new_tree(n)
    DNAbet = 'ATCG'
    ori_seq = ''.join(DNAbet[random.randint(0, 3)] for _ in range(n))
    tree.insert(tr, ori_seq, 1)
    cnt = 0
    tau = 6
    for _ in range(sim):
        seq = ''.join(DNAbet[random.randint(0, 3)] for _ in range(n))
        # Compute tree output
        result = tree.search(tr, seq, tau)
        result = min(result.distance, tau)

        # https://www.geeksforgeeks.org/edit-distance-dp-5/
        py_result = editDistance(ori_seq, seq, n, n)
        py_result = min(py_result, tau)

        # Compare
        if result == py_result:
            cnt += 1
    print('Validation: %d/%d success.'%(cnt, sim))
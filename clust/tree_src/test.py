import sys
sys.path.append('build')
import tree
import random
from tqdm import tqdm

def editDistance(A, B):
    N, M = len(A), len(B)
    # Create an array of size NxM
    dp = [[0 for i in range(M + 1)] for j in range(N + 1)]

    # Base Case: When N = 0
    for j in range(M + 1):
        dp[0][j] = j
    # Base Case: When M = 0
    for i in range(N + 1):
        dp[i][0] = i
    # Transitions
    for i in range(1, N + 1):
        for j in range(1, M + 1):
            if A[i - 1] == B[j - 1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j], # Insertion
                    dp[i][j-1], # Deletion
                    dp[i-1][j-1] # Replacement
                )

    return dp[N][M]


if __name__=='__main__':

    tr = tree.new_tree(8)
    inputs = []
    inputs.append('ATTGCATA')
    inputs.append('ATTGCATT') # 1
    inputs.append('ATTGCATA') # 0
    inputs.append('ATCGCATA') # 1
    inputs.append('ATGCATAT') # 2
    # inputs.append('ATTGCGAT') # 2
    
    clust_ind = 1
    for inp in inputs:
        out = tree.search(tr, inp, 6)
        out_quick = tree.quick_search(tr, inp, 6, 2)
        # print('%s %d %d'%(inp, out.label, out.distance))
        print('%s %d %d'%(inp, out_quick.label, out_quick.distance))
        if out.label == 0:
            tree.insert(tr, inp, clust_ind)
            clust_ind += 1
        
    n = 14
    sim = 10000
    tr = tree.new_tree(n)
    DNAbet = 'ATCG'
    ori_seq = ''.join(DNAbet[random.randint(0, 3)] for _ in range(n))
    tree.insert(tr, ori_seq, 1)
    cnt = 0
    tau = 6
    for i in tqdm(range(sim)):
        seq = ''.join(DNAbet[random.randint(0, 3)] for _ in range(n))
        # Compute tree output
        result = tree.search(tr, seq, tau)
        result_quick = tree.quick_search(tr, seq, tau, 3)
        result = min(result.distance, tau)
        # https://www.geeksforgeeks.org/edit-distance-dp-5/
        py_result = min(editDistance(ori_seq, seq), tau)

        # Compare
        if result == py_result:
            cnt += 1
        
    print('Validation: %d/%d success.'%(cnt, sim))
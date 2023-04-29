import sys
from clust import tree
import random
trie = tree.new_tree(21)
tree.insert(trie, 'GCCAGATACAAGCGCATACTG', 1)
tree.insert(trie, 'GCCAGATACAAAGCGCACGTG', 2)
#                          *AAGCGCATACT_
result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
# result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
# result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
# result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
# result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
print(result.label, result.distance)
DNAbet = 'ATCG'
for i in range(1000):
    ori_seq = ''.join(DNAbet[random.randint(0, 3)] for _ in range(21))
    tree.insert(trie, ori_seq, 1)
    result = tree.search(trie, 'GCCAGATACAAAGCGCATACT', 3)
    # result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
    # result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
    # result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
    # result = tree.quick_search(trie, 'GCCAGATACAAAGCGCATACT', 3, 3)
    print(result.label, result.distance)
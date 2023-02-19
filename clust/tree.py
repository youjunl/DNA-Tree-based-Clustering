"""Tree Structure Module

This module provides retrieval algorithms based on tree structures, 
including local greedy search algorithms.


"""
class Node:
    def __init__(self) -> None:
        self.node_nums = 4
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = -1

class Trie:
    def __init__(self):
        self.dna_to_num = {"A":0,"T":1,"G":2,"C":3}
        self.num_to_dna = {0:'A',1:'T',2:'G',3:'C'}
        self.node_nums = 4
        self.maxOptimDepth = 3
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = -1

    def insert(self, word: str,label:str) -> None:
        node = self
        for ch in word:
            ch = self.dna_to_num[ch]
            if not node.children[ch]:
                node.children[ch] = Node()
            node = node.children[ch]
        node.isEnd = label

    def delete(self,word:str):
        node = self
        for ch in word:
            ch = self.dna_to_num[ch]
            if not node.children[ch]:
                node.children[ch] = Node()
            node = node.children[ch]
        node.isEnd = -1

    def local_greedy_traversal(self, word, node, pos, error, maxErr):

        # Return True when traverse to the end of the tree
        if node.isEnd > 0:
            return node.isEnd
        
        # Return False when current error larger than thresold
        currentErr = sum(error)
        if currentErr > maxErr or pos >= len(word):
            return -1
        
        # Else continue traversal
        ch = word[pos]
        chNum = self.dna_to_num[ch]
        if node.children[chNum]:
            result = self.local_greedy_traversal(word, node.children[chNum], pos+1, error, maxErr)
            # Return if the traversal from this node is successful or fail
            return result

        # If there is not a node exist for word[pos]
        else:
            # Else get possible traversal nodes
            ins_list = []
            del_list = []
            sub_list = []
            tmp_list=[i for i in range(self.node_nums) if node.children[i]]
            maxOptimDepth = self.maxOptimDepth
            traverseNum = min(len(word)-pos-1, maxOptimDepth)         
            for chNum in tmp_list:
                #Deletion
                depth = 0
                tmp = node.children[chNum]
                while depth < traverseNum:
                    if not tmp.children[self.dna_to_num[word[pos+depth]]]:
                        break
                    tmp = tmp.children[self.dna_to_num[word[pos+depth]]]
                    depth += 1

                if depth == traverseNum:
                    del_list.append(chNum)

                #Insertion
                if len(ins_list) == 0:
                    depth = 0
                    tmp = node
                    while depth < traverseNum-1:
                        if not tmp.children[self.dna_to_num[word[pos+depth+1]]]:
                            break
                        tmp = tmp.children[self.dna_to_num[word[pos+depth+1]]]
                        depth += 1
                    
                    if depth == traverseNum-1:
                        ins_list.append(chNum)                     

                #Substitution
                depth = 0
                tmp = node.children[chNum]
                while depth < traverseNum:                                                
                    if not tmp.children[self.dna_to_num[word[pos+depth+1]]]:
                        break
                    tmp = tmp.children[self.dna_to_num[word[pos+depth+1]]]
                    depth += 1

                if depth == traverseNum:
                    sub_list.append(chNum)

                # Fix errors and continue traversal
                insErr, delErr, subErr = error
                # Insertion Fix
                if len(ins_list) > 0:
                    newWord = ''.join([word[:pos], word[pos-len(word)+1:], 'A'])
                    errorList = [insErr+1, delErr, subErr]
                    fixResult = self.local_greedy_traversal(newWord, node, pos+1, errorList, maxErr)

                    if fixResult > 0:
                        return fixResult

                # Deletion Fix
                for chNum in del_list:
                    newWord = ''.join([word[:pos], self.num_to_dna[chNum], word[pos-len(word):-1]])
                    errorList = [insErr, delErr+1, subErr]
                    fixResult = self.local_greedy_traversal(newWord, node.children[chNum], pos, errorList, maxErr)

                    if fixResult > 0:
                        return fixResult
                    
                # Substitution Fix
                for chNum in sub_list:
                    newWord = ''.join([word[:pos], self.num_to_dna[chNum], word[pos-len(word)+1:]])
                    errorList = [insErr, delErr, subErr+1]
                    fixResult = self.local_greedy_traversal(newWord, node.children[chNum], pos+1, errorList, maxErr)

                    if fixResult > 0:
                        return fixResult
        
            return -1
            
    def match(self, word, max_value, mode = 0):
        node = self
        result = self.local_greedy_traversal(word, node, 0, [0, 0, 0], max_value)
        if type(result) == list:
            return -1
        return result
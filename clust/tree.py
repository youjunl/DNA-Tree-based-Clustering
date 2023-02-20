class Node:
    def __init__(self) -> None:
        self.node_nums = 4
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = False   
        self.singleBranch = True
class Trie:
    def __init__(self):
        self.dna_dict = {"A":0,"T":1,"G":2,"C":3}
        self.node_nums = 4
        self.maxOptimDepth = 3
        self.minTerminate = 9
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = False
        self.singleBranch = False

    def insert(self, word: str,label:str) -> None:
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Node()
                node = node.children[ch]
            else:
                node = node.children[ch]
                node.singleBranch = False
        node.isEnd = label

    def delete(self,word:str):
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = False

    def fuzz_align_normal(self, dnaInfo):
        word, e = dnaInfo
        currentErr = sum(e)
        node = self 
        tmp_list=[]
        sub_list=[]
        ins_list=[]
        del_list=[]
        len_=len(word)

        for num, ch in enumerate(word) :
            chNum = self.dna_dict[ch]
            if False:#node.singleBranch and num > self.minTerminate:
                tmp_list=[i for i in range(self.node_nums) if node.children[i]]
                node = node.children[tmp_list[0]]

            elif not node.children[chNum]:
                # Collect existed node
                tmp_list=[i for i in range(self.node_nums) if node.children[i]]                
                maxOptimDepth = self.maxOptimDepth
                traverseNum = min(len_-num-1, maxOptimDepth)
                for k in tmp_list:
                    #Deletion
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:
                        if not tmp.children[self.dna_dict[word[num+depth]]]:
                            break
                        tmp = tmp.children[self.dna_dict[word[num+depth]]]
                        depth += 1

                    if depth == traverseNum:
                        del_list.append(k)

                    #Insertion
                    if len(ins_list) == 0:
                        depth = 0
                        tmp = node
                        while depth < traverseNum-1:
                            if not tmp.children[self.dna_dict[word[num+depth+1]]]:
                                break
                            tmp = tmp.children[self.dna_dict[word[num+depth+1]]]
                            depth += 1
                        
                        if depth == traverseNum-1:
                            ins_list.append(k)                     

                    #Substitution
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:                                                
                        if not tmp.children[self.dna_dict[word[num+depth+1]]]:
                            break
                        tmp = tmp.children[self.dna_dict[word[num+depth+1]]]
                        depth += 1

                    if depth == traverseNum:
                        sub_list.append(k)

                return [num, ins_list, del_list, sub_list]

            else:
                node = node.children[chNum]

        if type(node.isEnd) == int:
            return node.isEnd

        elif type(node.isEnd) == bool:
            tmp_list=[i for i in range(self.node_nums) if node.children[i]]
            return [num, [], tmp_list, []]

    def fuzz_align_sensitive(self, dnaInfo):
        word, e = dnaInfo
        currentErr = sum(e)
        node = self 
        tmp_list=[]
        sub_list=[]
        ins_list=[]
        del_list=[]
        len_=len(word)

        for num, ch in enumerate(word) :
            chNum = self.dna_dict[ch]
            if not node.children[chNum]:
                # Collect existed node
                tmp_list=[i for i in range(self.node_nums) if node.children[i]]                
                maxOptimDepth = self.maxOptimDepth
                traverseNum = min(len_-num-1, maxOptimDepth)
                for k in tmp_list:
                    #Deletion
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:
                        if not tmp.children[self.dna_dict[word[num+depth]]]:
                            break
                        tmp = tmp.children[self.dna_dict[word[num+depth]]]
                        depth += 1

                    if depth == traverseNum or num - 2*(currentErr - traverseNum + depth) > 0:
                        del_list.append(k)

                    #Insertion
                    if len(ins_list) == 0:
                        depth = 0
                        tmp = node
                        while depth < traverseNum-1:
                            if not tmp.children[self.dna_dict[word[num+depth+1]]]:
                                break
                            tmp = tmp.children[self.dna_dict[word[num+depth+1]]]
                            depth += 1
                        
                        if depth == traverseNum-1 or num - 2*(currentErr - traverseNum + depth) > 2:
                            ins_list.append(k)                     

                    #Substitution
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:                                                
                        if not tmp.children[self.dna_dict[word[num+depth+1]]]:
                            break
                        tmp = tmp.children[self.dna_dict[word[num+depth+1]]]
                        depth += 1

                    if depth == traverseNum or num - 2*(currentErr - traverseNum + depth) > 2:
                        sub_list.append(k)

                return [num, ins_list, del_list, sub_list]

            else:
                node = node.children[chNum]

        if type(node.isEnd) == int:
          return node.isEnd

        elif type(node.isEnd) == bool:
          tmp_list=[i for i in range(self.node_nums) if node.children[i]]
          return [num, [], tmp_list, []]

    def fuzz_fin(self, word, max_value, mode = 0):
        queue = []
        queue.append([word,[0, 0, 0]])        
        output=["",[100, 100, 100]]
        minOutputSum = 300
        num2dnaDict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
        adjacent = []
        adjacentList = []
        if mode == 0:
            alignFunc = self.fuzz_align_normal
        else:
            alignFunc = self.fuzz_align_sensitive
        while True :
            if sum(output[1]) == 0 or queue == []:
                break

            word, errorNums = queue.pop(0)
            errorSum = sum(errorNums)
            if errorSum > max_value or errorSum > minOutputSum:
              # Current number of error is larger than the threshold
                continue

            result = alignFunc([word, errorNums])

            if type(result) == int:
                adjacent.append([result, errorNums])
                adjacentList.append(result)
                errorSum = sum(errorNums)
                if errorSum < minOutputSum :
                    output=[result, errorNums]
                    minOutputSum = errorSum

            elif result[0] == len(word)-1:
                for i in range(len(result[1])):
                    chNum = result[1][i]
                    k = word[:result[0]]+num2dnaDict[chNum]
                    errorNums[1] += 1
                    queue.append([k,errorNums])

            else:
                pos, ins_list, del_list, sub_list = result
                insErr, delErr, subErr = errorNums

                # Insertion Fix
                if len(ins_list) > 0:
                    k = word[:result[0]]+word[result[0]-len(word)+1:]+'A'
                    errorList = [insErr+1, delErr, subErr]
                    queue.append([k,errorList])

                # Deletion Fix
                for chNum in del_list:
                    k = word[:pos]+num2dnaDict[chNum]+word[pos-len(word):]
                    k = k[:len(word)]
                    errorList = [insErr, delErr+1, subErr]
                    queue.append([k,errorList])

                # Substitution Fix
                for chNum in sub_list:
                    k = word[:pos]+num2dnaDict[chNum]+word[pos-len(word)+1:]
                    errorList = [insErr, delErr, subErr+1]
                    queue.append([k,errorList])
                    
        output.append(adjacent)
        return output
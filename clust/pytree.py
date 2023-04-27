class Node:
    def __init__(self) -> None:
        self.node_num = 4
        self.children = [None for _ in range(self.node_num)]
        self.isEnd = False
        self.singleBranch = True

class Trie:
    def __init__(self):
        self.dna_dict = {"A": 0, "T": 1, "G": 2, "C": 3}
        self.dna_dict_inv = {0: "A", 1: "T", 2: "G", 3: "C"}
        self.node_nums = 4
        self.maxOptimDepth = 2
        self.minTerminate = 14
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = False
        self.singleBranch = False

    def insert(self, word: str, label: int) -> None:
        node = self
        dict = self.dna_dict
        for i, ch in enumerate(word):
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Node()
                node = node.children[ch]
            else:
                node = node.children[ch]
                node.singleBranch = False
        node.isEnd = label

    def delete(self, word: str):
        node = self
        dict = self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = False

    def align(self, word, node, depth, max_value, memory, prev=''):
        # Exceed edit threshold
        if memory[max_value] > max_value:
            return False

        # Traverse to the end, return the cluster index and error number
        if node.isEnd != False:
            return [node.isEnd, memory[max_value], prev]

        if depth > len(word) - 1:
            return False

        chNum = self.dna_dict[word[depth]]

        # Update memory
        if depth < max_value:
            memory[max_value - depth - 1] = depth + 1
            memory[max_value + depth + 1] = depth + 1
        inp = word[max(0, depth + 1 - max_value):depth + 1]
        # If node exist
        if node.children[chNum]:
            out = prev
            out += self.dna_dict_inv[chNum]
            new_memory = self.update_memory(memory, inp, out, max_value)
            result = self.align(word, node.children[chNum], depth + 1, max_value, new_memory, out)
            # If traverse fail, try other nodes
            if result == False:
                for tmp_chNum, tmp_node in enumerate(node.children):
                    if tmp_node and tmp_chNum != chNum:
                        out = prev
                        out += self.dna_dict_inv[tmp_chNum]
                        new_memory = self.update_memory(memory, inp, out, max_value)
                        result = self.align(word, tmp_node, depth + 1, max_value, new_memory, out)
                        if result:
                            return result # Find a matching branch
                return False  # Cannnot find a branch
            else:
                return result  # Found matching branch
        else:
            for tmp_chNum, tmp_node in enumerate(node.children):
                if tmp_node:
                    out = prev
                    out += self.dna_dict_inv[tmp_chNum]
                    new_memory = self.update_memory(memory, inp, out, max_value)
                    result = self.align(word, tmp_node, depth + 1, max_value, new_memory, out)
                    if result != False:
                        return result
            return False  # Cannnot find a branch

    def update_memory(self, memory, inp, out, max_value):
        new_memory = memory.copy()
        if len(out) >= max_value:
            new_memory[0] = min(memory[0] + int(out[-1] != inp[0]), memory[1] + 1)
            new_memory[-1] = min(memory[-1] + int(out[-max_value] != inp[-1]), memory[-2] + 1)
            for i in range(1, max_value):
                delta = int(out[-1] != inp[i - 1])
                new_memory[i] = min(new_memory[i - 1] + 1, memory[i] + delta, memory[i + 1] + 1)
                j = -i - 1
                delta = int(out[-i] != inp[-1])
                new_memory[j] = min(new_memory[j + 1] + 1, memory[j] + delta, memory[j - 1] + 1)

        else:
            offset = max_value - len(out)
            for i in range(offset + 1, max_value):
                delta = int(out[-1] != inp[i - offset - 1])
                new_memory[i] = min(new_memory[i - 1] + 1, memory[i] + delta, memory[i + 1] + 1)
                j = -i - 1
                delta = int(out[i - offset - 1] != inp[-1])
                new_memory[j] = min(new_memory[j + 1] + 1, memory[j] + delta, memory[j - 1] + 1)

        delta = int(out[-1] != inp[-1])
        new_memory[max_value] = min(memory[max_value] + delta, memory[max_value - 1] + 1, memory[max_value + 1] + 1)
        return new_memory

    def clover_fuzz_align(self, word):
        """Horizontal drift function

        Args:
            word: Sequence of fuzzy retrieval.

        return:
            Returns a list with the positions that need to be drifted
            laterally and the nodes that can be drifted laterally.

        """
        node = self
        num = 0
        tmp_list = []
        fin_list = []
        len_ = len(word)
        for ch in word:

            ch = self.dna_dict[ch]
            if not node.children[ch]:

                for i in range(self.node_nums):
                    if not node.children[i]:
                        pass
                    else:
                        tmp_list.append(i)
                # print(word,ch,list,num)
                if num + 2 < len_:
                    for k in tmp_list:
                        if not node.children[k].children[self.dna_dict[word[num + 1]]]:
                            pass
                        elif not node.children[k].children[self.dna_dict[word[num + 1]]].children[self.dna_dict[word[num + 2]]]:
                            pass
                        else:
                            fin_list.append(k)
                    return [num, fin_list]
                if num + 2 == len_:
                    for k in tmp_list:
                        if not node.children[k].children[self.dna_dict[word[num + 1]]]:
                            pass
                        else:
                            fin_list.append(k)
                    return [num, fin_list]
                if num + 1 == len_:
                    for k in tmp_list:
                        fin_list.append(k)
                        return [num, fin_list]
            else:
                node = node.children[ch]
            num = num + 1

        return node.isEnd

    def clover_fuzz_fin(self, word, max_value):
        """Fuzzy search with horizontal drift.

        Args:
            word: str,Sequence of search
            max_value: int,The maximum number of horizontal drifts.

        return:
            Returns a list, the first element of which is the index of the final matched
            core sequence, and the second element is the number of horizontal drifts.
        """

        tmp_list = [[word, 0]]

        fin_list = ["", 1000]
        error_list = []
        while True:
            if tmp_list == [] or fin_list[1] == 0:
                break
            dna = tmp_list[0]

            if dna[1] > max_value:
                break
            del tmp_list[0]
            a = self.clover_fuzz_align(dna[0])
            error_list.append([dna, a])
            if type(a) != list:
                if dna[1] < fin_list[1]:
                    fin_list = [a, dna[1]]
            elif self.clover_fuzz_align(word)[1] == []:
                break
            elif a[0] == len(dna[0]) - 1:
                for i in range(len(a[1])):
                    chNum = a[1][i]
                    k = dna[0][:int(a[0])] + self.dna_dict_inv[chNum]
                    tmp_list.append([k, dna[1] + 1])
            else:
                for i in range(len(a[1])):
                    chNum = a[1][i]
                    k = dna[0][:int(a[0])] + self.dna_dict_inv[chNum] + dna[0][int(a[0]) - len(dna[0]) + 1:]
                    tmp_list.append([k, dna[1] + 1])
        fin_list.append(error_list)
        return fin_list[:2]

    def search(self, inData, max_value, tree_depth, train=False):
        word = inData[:tree_depth]
        if not train:
            # First stage: fast traversal
            result = self.clover_fuzz_fin(word, max_value)
            return result

        else:
            # Sencond stage: using edit distance
            node = self
            memory = [1000 for _ in range(int(2 * max_value + 1))]
            memory[max_value] = 0
            depth = 0
            result = self.align(word, node, depth, max_value, memory, prev='')
            if result == False:
                return ["", 1000]
            else:
                return result

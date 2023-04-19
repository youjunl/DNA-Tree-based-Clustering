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
        self.prev_num = 2
        self.maxOptimDepth = 4
        self.minTerminate = 14
        self.children = [None for _ in range(self.node_nums)]
        self.isEnd = False
        self.singleBranch = False

    def insert(self, word: str, label: str) -> None:
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

    def fuzz_align(self, word, node, depth, max_value, memory, prev=''):
        # Exceed edit threshold
        if min(memory) > max_value:
            return False

        # Traverse to the end, return the cluster index and error number
        if node.isEnd != False:
            return [node.isEnd, memory[max_value]]

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
            result = self.fuzz_align(word, node.children[chNum], depth + 1, max_value, new_memory, out)
            # If traverse fail, try other nodes
            if result == False:
                for tmp_chNum, tmp_node in enumerate(node.children):
                    if tmp_node and tmp_chNum != chNum:
                        out = prev
                        out += self.dna_dict_inv[tmp_chNum]
                        new_memory = self.update_memory(memory, inp, out, max_value)
                        result = self.fuzz_align(word, tmp_node, depth + 1, max_value, new_memory, out)
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
                    result = self.fuzz_align(word, tmp_node, depth + 1, max_value, new_memory, out)
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

    def fuzz_fin(self, inData, max_value, tree_depth=16):
        word = inData[:tree_depth]
        node = self
        memory = [1000 for _ in range(int(2 * max_value + 1))]
        memory[max_value] = 0
        depth = 0
        result = self.fuzz_align(word, node, depth, max_value, memory, prev='')
        if result == False:
            return ["", 1000]
        else:
            return result

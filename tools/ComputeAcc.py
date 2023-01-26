import os
import matplotlib.pyplot as plt
from collections import defaultdict
import getopt, sys

args = ['origin', 'f1', 'f2']

def compute(fileIn, clustNum, gamma):
    hashmap = defaultdict(int)
    falsehashmap = defaultdict(int)
    acc = [0] * len(gamma)
    with open(fileIn, 'r') as f:
        for text in f.readlines():
            text = text.strip()
            content = text.split(',')
            origin = content[0]         
            new = content[1]
            if new == origin:
                hashmap[new] += 1      
            else:
                falsehashmap[new] += 1
        f.close()
           
    for i in range(len(gamma)):
        cnt = 0
        for c in hashmap:
            if hashmap[c] >= gamma[i]*clustNum[c] and falsehashmap[c] == 0:
                cnt += 1
        acc[i] = cnt / len(clustNum)   
    return gamma, acc

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], '', '')

    if len(args) != 3:
        print('Wrong command')
        sys.exit(2)

    originFile = args[0]
    outputfile1 = args[1]
    outputfile2 = args[2]

    gamma = [i*0.02 for i in range(20, 51)]
    acc = [0] * len(gamma)
    clustNum = defaultdict(int)
    
    with open(originFile, 'r') as f:
        for text in f.readlines():
            origin = text.split(' ')[0]
            clustNum[origin] += 1

    gamma, acc = compute(outputfile1, clustNum, gamma)
    gamma, accNew = compute(outputfile2, clustNum, gamma)

    with open('CompareResult.txt', 'w') as f:
        f.truncate(0)
        for g in gamma:
            f.write('%.2f, '%g)
        f.write('\n')
        for g in acc:
            f.write('%.2f, '%g)
        f.write('\n')
        for g in accNew:
            f.write('%.2f, '%g)

    print('Gamma:')
    print(gamma)
    print('File1 Acc:')
    print(acc)
    print('File2 Acc:')
    print(accNew)

    plt.figure()    
    plt.plot(gamma, acc, 'k')
    plt.plot(gamma, accNew, 'r--')
    plt.show()
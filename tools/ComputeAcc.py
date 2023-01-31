import os
import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

args = ['origin', 'f1', 'f2']

def compute(fileIn, clustNum, gamma):
    nClust = len(clustNum.keys())

    # Read result
    results = []
    with open(fileIn, 'r') as f:
        for text in tqdm(f.readlines()):
            text = text.strip()
            content = text.split(',')
            tag = int(content[0])         
            cluster = int(content[1])
            results.append((tag, cluster))
    results = sorted(results, key=lambda k: k[1])

    maxClustNum = results[-1][1] if len(results) > 0 else 0
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append(results[i][0])

    clusters = [c for c in clusters if c != []]
    
    score = [0] * max(clustNum.keys())
    for cluster in tqdm(clusters):
        tagRef = cluster[0]
        # check identity of tags
        valid = True
        for tag in cluster:
            if tag != tagRef:
                valid = False
                break
        if valid is False:
            continue
        
        score[tagRef-1] = max(score[tagRef-1], len(cluster))

    # Compute accuracy
    gamma = [i*0.02 for i in range(20, 51)]
    acc = [0] * len(gamma)
    
    for i, g in enumerate(gamma):
        cnt = 0
        for tag in clustNum.keys():
            if score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust
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

    print('Building references')
    with open(originFile, 'r') as f:
        for text in tqdm(f.readlines()):
            origin = int(text.split(' ')[0])
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
import core
from tqdm import tqdm

if __name__ == '__main__':
    infile = 'testdata/toClustMid.txt'
    inData = []
    cnt = 0
    f = open(infile, 'r').readlines()
    clustNum = {}
    for text in tqdm(f):
        text = text.split()
        if 'N' not in text[1]:
            tag = int(text[0])
            seq = text[1]
            cluster = core.Sequence(seq, tag, cnt)
            inData.append(cluster)
            if int(tag) not in clustNum.keys():
                clustNum[int(tag)] = 1
            else:
                clustNum[int(tag)] += 1
            
            cnt += 1

    # Run clustering
    ans = core.compute_comm(inData)
    print(len(ans))

    # Sort ans
    ans = sorted(ans, key=lambda k: (k[0].tag, k[0].cluster))
    with open('output_test.txt', 'w') as f:
        f.truncate(0)
        for clust in ans:
            for seq in clust:
                f.write('%d,%d\n'%(seq.tag, seq.cluster))

    # Compute Accuaracy according to the tag
    nClust = len(clustNum.keys())
    print('Number of clusters %d'%nClust)

    score = [0] * max(clustNum.keys())
    for clust in ans:
        clustTagsRef = clust[0].tag
        # check identity of tags
        valid = True
        for seq in clust:
            if seq.tag != clustTagsRef:
                valid = False
                break
        if valid is False:
            continue
        
        score[clustTagsRef-1] = max(score[clustTagsRef-1], len(clust))

    gamma = [i*0.02 for i in range(20, 51)]
    acc = [0] * len(gamma)
    
    for i, g in enumerate(gamma):
        cnt = 0
        for tag in clustNum.keys():
            if score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust
    
    print('Gamma:')
    print(gamma)
    print('File1 Acc:')
    print(acc)
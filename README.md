# DNA Tree based Clustering

## How to use
For benchmarking these algorithms, a labeled dataset in TXT format is needed. And each line should be like:
```shell
[Label] [Read]
```
For example:
```shell
1 ATAAGGG
1 AAAAGGG
1 AAAAGGG
2 GGACCTA
2 GGACCTA
...
```

Run the clustering with command:
```shell
python -m clust.main -I [input file] -O [output file] -L [options] --no-tag
```
For example:
```shell
# test data
python -m clover.main -I testdata/toClust.txt -O output_file_1 -L 152 -P 0 --no-tag
python -m clust.main -I testdata/toClust.txt -O output_file_2 -L 152 -P 0 --no-tag
```

The output of the clustering result consists of original label and the label of cluster assigned, 
```shell
[Label],[Label of cluster]
```
for example:
```shell
1,1
1,1
1,1
2,2
2,2
...
```
For comparing the clustering result, two different metrics can be computed, with commands:
```shell
python tools/computeAcc.py [Labeled dataset] [Cluster result 1] [Cluster result 2] ... [Output file]
python tools/computePur.py [Labeled dataset] [Cluster result 1] [Cluster result 2] ... [Output file]
```
For example:
```shell
python tools/computeAcc.py testdata/toClust.txt output_file_1.txt output_file_2.txt CompareResults.txt
```
# DNA Tree based Clustering
## Introduction
This project contains code and benchmark results used in my Master thesis "Clustering for DNA Storage". The idea was developed from the tree based clustering algorithm [Clover](https://github.com/Guanjinqu/Clover). And this project provides both Python and C++ implementations for our clustering method.

## Usage
### Clustering
A read file in TXT format is needed for clustering, and each line should be like:
```shell
[Index] [Read]
```
An example of the read file:
```shell
1 ATAAGGG
2 AAAAGGG
3 AAAAGGG
4 GGACCTA
5 GGACCTA
...
```

Run the clustering with command:
```shell
python -m clust.main -I [input file] -O [output file]
```
For example:
```shell
# test data
python -m clust.main -I testdata/toClust.txt -O output_file
```

The clustering result consists of original indexes and the label of cluster assigned, 
```shell
[Index],[Label of cluster]
```
An example of clustering result:
```shell
1,1
2,1
3,1
4,2
5,2
...
```
### Benchmark
For comparing the clustering result, a TXT file that indeicates accuracte clustering index is need:
```shell
[Index],[Label of cluster]
```
An example of :
```shell
1,1
2,1
3,1
4,2
5,2
...
```

Two different metrics can be computed with commands:
```shell
python tools/computeAcc.py [Accurate indexes] [Cluster result 1] [Cluster result 2] ... [Output file]
python tools/computePur.py [Accurate indexes] [Cluster result 1] [Cluster result 2] ... [Output file]
```
For example:
```shell
python tools/computeAcc.py testdata/toClust.txt output_file_1.txt output_file_2.txt CompareResults.txt
```

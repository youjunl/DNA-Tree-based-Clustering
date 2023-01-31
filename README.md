# DNA Tree based Clustering
 
## How to use
python -m clover.main -I testdata/toClust.txt -O output_file_1 -L 152 -P 0 --no-tag
python -m clust.main -I testdata/toClust.txt -O output_file_2 -L 152 -P 0 --no-tag
python tools/computeAcc.py testdata/toClust.txt output_file_1.txt output_file_2.txt

## For testing on small size test data
python -m clover.main -I testdata/toClustSmall.txt -O output_file_1 -L 152 -P 0 --no-tag
python -m clust.main -I testdata/toClustSmall.txt -O output_file_2 -L 152 -P 0 --no-tag
python tools/computeAcc.py testdata/toClustSmall.txt output_file_1.txt output_file_2.txt
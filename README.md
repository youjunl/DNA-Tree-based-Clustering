# DNA Tree based Clustering
 
## How to use
<<<<<<< HEAD
python -m clover.main -I testdata/toClust.txt -O output_file_1 -L 152 -P 0 --no-tag
python -m clust.main -I testdata/toClust.txt -O output_file_2 -L 152 -P 0 --no-tag
python tools/computeAcc.py testdata/toClust.txt output_file_1.txt output_file_2.txt

## For testing on small size test data
python -m clover.main -I testdata/toClustSmall.txt -O output_file_1 -L 152 -P 0 --no-tag
python -m clust.main -I testdata/toClustSmall.txt -O output_file_2 -L 152 -P 0 --no-tag
python tools/computeAcc.py testdata/toClustSmall.txt output_file_1.txt output_file_2.txt
=======
python -m clover.main -I toClust.txt -O output_file -L 152 -P 0 --no-tag

python -m clust.main -I toClust.txt -O output_file_new -L 152 -P 0 --no-tag

python computeAcc.py
>>>>>>> ef2568e431a90a0bbce784f7c6e1fc949676a8d8

# Engorgio: An Arbitrary-Precision Unbounded-Size Hybrid  Encrypted Database via Quantized Fully Homomorphic Encryption
This code is the implementation of the paper Engorgio: An Arbitrary-Precision Unbounded-Size Hybrid  Encrypted Database via Quantized Fully Homomorphic Encryption.

## Requirements

```
git 
gcc >= 10
cmake >= 3.16
GMP 6.2.0 
```

## Building Engorgio
You can build the Engorgio (out-of-source) for your machine by executing the following commands:
```
mkdir build
cd build
cmake ..
make
```


Following the instructions above, output binaries will be generated in the `build/bin/` directory. You can run these binaries by:
```
$./bin/comparison_test
$./bin/sort_test
$./bin/sync_test
$./bin/topk_test
$./bin/relational_query_test
$./bin/vectorized_query_test_1
$./bin/vectorized_query_test_2
$./bin/hybrid_query_test_1
$./bin/hybrid_query_test_2
```

## Examples
### Homomorphic Comparison
- codes `test/comparison_test.cpp`  
- output binary `build/bin/comparison_test` 
- This demo shows the homomorphic comparison operator `ArbQuantComp` in Engorgio. Execute the binary file `comparison_test`, we will get the latencry of Engorgio in fig 5a and fig 5b. 

### Homomorphic Sorting
- codes `test/sort_test.cpp`  
- output binary `build/bin/sort_test` 
- This demo shows the homomorphic sorting operator `HomSort` in Engorgio. Execute the binary file `sort_test`, we will get the latencry of Engorgio in fig 5c. Although we have reduced the number of experiments, it still takes a long time (20+ hours) to complete the experiment.

### Homomorphic Synchronization
- codes `test/sync_test.cpp`  
- output binary `build/bin/sync_test` 
- This demo shows the homomorphic synchronization operator `HomSync` in Engorgio. Execute the binary file `sync_test`, we will get the latencry of Engorgio in fig 5d. 

### Homomorphic Top-k
- codes `test/topk_test.cpp`  
- output binary `build/bin/topk_test` 
- This demo shows the homomorphic Topk operator `HomTopK` in Engorgio. Execute the binary file `topk_test`, we will get the latencry of Engorgio in fig 5e. 

### Relational Query Evaluation
- codes `test/relational_query_test.cpp`  
- output binary `build/bin/relational_query_test` 
- This demo evaluates relational queries on an encrypted database with varying numbers of rows. Execute the binary file `relational_query_test`, we will get the latencry of Engorgio in fig 6a and fig 6b. 

### Vectorized Query Evaluation
- codes `test/vectorized_query_test_1.cpp` and `test/vectorized_query_test_2.cpp`
- output binary `build/bin/vectorized_query_test_1` and `build/bin/vectorized_query_test_2`
- This demo evaluates vectorized queries on an encrypted database with varying numbers of rows. Execute the binary files `vectorized_query_test_1` and `vectorized_query_test_2`, we will get the latencry of Engorgio in fig 6c and fig 6d. Although we have reduced the number of experiments, it still takes a long time (20+ hours) to complete the experiment.

### Hybrid Query Evaluation
- codes `test/hybrid_query_test_1.cpp` and `test/hybrid_query_test_2.cpp` 
- output binary `build/bin/hybrid_query_test_1` and `build/bin/hybrid_query_test_2`
- This demo evaluates hybrid queries on an encrypted database with varying numbers of rows. Execute the binary files `hybrid_query_test_1` and `hybrid_query_test_2`, we will get the latencry of Engorgio in fig 6e.  Although we have reduced the number of experiments, it still takes a long time (20+ hours) to complete the experiment.

## Citation
```
@inproceedings{engorgio,
 author = {Song Bian and Haowen Pan and Jiaqi Hu and Zhou Zhang and Yunhao Fu and Jiafeng Hua and Yi Chen and Bo Zhang and Yier Jin and Jin Dong and Zhenyu Guan},
 title = {Engorgio: An Arbitrary-Precision Unbounded-Size Hybrid Encrypted Database via Quantized Fully Homomorphic Encryption},
 booktitle = {{USENIX} Security Symposium ({USENIX} Security)},
 year = {2025}
}
```
# Description
Technical report and implementation code for the paper "Efficient Algorithm for Budgeted Adaptive Influence
Maximization: An Incremental RR-set Update Approach" in SIGMOD 24

# Guideline for running code
0. Compile code: 
```
   cd code
   make
```

1.  format graph and generate realization
```
   ./format_graph dataset/Epinions
   ./gene_real dataset/Epinions 20
```


2. Run sub-policies on dataset Epinions
```shell
# Running CaGreedy with incremental update strategy under degree-based cost model
./baim -dataset dataset/Epinions -model IC -seedfile seed -time 20 -budget 100 -algo 1

# Running CaGreedy with incremental update strategy under random cost model
./baim -dataset dataset/Epinions -model IC -seedfile seed -time 20 -budget 100 -algo 1 -randomcost 1

# Running CaGreedy without incremental update strategy under degree-based cost model
./baim -dataset dataset/Epinions -model IC -seedfile seed -time 20 -budget 100 -algo 1 -reuse 0

# Running MIS under degree-based cost model
./baim -dataset dataset/Epinions -model IC -seedfile seed -time 20 -budget 100 -algo 3

# Running IMAGE with epsilon 0.5 under degree-based cost model
./baim -dataset dataset/Epinions -model IC -seedfile seed -time 20 -budget 100 -algo 2 -beta 0.5
```

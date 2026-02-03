conda activate julia
julia wrap_gibbs.jl test
directory structure (after running example)
├── data
│   └── test_gibbs_info.RDS   ### readRDS this to see the input data format.  
├── output
│   ├── test_branch_VAFs.txt.gz
│   └── test_posterior_VAFs.txt.gz
├── README.txt
├── src
│   └── Deep_seq_tree_GS.jl
└── wrap_gibbs.jl

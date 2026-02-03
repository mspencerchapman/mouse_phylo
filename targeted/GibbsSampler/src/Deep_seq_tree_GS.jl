using Phylo
using RCall
using DataFrames
using Distributions
using Random
using CodecZlib

mutable struct Mutation
    ID::String
    obs_var::Int64
    obs_ref::Int64
    obs_depth::Int64
    seq_error::Float64
    posterior_VAF::Array{Float64, 1}
end


# Function for splitting branches of tree into each mutation with deep sequencing data
function split_tree_by_mut(tree; FUN = x->x.obs_var / x.obs_depth)
    new_tree = RootedTree()
    intensity = Vector{Float64}()
    for i in traversal(tree, preorder)
        curr_nodename = getnodename(tree, i)
        if isroot(tree,i)
            createnode!(new_tree, curr_nodename)
            push!(intensity, 0)
            continue
        end
        muts_on_i = getbranchdata(tree, getinbound(tree,i))["Muts"]
        old_br_l = getlength(tree, getinbound(tree,i))
        parent_nodename = getnodename(tree, getparent(tree,i))
        if length(muts_on_i) == 0
            push!(intensity, 0)
            createnode!(new_tree, curr_nodename)
            createbranch!(new_tree, parent_nodename, curr_nodename, old_br_l)
            continue
        end
        new_br_l = old_br_l / length(muts_on_i)
        vafs = sort(map(x -> x => FUN(x), muts_on_i), by=x->x[2], rev=true)
        for (j,k) in enumerate(vafs)
            push!(intensity, k[2])
            new_nodename = j == length(muts_on_i) ? curr_nodename : "$curr_nodename $(k[1].ID)"
            createnode!(new_tree, new_nodename)
            createbranch!(new_tree, parent_nodename, new_nodename, new_br_l)
            parent_nodename = new_nodename
        end
    end
    # plot(new_tree, treetype = :dendrogram, line_z = intensity, linecolor = :RdYlBu,
    #        linewidth = 5, showtips = false)
    return (new_tree, intensity)
end


function deep_seq_GS(tree, start_VAFs, end_VAFs, muts, iter, burn_in, thin;
    scale_pm = 100, min_VAF = 10.0^-10)
    # Create the variables for tracking output of Gibbs sampler and other intermediates
    curr_VAFs = Dict{LinkBranch, Dict{String, Float64}}()
    br_starts = Dict{LinkBranch, Float64}()
    br_ends = Dict{LinkBranch, Float64}()
    
    # Set initial conditions - split values for internal nodes
    for i in traversal(tree, postorder)
        if isroot(tree,i)
            for j in getoutbounds(tree,i)
                br_starts[j] = 0.5
            end
        elseif isleaf(tree,i)
            br_ends[getinbound(tree,i)] = min_VAF
        else
            for j in getoutbounds(tree,i)
                br_starts[j] = br_ends[j] + min_VAF
            end
            br_ends[getinbound(tree,i)] = sum(map(x->br_starts[x],getoutbounds(tree,i))) + min_VAF
        end
    end
    for i in branchiter(tree)
        curr_VAFs[i] = Dict{String, Float64}()
        for j in muts[i]
            curr_VAFs[i][j.ID] = max(min_VAF, (br_starts[i] + br_ends[i])/2)
        end
    end
    
    # Go through iterations
    for i in 1:iter
        # Block 1 - update VAFs of individual mutations
        # Proceed branch by branch
        for j in branchiter(tree)
            for k in muts[j]
                # Get proposed new VAF from proposal distribution
                old_VAF = curr_VAFs[j][k.ID]
                new_VAF = rand(truncated(Beta(old_VAF*scale_pm/(1-old_VAF),scale_pm), br_ends[j], br_starts[j]))
                if new_VAF<=0
                    new_VAF=min_VAF
                end
                # Calculate the acceptance ratio
                old_p = old_VAF + k.seq_error - 2 * old_VAF * k.seq_error
                new_p = new_VAF + k.seq_error - 2 * new_VAF * k.seq_error
                v1=new_VAF * scale_pm / (1-new_VAF)
                v2=old_VAF * scale_pm / (1-old_VAF)
                #if v2<=0
                #    println("v2:",v2," ",old_VAF)
                #end
                log_AR = logpdf(Binomial(k.obs_depth, new_p), k.obs_var) -
                        logpdf(Binomial(k.obs_depth, old_p), k.obs_var) +
                        logpdf(truncated(Beta(v1, scale_pm), br_ends[j], br_starts[j]), old_VAF) -
                        logpdf(truncated(Beta(v2, scale_pm), br_ends[j], br_starts[j]), new_VAF)

                # Test the acceptance ratio
                if rand(Uniform(0,1)) <= exp(log_AR)
                    if new_VAF< min_VAF
                        println("Forcing minimum VAF for:")
                        print(k)
                        println("")
                        new_VAF=min_VAF
                    end
                    curr_VAFs[j][k.ID] = new_VAF
                end
            end
        end
        # Block 2 - update branch start and end VAFs
        for j in traversal(tree, preorder)
            if isroot(tree,j) || isleaf(tree,j)
                continue
            end
            max_VAFs_out = Vector{Float64}()
            for k in getoutbounds(tree, j)
                push!(max_VAFs_out, isempty(curr_VAFs[k]) ? br_ends[k] : maximum(values(curr_VAFs[k])))
            end
            min_VAF_in = isempty(values(curr_VAFs[getinbound(tree,j)])) ?
                            br_starts[getinbound(tree,j)] : minimum(values(curr_VAFs[getinbound(tree,j)]))
            unallocated_VAF = sort(rand(Uniform(0, min_VAF_in - sum(max_VAFs_out)), length(max_VAFs_out)+2))
            unallocated_VAF = diff([unallocated_VAF;1])
            br_ends[getinbound(tree,j)] = min_VAF_in - unallocated_VAF[1]
            for (k,l) in enumerate(getoutbounds(tree,j))
                br_starts[l] = max_VAFs_out[k] + unallocated_VAF[k+1]
            end
        end

        # Now update posterior samples data if needed
        if i % thin == 0
            println(i)
            flush(stdout)
            if i > burn_in
                for j in branchiter(tree)
                    push!(start_VAFs[j], br_starts[j])
                    push!(end_VAFs[j], br_ends[j])
                    for k in muts[j]
                        push!(k.posterior_VAF, curr_VAFs[j][k.ID])
                    end
                end
            end
        end
    end
    return start_VAFs, end_VAFs, muts
end


function truncate_tree(tree, max_height)
    new_tree = RootedTree()
    for i in traversal(tree)
        curr_nodename = getnodename(tree, i)
        if isroot(tree,i)
            createnode!(new_tree, curr_nodename)
            continue
        end

        old_br_l = getlength(tree, getinbound(tree,i))
        parent_nodename = getnodename(tree, getparent(tree,i))

        if getheight(tree, i) <= max_height
            # Node contained within range to be retained
            createnode!(new_tree, curr_nodename)
            createbranch!(new_tree, parent_nodename, curr_nodename, old_br_l,
                            data = getbranchdata(tree, getinbound(tree,i)))
        elseif getheight(tree, parent_nodename) <= max_height
            # Node in range to be truncated but parent retained
            createnode!(new_tree, curr_nodename)
            createbranch!(new_tree, parent_nodename, curr_nodename,
                            max_height - getheight(tree, parent_nodename),
                            data = getbranchdata(tree, getinbound(tree,i)))
        #else
            # Node in range to be truncated and parent already truncated
        end
    end

    return new_tree
end


# Function for splitting branches of donor and recipient pairs into each mutation with deep sequencing data
function split_tree_pairs_by_mut(tree_dict; SORTFUN = x->median(x.posterior_VAF),
                                    PLOTFUN = (x,y)->median(y.posterior_VAF) - median(x.posterior_VAF))
    new_tree = RootedTree()
    intensity = Vector{Float64}()
    for ix in eachindex(traversal(tree_dict["Donor"], preorder))
        i = traversal(tree_dict["Donor"], preorder)[ix]
        h = traversal(tree_dict["Recipient"], preorder)[ix]
        curr_nodename = getnodename(tree_dict["Donor"], i)
        if isroot(tree_dict["Donor"],i)
            createnode!(new_tree, curr_nodename)
            push!(intensity, 0)
            continue
        end
        muts_on_i = getbranchdata(tree_dict["Donor"], getinbound(tree_dict["Donor"],i))["Muts"]
        muts_on_i_recip = getbranchdata(tree_dict["Recipient"], getinbound(tree_dict["Recipient"],h))["Muts"]
        old_br_l = getlength(tree_dict["Donor"], getinbound(tree_dict["Donor"],i))
        parent_nodename = getnodename(tree_dict["Donor"], getparent(tree_dict["Donor"],i))
        if length(muts_on_i) == 0
            push!(intensity, 0)
            createnode!(new_tree, curr_nodename)
            createbranch!(new_tree, parent_nodename, curr_nodename, old_br_l)
            continue
        end
        new_br_l = old_br_l / length(muts_on_i)

        # Define ordering of mutations
        # Pick the donor or recipient pair by the greatest range
        donor_FUN = map(x -> x => SORTFUN(x), muts_on_i)
        recip_FUN = map(x -> x => SORTFUN(x), muts_on_i_recip)

        sort_donor_FUN = sortperm(donor_FUN, by=x->x[2], rev=true)
        sort_recip_FUN = sortperm(recip_FUN, by=x->x[2], rev=true)
        sortorder = maximum(last.(donor_FUN)) - minimum(last.(donor_FUN)) >
                        maximum(last.(recip_FUN)) - minimum(last.(recip_FUN)) ?
                        sort_donor_FUN : sort_recip_FUN

        # Now run the function for generating the intensity values on the sorted mutations
        pair_FUN = map((x,y) -> x => PLOTFUN(x,y), muts_on_i[sortorder], muts_on_i_recip[sortorder])

        for (j,k) in enumerate(pair_FUN)
            push!(intensity, k[2])
            new_nodename = j == length(muts_on_i) ? curr_nodename : "$curr_nodename $(k[1].ID)"
            createnode!(new_tree, new_nodename)
            createbranch!(new_tree, parent_nodename, new_nodename, new_br_l)
            parent_nodename = new_nodename
        end
    end
    # plot(new_tree, treetype = :dendrogram, line_z = intensity, linecolor = :RdYlBu,
    #        linewidth = 5, showtips = false)
    return (new_tree, intensity)
end


function write_GS_output(tree, GS_output, dir, file_stem, node_to_branch)
    # Reverse node_to_branch Dict
    branch_to_node = Dict(j => i for (i,j) in node_to_branch)

    # Write posterior samples of VAFs for each mutation to file
    # File format is Node_assignment \t mutation_ID \t posterior_VAFs (comma-separated)
    vaf_out = CodecZlib.GzipCompressorStream(open("$dir/$(file_stem)_posterior_VAFs.txt.gz", "w"); level = 1)
    println(vaf_out, "Node_assignment\tmutation_ID\tPosterior_VAFs")
    for (i,j) in GS_output[3]
        for k in j
            println(vaf_out,
                join([branch_to_node[i], k.ID, join(k.posterior_VAF, ",")], "\t"))
        end
    end
    close(vaf_out)

    # Write branch start and end VAFs
    # File format is Node_assignment \t Type (either 'Top_VAF' or 'Bottom_VAF') \t posterior_VAFs (comma-separated)
    branch_vaf_out = CodecZlib.GzipCompressorStream(open("$dir/$(file_stem)_branch_VAFs.txt.gz", "w"); level = 1)
    println(branch_vaf_out, "Node_assignment\tType\tPosterior_VAFs")
    for (i,j) in GS_output[1]
        println(branch_vaf_out,
            join([branch_to_node[i], "Top_VAF", join(j, ",")], "\t"))
        println(branch_vaf_out,
            join([branch_to_node[i], "Bottom_VAF", join(GS_output[2][i], ",")], "\t"))
    end
    close(branch_vaf_out)

end


"""
read_GS_output!(tree, mut_vafs_file, branch_vafs_file, start_VAF, end_VAF, muts)

Reads data from files produced by `write_GS_output()` and places them into `tree` data

`tree` is a phylogenetic tree with the branch data partially filled in
    (ie mutations already assigned to branches) - this function just adds the
    posterior VAFs etc
`muts` is populated Dict of branches to array of mutations (missing posterior VAFs but otherwise complete)
`start_VAF` and `end_VAF` are empty Dicts of branches to arrays
"""
function read_GS_output!(tree, mut_vafs_file, branch_vafs_file, start_VAF, end_VAF, muts)
    # Define a mapping for the nodes assigned in the VAF files to the branch names of tree
    node_to_branch = Dict{Int64, LinkBranch}()
    for i in 1:nleaves(tree)
        node_to_branch[i] = getinbound(tree, getleaves(tree)[i])
    end
    for i in (nleaves(tree)+2):nnodes(tree) # Skipping root, which in R is always nleaves+1
        node_to_branch[i] = getinbound(tree, "Node $i")
    end

    mutation_locator = Dict{Array{Any,1}, Int}()
    for i in branchiter(tree)
        for (j,k) in enumerate(muts[i])
            mutation_locator[[i; k.ID]] = j
        end
    end

    # Read in data from mutation VAFs file
    vaf_in = CodecZlib.GzipDecompressorStream(open(mut_vafs_file, "r"))
    for line in eachline(vaf_in)
        if startswith(line, "Node_assignment") # Skip header line
            continue
        end
        line_spl = split(line, '\t')
        curr_branch = node_to_branch[parse(Int, line_spl[1])]
        curr_pos = mutation_locator[[curr_branch, line_spl[2]]]
        muts[curr_branch][curr_pos].posterior_VAF =
            parse.(Float64, split(line_spl[3], ','))
    end

    # Read in data from branch VAFs file
    branch_vaf_in = CodecZlib.GzipDecompressorStream(open(branch_vafs_file, "r"))
    for line in eachline(branch_vaf_in)
        if startswith(line, "Node_assignment") # Skip header line
            continue
        end
        line_spl = split(line, '\t')
        curr_branch = node_to_branch[parse(Int, line_spl[1])]
        if line_spl[2] == "Top_VAF"
            start_VAF[curr_branch] = parse.(Float64, split(line_spl[3], ','))
        elseif line_spl[2] == "Bottom_VAF"
            end_VAF[curr_branch] = parse.(Float64, split(line_spl[3], ','))
        end
    end
end

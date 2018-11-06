
using JLD2
using Distributions
using StatsBase

"""
    COOccurence(coo_matrix, labels)
"""
type COOccurence
    "coo_matrix"
    coo_matrix::SparseMatrixCSC{Float64,Int64}
    "ind2labels"
    ind2labels::Dict{Int, String}
end

# Note: Sample input occurrance files are those generated with PubMedMiner.jl
# https://github.com/bcbi/PubMedMiner.jl.git
function occur2coo(occur, label2ind)
    #invert map
    ind2label = Dict()
    for key  in keys(label2ind)
        val = label2ind[key]
        ind2label[val] = key
    end
    # labels = Array{String}(length(ind2label))
    # for i in range(1,length(ind2label))
    #     labels[i] = ind2label[i]
    # end

    coo_m = occur* occur.'
    COO = COOccurence(coo_m, ind2label)
    return COO
end

function read_sp_occur(file="./occur_sp.jdl")
    fin  = jldopen(file, "r")
    occur = read(fin, "occur")
    return occur
end

function read_occur_label_dict(file="/mesh2ind.jdl")
    file  = jldopen(file, "r")
    label2ind = read(file, "mesh2ind")
    return label2ind
end

"""
    corrcoef(occur)

Compute correlation coeffiecients for a occurence data matrix
###Arguments

* `occur`: Occurrence data matrix
"""
function corrcoef(occur; vardim = 1)
    # compute the covariance matrix
    C = cov(occur, vardim)

    # the correlation coefficients are given by
    # C_{i,j} / sqrt(C_{ii} * C_{jj})
    d = diag(C)
    CiiCjj = zeros(size(C))
    BLAS.ger!(1.0, d, d, CiiCjj) #CiiCjj += d*d'

    #should we devectorize?
    coeffs = C./sqrt.(CiiCjj)

    return coeffs

end

"""
    pmi(coo_matrix)

Point-wise mutual information from co-occurrence  matrix

Point-wise mutual information
log p(x,y)/p(x)p(y)

A measure the takes into account chance by comparing the probability of
joint event to the product of the individual events


###Outcome

* `pmi`:PMI matrix (symmetric matrix)
"""
function pmi_mat(coo_matrix, N)
    d = diag(coo_matrix)

    #Compute conditionals P(X,Y)/P(X)
    cooccurrence_diagonal = diag(coo_matrix)
    conditional = coo_matrix./cooccurrence_diagonal

    # PMI (point-wise mutual information) P(X,Y)/P(X)P(Y)
    pmi_full = log.((N .* conditional') ./ cooccurrence_diagonal)
    # pmi = LowerTriangular(pmi_full)
end

"""
    chi2(coo_matrix)

Point-wise Chi square statistic from co-occurrence  matrix

Chi square statistic to test whether the occurence of one
feature is independent of another, e.i
Null hypothesis = p(x,y) = p(x)p(y)

###Arguments
* `occur`: Occurrance data matrix: Observations are rows, columns are features

###Outcome

* `chi2`:Matrix of Chi Square Statistics  (symmetric matrix)
"""
function chi2_mat(occur; min_freq=5)

    nvars = size(occur, 2)
    chi2 = LowerTriangular(zeros(nvars, nvars))
    for j1=1:nvars
        x = occur[:,j1]
        for j2=j1+1:nvars
            y = occur[:,j2]
            chi2[j2, j1] = chi2_statistic(x,y,span(occur); min_freq=min_freq)
        end
    end
    return chi2
end

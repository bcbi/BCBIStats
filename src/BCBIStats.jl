module BCBIStats

#CO-Ocurrence module
module COOccur

export power_divergence_statistic,
       chi2_statistic,
       likelihood_ratio,
       mod_likelihood_ratio,
       neyman_statistic

include("COOccur/contingency_statistics.jl")

export COOccurence,
       occur2coo,
       read_sp_occur,
       read_occur_label_dict,
       corrcoef,
       pmi

include("COOccur/COOccur.jl")



end



end # module

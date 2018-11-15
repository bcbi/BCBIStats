using BCBIStats.COOccur
using StatsBase

#Co-occurrence matrix
occur_file = "/Users/isa/Dropbox/BrownAgain/Projects/BCBI/test/occur_sp.jdl"
labels_file = "/Users/isa/Dropbox/BrownAgain/Projects/BCBI/test/mesh2ind.jdl"

#Contingency
x=[1,0,0,1,0,1,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,1,0,1]
y=[1,1,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0]

occur = [x y]
label2ind = Dict("x" => 1, "y" => 2)

try
  global occur = read_sp_occur(occur_file)
  global label2ind = read_occur_label_dict(labels_file)
catch
  nothing
end

COO = occur2coo(occur, label2ind)

corr = corrcoef(occur)
pmi = pmi_mat(COO.coo_matrix, sum(occur))



d = counts(x,y, span(occur))
stat = chi2_statistic(x,y,span(occur))
@test isapprox(stat,  0.4057853910795092)
stat = chi2_statistic(d)
@test isapprox(stat, 0.4057853910795092)
stat = chi2_statistic(d, min_freq=6)
@test isnan(stat)

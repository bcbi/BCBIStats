using BCBIStats.COOccur
using StatsBase

#Co-occurrence matrix
occur_file = "/Users/isa/Dropbox/BrownAgain/Projects/BCBI/test/occur_sp.jdl"
labels_file = "/Users/isa/Dropbox/BrownAgain/Projects/BCBI/test/mesh2ind.jdl"

occur = read_sp_occur(occur_file)
label2ind = read_occur_label_dict(labels_file)

COO = occur2coo(occur, label2ind)

corr = corrcoef(occur)
pmi = pmi(COO.coo_matrix, sum(occur))

#Contingency
x=[1,0,0,1,0,1,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,1,0,1]
y=[1,1,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0]

occur = [x y]

d = counts(x,y, span(occur))
stat = chi2_statistic(x,y,span(occur))
@test_approx_eq stat  0.4057853910795092
stat = chi2_statistic(d)
@test_approx_eq stat  0.4057853910795092
stat = chi2_statistic(d, min_freq=6)
@test isnan(stat)

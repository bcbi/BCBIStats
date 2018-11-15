# power divergence statistic for testing goodness of fit
# Cressie and Read 1984; Read and Cressie 1988
# lambda  =  1: Pearson's Chi-square statistic
# lambda ->  0: Converges to Likelihood Ratio test stat
# lambda -> -1: Converges to minimum discrimination information statistic (Gokhale, Kullback 1978)
# lambda =  -2: Neyman Modified chi-squared statistic (Neyman 1949)
# lambda = -.5: Freeman-Tukey statistic (Freeman, Tukey 1950)
# Mirror functionality of https://github.com/JuliaStats/HypothesisTests.jl/blob/master/src/power_divergence.jl

# Under regularity conditions, their asymptotic distributions are all the same (Drost 1989)
# Chi-squared null approximation works best for lambda near 2/3
function power_divergence_statistic(x::Matrix{T};
                                                lambda=1, min_freq=5) where T<:Integer

  nrows, ncols = size(x)
  n = sum(x)

  #validate
  any( x .< 0) || any( .!isfinite.(x)) ? error("all entries must be nonnegative and finite") : nothing
  n == 0 ? error("at least one entry must be positive") : nothing
  isfinite(nrows) && isfinite(ncols) && isfinite(nrows*ncols) ? nothing : error("Invalid number of rows or columns")



  if nrows > 1 && ncols > 1
    rowsums = mapslices(sum, x, dims=2)
    colsums = mapslices(sum, x, dims=1)
    df = (nrows-1)*(ncols-1)
    xhat = rowsums * colsums / n
    if( minimum(x) < min_freq )
        @warn "Min frequecy requirement violated - returning NaN:" minimum(x)
        return NaN
    end
   else
    error("Number of rows and columns must be greater than 1")
  end

  stat = 0
  if lambda == 0
    for i in 1:length(x)
      stat += x[i]*(log(x[i]) - log(xhat[i]))
    end
    stat *= 2
  elseif lambda == -1
    for i in 1:length(x)
      stat += xhat[i]*(log(xhat[i]) - log(x[i]))
    end
    stat *= 2
  else
    for i in 1:length(x)
      stat += x[i]*( (x[i]/(xhat[i]))^lambda - 1)
    end
    stat *= 2/(lambda*(lambda+1))
  end

  return stat
end


#Chi2 Statictic
function chi2_statistic(x::AbstractMatrix{T}; min_freq=5) where T<:Integer
  power_divergence_statistic(x,lambda=1.0, min_freq=min_freq)
end

function chi2_statistic(x::AbstractVector{T}, y::AbstractVector{T},
                        levels::UnitRange{T}; min_freq=5) where T<:Integer
  d = counts(x,y,levels)
  power_divergence_statistic(d,lambda=1.0, min_freq=min_freq)
end

#Likelihood ratio (G2) Statictic
function likelihood_ratio(x::AbstractMatrix{T}; min_freq=5) where T<:Integer
  power_divergence_statistic(x,lambda=0, min_freq=min_freq)
end

function likelihood_ratio(x::AbstractVector{T}, y::AbstractVector{T},
                                      levels::UnitRange{T}; min_freq=5) where T<:Integer
  d = counts(x,y,levels)
  power_divergence_statistic(d,lambda=0, min_freq=min_freq)
end

#Minimum discrimination information
function mod_likelihood_ratio(x::AbstractMatrix{T}; min_freq=5) where T<:Integer
  power_divergence_statistic(x,lambda=-1, min_freq=min_freq)
end

function mod_likelihood_ratio(x::AbstractVector{T}, y::AbstractVector{T},
                                          levels::UnitRange{T}; min_freq=5) where T<:Integer
  d = counts(x,y,levels)
  power_divergence_statistic(d,lambda=-1, min_freq=min_freq)
end

#Neyman's Statistic
function neyman_statistic(x::AbstractMatrix{T}; min_freq=5) where T<:Integer
  power_divergence_statistic(x,lambda=-2, min_freq=min_freq)
end

function neyman_statistic(x::AbstractVector{T}, y::AbstractVector{T},
                                          levels::UnitRange{T}; min_freq=5) where T<:Integer
  d = counts(x,y,levels)
  power_divergence_statistic(d,lambda=-2, min_freq=min_freq)
end

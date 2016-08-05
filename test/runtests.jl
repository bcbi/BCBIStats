using BCBIStats
using Base.Test

# write your own tests here
@test 1 == 1

my_tests = ["cooccur.jl"]

println("Running tests:")

for my_test in my_tests
    @printf " * %s\n" my_test
    include(my_test)
end

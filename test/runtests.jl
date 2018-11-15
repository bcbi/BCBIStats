using BCBIStats
using Test

my_tests = ["cooccur.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * %s\n", my_test)
    include(my_test)
end

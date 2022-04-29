using DifferentialAlgebra.Test
using DifferentialAlgebra.TestSetExtensions

using DifferentialAlgebra

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end


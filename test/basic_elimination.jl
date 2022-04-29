@testset "Basic elimination" begin
    coef_ring, (a, b, c, dd) = DifferentialAlgebra.Nemo.PolynomialRing(DifferentialAlgebra.Nemo.QQ, ["a", "b", "c", "d"])
    F = DifferentialAlgebra.AbstractAlgebra.FractionField(coef_ring)

    R, (x, y) = DifferentialPolynomialRing(F, ["x", "y"])
    a, b, c, dd = map(x -> R(x), [a, b, c, dd])

    eqs = [
        d(x) - a * x + b * x * y,
        d(y) + c * y - dd * x * y
    ]
    for i in 1:2
        push!(eqs, d(eqs[i]))
    end
    Rsing, pls = get_singular_polys(eqs, :lex, lex_ranking(["y", "x"]))
    gb = DifferentialAlgebra.Singular.gens(DifferentialAlgebra.Singular.std(DifferentialAlgebra.Singular.Ideal(Rsing, pls)))
    ioequation = a * dd * from_singular(gb[1], R)
    @test ioequation == a * dd * x^3 - dd * x^2 * d(x) - a * c * x^2 + c * x *d(x, 1) + x * d(x, 2) - d(x, 1)^2
end

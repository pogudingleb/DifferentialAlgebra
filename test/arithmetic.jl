using LinearAlgebra

@testset "Basic arithmetic" begin
    R, (x, y) = DifferentialPolynomialRing(DifferentialAlgebra.Nemo.QQ, ["x", "y"])
    @test x * d(y) + y * d(x) - d(x * y) == 0
    @test d(x^2) == 2 * x * d(x)

    R, (x,) = DifferentialPolynomialRing(DifferentialAlgebra.Nemo.QQ, ["x"])
    f = d(x)^2 - x
    g = d(x, 2) - 1//2
    @test diffred(g, f)[1] == 0

    f = d(x)^2 - 4 * x^3 - 1
    g = d(x, 2) - 6 * x^2
    @test diffred(g, f)[1] == 0

    f = d(x)^2 - x^2
    g_a = d(x, 2) - x
    g_b = det(wronskian([x, d(x)]))
    @test diffred(g_a, f)[1] == 0

    R, (t, y1, y2) = DifferentialPolynomialRing(DifferentialAlgebra.Nemo.QQ, ["t", "y1", "y2"])
    @test det(wronskian([t * y1, t * y2])) - t^2 * det(wronskian([y1, y2])) == 0

    R, (t, y1, y2, y3) = DifferentialPolynomialRing(DifferentialAlgebra.Nemo.QQ, ["t", "y1", "y2", "y3"])
    @test det(wronskian([t * y1, t * y2, t * y3])) - t^3 * det(wronskian([y1, y2, y3])) == 0

    R, (x,) = DifferentialPolynomialRing(DifferentialAlgebra.Nemo.QQ, ["x"])
    @test numerator(evaluate(d(x, 2) - x, [1 // x])) == -d(x, 2) * x + 2 * d(x, 1)^2 - x^2
end

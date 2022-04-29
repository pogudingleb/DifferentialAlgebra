module DifferentialAlgebra

import Base

using AbstractAlgebra
using Nemo
using Singular

# for testing
using Test
using TestSetExtensions

export DifferentialPolynomialRing, DiffPolyRing, DiffPoly, d;
export evaluate;
export order, separant, leader, diffred;
export wronskian;
export get_singular_polys, from_singular, lex_ranking;

include("utils.jl")
include("abstract_types.jl")
include("differential_polynomials.jl")
include("algebraic_ideals.jl")

end

# max x₁+x₂ s.t.
#  x₁ + 2x₂ ≤ 4
# 4x₁ + 2x₂ ≤ 12
# -x₁ +  x₂ ≤ 1
# and xᵢ ≥ 0 ∀ i
a = [
     1  2  1  0  0  4
     4  2  0  1  0 12
    -1  1  0  0  1  1
     1  1  0  0  0  0
]
# should give 10/3 with x = (8/3, 2/3)

simplex_max(a::Array{<:Integer}) =
    simplex_max(Rational.(a))

simplex_max(a) = _simplex(a, findmax, x -> x .> 0)
simplex_min(a) = _simplex(a, findmin, x -> x .< 0)

"""
solve a linear problem with the basic simplex algorithm
TODO doesn't handle unbounded solution or empty solution spaces
TODO currently the initial matrix must represent a feasible solution
    => use simplex_dualpass
# Arguments
- `a`: the problem as a matrix
"""
function _simplex(a, objective, filter)
    nrow,ncol = size(a) .- 1
    base = collect(ncol-nrow+1:ncol)
    # assume qu'on a l'identité à droite du tableau, donc la base c'est juste
    # les n dernières variables avec n = nb de contraintes (oui on assume aussi
    # que c'est que des contraintes ≤)
    while any(filter(getcosts(a)))
        _,new,old = simplex_pass!(a, objective, filter)
        base[old] = new # maj la base courante
        println("colonne #$new rentre, pivot sur ligne #$old, nouvelle base=$base")
    end
    println("Coût final = $(-a[end])")
    values = a[1:end-1,end]
    solution = zeros(eltype(a), ncol)
    solution[base] = values
    print("Solution = ") ; flush(STDOUT)
    display(solution)
    -a[end], base, values, a
end

getcosts(a) = a[end,1:end-1]
debug(obj::AbstractString) = print(obj)
debug(obj) = (flush(STDOUT); display(obj); println())
debug(objs...) = for obj in objs debug(obj) end
indf(indfunc, arr, f) = indfunc(f(e)?e:typemax(e) for e in arr)
indminpos(itr) = indmin(e>0?e:NaN for e in itr)

simplex_pass(a, objective, filter) = simplex_pass!(copy(a), objective, filter)
# marche avec findmax TODO tester avec findmin
function simplex_pass!(a::Array{<:Real}, heuristic, filter)
    costs = getcosts(a)
    extr, idx = findmax(costs) ; @assert filter(extr) # variable à entrer dans la base
    p = a[1:end-1,end] ./ a[1:end-1,idx]        # pour trouver le pivot
    ip = indf(indmin, p, filter)                # pivot, ⟹ variable à sortir
    a[ip,:] ./= a[ip,idx]                       # normalise la ligne du pivot
    temp = a[:,idx] * a[ip,:]'                  # on fait entrer la nouvelle variable
    temp[ip,:] = 0                              #
    a .-= temp                                  #
    a, idx, ip
end

function simplex_dualpass!(a, f)
    a,idx,ip = simplex_pass!(a', f)
    a',idx,ip
end

@time x,zup,zop = simplex_max(a)
println()
@time simplex_min(a)

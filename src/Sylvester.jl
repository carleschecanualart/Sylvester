module Sylvester

## problems in windows for polymake

using DynamicPolynomials
using EigenvalueSolver
using LinearAlgebra
using Polymake
using Polyhedra
using IterTools
using AlgebraicSolvers

## this procedure provides the Sylvester forms of a given degree in the dense case

function getSylvesterFormsDense(f,x, nu)

    n = length(x)
    
    if (nu < 0) return fill(0*x[1] + 0.0 * im, 0) end

    sylvesterFormList = fill(0*x[1] + 0.0 * im, length(monomials(x, 0:nu)))

    iter = 1

    for mono in monomials(x, 0:nu)
        
        sylvesterForm = fill(0*x[1] + 0.0 * im, n + 1, n + 1)

        ## Here we build the Sylvester form associated with this monomial

        for i = 1:n+1
                        
            F = f[i]      
        
            for j = 1:n
                    
                d = maximum(degree.(terms(F)))
            
                exp = exponents.(mono)[j] 
            
                degreeSylvesterMonomials = d - exp - 1
            
                for sylvesterMono in monomials(x, 0:degreeSylvesterMonomials)
                
                    sylvesterTerm = coefficient(F, x[j]^(exp+1)*sylvesterMono)*sylvesterMono
                
                    sylvesterForm[i,j] += sylvesterTerm
                
                    F -= x[j]^(exp+1)*sylvesterTerm
                
                end
                    
            end
        
            sylvesterForm[i,n+1] = F
        
        end
        
        sylvesterFormList[iter] = det(sylvesterForm) 
    
        iter = iter + 1
    
    end
    
    return sylvesterFormList
    
end

## returns the elimination matrix and the list of monomials of minimal degree

function getResDense(f,x)
    
    n = length(x)
    
    labelPolynomials = [1:(length(f));]
    
    ds = [maximum(degree.(terms(f[i]))) for i in labelPolynomials]
    
    nu = minimum(ds) - 1
    
    alpha = sum([sort(ds)[length(ds) - i] for i = 0:n]) - n - minimum(ds)

    # monomials of that degree alpha
    
    allSubsets = collect(subsets(labelPolynomials, length(x) + 1))

    Σ = monomials(x, 0:alpha)

    numRows = [ length(monomials(x, 0:(alpha-ds[i]))) for i = 1:(length(ds))]

    macaulayMulti = [ monomials(x, 0:(alpha-ds[i])) for i = 1:(length(ds))]
    
    degreeSylv = [sum([ds[i] for i in subs]) - n - 1 - alpha for subs in allSubsets] 
    
    lengthSylv = [0 for subs in allSubsets]
    
    for i in 1:length(allSubsets)
        if (degreeSylv[i] >= 0) lengthSylv[i] = length(monomials(x, 0:(degreeSylv[i]))) end
    end
    
    σ = sum(numRows) + sum(lengthSylv)

    iter = 1;

    res = fill(0.0 + 0.0 * im, σ, length(Σ))

    mapping = Dict(Σ .=> 1:length(Σ))

    for i = 1:length(ds)
        for j = 1:length(macaulayMulti[i])            
            macaulayPoly = macaulayMulti[i][j] * f[i]
            #println(macaulayPoly)
            res[iter,map(k -> mapping[k], monomials(macaulayPoly))] = coefficients(macaulayPoly)

            iter = iter + 1
        end
    end
    
    for subs in allSubsets
        
        fsubs = [f[i] for i in subs]
        
        nup = sum([ds[i] for i in subs]) - n - 1 - alpha
    
        sylvesterFormList = getSylvesterFormsDense(fsubs,x, nup)

        for i = 1:length(sylvesterFormList)            
            sylv = sylvesterFormList[i]
            res[iter, map(k -> mapping[k], monomials(sylv))] = coefficients(sylv)

            iter = iter + 1
        end
        
    end
    
    return res, Σ
    
end

### this procedure provides the Sylvester forms of a given degree in the multihomogeneous case

function getSylvesterFormsMultiDense(f,vargroups, nu)
    
    if (all(nu .>= 0)) 

        sylvesterFormList = fill(0*x[1] + 0.0 * im, length(EigenvalueSolver.getMultiMonomials(vargroups, nu)))

        iter = 1

        for mono in EigenvalueSolver.getMultiMonomials(vargroups, nu)
        
            sylvesterForm = fill(0*x[1] + 0.0 * im, n + 1, n + 1)

        ## Here we build the Sylvester form associated with this monomial

            for i = 1:n+1
                        
                F = f[i]      
        
                for j = 1:n
                    
                    d = maximum(degree.(terms(F)))
            
                    exp = exponents.(mono)[j] 
            
                    degreeSylvesterMonomials = d - exp - 1
            
                    for sylvesterMono in monomials(x, 0:degreeSylvesterMonomials)
                
                        sylvesterTerm = coefficient(F, x[j]^(exp+1)*sylvesterMono)*sylvesterMono
                
                        sylvesterForm[i,j] += sylvesterTerm
                
                        F -= x[j]^(exp+1)*sylvesterTerm
                
                    end
                    
                end
        
                sylvesterForm[i,n+1] = F
        
            end
        
            sylvesterFormList[iter] = det(sylvesterForm) 
    
            iter = iter + 1
    
        end
    
        return sylvesterFormList

    else return fill(0*x[1] + 0.0 * im, 0) end
    
end

### returns matrix and monomials of the multi-homogeneous case

function getResMultiDense(f, vargroups, varsize, ds, nefsubs)
    
    n = length(vargroups)
    
    labelPolynomials = [1:(length(f));]
    
    alpha = sum([ds[i,:]  for i in nefsubs])'  - varsize' - minimum(ds; dims=1)
    
    println(alpha)

    # monomials of that degree alpha
    
    allSubsets = collect(subsets(labelPolynomials, length(x) + 1))

    Σ = EigenvalueSolver.getMultiMonomials(vargroups, alpha)
        
    numRows = [ length(EigenvalueSolver.getMultiMonomials(vargroups, alpha - ds[i,:]')) for i = 1:(size(ds)[1])]
    
    macaulayMulti = [ EigenvalueSolver.getMultiMonomials(vargroups, alpha - ds[i,:]') for i = 1:(size(ds)[1])]
        
    degreeSylv = [sum([ds[i,:] for i in subs])' - varsize' - alpha - [1 1] for subs in allSubsets]
    
    lengthSylv = [0 for subs in allSubsets]
    
    for i in 1:length(allSubsets)
        if (all(degreeSylv[i] .>= 0))
                lengthSylv[i] = length(EigenvalueSolver.getMultiMonomials(vargroups, degreeSylv[i])) 
        end
    end
    
    σ = sum(numRows) + sum(lengthSylv)
    
    println(σ)

    iter = 1;

    res = fill(0.0 + 0.0 * im, σ, length(Σ))

    mapping = Dict(Σ .=> 1:length(Σ))
    
    println(mapping)

    for i = 1:size(ds)[1]
        for j = 1:length(macaulayMulti[i])            
            macaulayPoly = macaulayMulti[i][j] * f[i]
            #println(macaulayPoly)
            res[iter,map(k -> mapping[k], monomials(macaulayPoly))] = coefficients(macaulayPoly)

            iter = iter + 1
        end
    end
    
    for subs in allSubsets
        
        fsubs = [f[i] for i in subs]
        
        nup = sum([ds[i,:] for i in subs])' - varsize' - alpha - [1 1]
    
        sylvesterFormList = getSylvesterFormsMultiDense(fsubs,vargroups, nup)

        for i = 1:length(sylvesterFormList)            
            sylv = sylvesterFormList[i]
            res[iter, map(k -> mapping[k], monomials(sylv))] = coefficients(sylv)

            iter = iter + 1
        end
        
    end
    
    return res, Σ
    
end


function returnPoly(A,b)
    
    P = polyhedron(hrep(A,b))
    
    pts = points(vrep(P))
    npts = npoints(vrep(P))

    poly = fill(0, npts, n)
    
    for j in 1:n
        k = 1
        for p in pts
            poly[k,j] = convert(Int64, floor(p[j]))
            k = k + 1
        end
    end
    return poly
end

### list of sylvester forms in the sparse case -- this checks sigma-positive itself.

function getSylvesterFormsSparse(f, x, A, nu)
    
    n = length(x)
    
    iter = 1
    
    monomialListNu = getAllMonomials(A, nu, x)
    
    println(exponents.(monomialListNu))
    
    nForms = size(monomialListNu)[1]
    
    sylvesterFormList = fill(0.0*x[1] + 0.0, nForms)

    for mono in monomialListNu
        
       sylvesterForm = fill(0.0*x[1] + 0.0, n + 1, n + 1)

        ## Here we build the Sylvester form associated with this monomial

        for i = 1:n+1
                        
            F = f[i]
                    
            for j = 1:n
                                        
                d = maximum(degree.(terms(F)))
            
                exp = exponents.(mono)[j] 
                        
                degreeSylvesterMonomials = d - exp - 1
                            
                for sylvesterMono in monomials(x, 0:degreeSylvesterMonomials)
                                                        
                    sylvesterTerm = coefficient(F, x[j]^(exp+1)*sylvesterMono)*sylvesterMono
                                    
                    sylvesterForm[i,j] += sylvesterTerm
                
                    F -= x[j]^(exp+1)*sylvesterTerm
                
                end
                    
            end
        
            sylvesterForm[i,n+1] = F[1]
        
        end
        
        sylvesterFormList[iter] = det(sylvesterForm) 
    
        iter = iter + 1
    
    end
    
    return sylvesterFormList
    
end

## return monomials in the sparse case

function getAllMonomials(A, b, x)
    poly = returnPoly(A, b)
    f = EigenvalueSolver.getRandomSystem_unmixed(x, poly, [1])[1]
    P = EigenvalueSolver.newtonPolytope(f, x)
    latticePoints = EigenvalueSolver.getLatticePoints(P)
    Σ = [prod.(eachrow(x'.^(latticePoints[i,:]')))[1] for i = 1:size(latticePoints)[1]]
    return Σ
end

## return the matrix in the sparse case

function getResSparse(f, x, A, b, λ, nefsubs)
    
    n = length(x)
    
    picardRank = size(A)[1] - n
    
    labelPolynomials = [1:(length(f));]
    
    canonical = hcat(Matrix(I, n, n) ,A[(n+1):size(A)[1],:])*ones(size(A)[1])
    
    nu = minimum([b[i,(n+1):(n+picardRank)] for i in nefsubs]; dims = 1)[1] - [1 1]'
    
    println(nu)
      
    alpha = Int.(hcat(zeros(Int64, 1, n), sum([b[i,(n+1):(n+picardRank)] for i in nefsubs])' - canonical' - nu'))
    
    # monomials of that degree alpha
    
    allSubsets = collect(subsets(labelPolynomials, n + 1))

    Σ = getAllMonomials(A, alpha[1,:], x)
        
    numRows = [ length(getAllMonomials(A, alpha[1,:] - b[1,:], x)) for i = 1:(size(b)[1])]
    
    macaulayMulti = [ getAllMonomials(A, alpha[1,:] - b[1,:], x) for i = 1:(size(b)[1])]
            
    degreeSylv = [Int.(hcat(zeros(Int64, 1, n), sum([b[i,(n+1):(n+picardRank)] for i in subs])' - canonical')) - alpha for subs in allSubsets]  
    
    lengthSylv = [0 for subs in allSubsets]
    
    for i in 1:length(allSubsets)
        
        if (all(degreeSylv[i] .>= 0))
            
                lengthSylv[i] = length(getAllMonomials(A, degreeSylv[i][1,:], x)) 
        end
    end 
        
    σ = sum(numRows) + sum(lengthSylv)
        
    iter = 1;
    
    mapping = Dict(Σ .=> 1:length(Σ))
    
    res = fill(0.0 + 0.0 * im, σ, length(Σ))

    for i = 1:size(ds)[1]
        
        for j = 1:length(macaulayMulti[i])   
            
            macaulayPoly = macaulayMulti[i][j] * f[i]

            res[iter,map(k -> mapping[k], monomials(macaulayPoly))] = coefficients(macaulayPoly)

            iter = iter + 1
        end
    end
    
    for subs in allSubsets
        
        fsubs = [f[i] for i in subs]
        
        nup = hcat(zeros(Int64, 1, n), minimum([b[i,(n+1):(n+picardRank)] for i in subs]; dims = 1)[1]' - [1 1])
    
        sylvesterFormList = getSylvesterFormsSparse(fsubs, x, A, nup[1,:])

        for i = 1:length(sylvesterFormList)   
            
            sylv = sylvesterFormList[i]

            res[iter, map(k -> mapping[k], monomials(sylv))] = coefficients(sylv)

            iter = iter + 1
        end
    end
    
    return res, Σ
    
end

end



"""


     errormap=DIVAnd_diagapp(P,pmn,len,sv;wheretocalculate=fill(true,size(pmn[1])),Hobs=(),Rmatrix=(),Binv=false,iterB=0)

Calculates an appriximation to the error map exploiting the fact that inv(P) is available via a cholesky decomposition Q'C*C'Q where Q is a permutation matrix. The diagonal component i  is then z'*z where z= inv(Q'C)  ei where ei is normally an array with zeros except in position i where it has a value of 1. Here we exploit that z in that case have only values in specific locations and hence if there is another 1 in ei at a very different location, we can make two calculation at the price of one by summing up the relevant parts of z'z only. And of course you can repeat the reasoning and covering the domain with well placed ones. Advantage compared to AEXERR: works with any number of data points. Should be particularly efficient in higher dimensions and situations where len are small compared to the domain size.
# WARNING: if you provide the wheretocalculate array it WILL BE CHANGED IN PLACE



# Input

* `P` : is the covariance matrix found in the structure s from the analysis s.P`
* `pmn` : are the metrics of the grid
* `len` : the length scales used
* `sv` : the statevector from structure s : s.sv
* `wheretocalculate` : a boolean array of the same dimensions as pmn[1] specifying if in a point the error is to be calculated. Default is everywhere
*                        usefull if you are only interested in part of the error field or when calculating tiles with overlaps
* `Hobs`: s.obsconstrain.H only needed if Binv is true
* `Rmatrix` : epsilon2 only needed if Binv is true
* `Binv` : boolean forcing the calculation of diag(B) by an iterative correction to the estimate of diag(inv(P))
* `iterB` : number of iterations in case B is calculated, the higher epsilon2 the better the iterations converge, so something like 100/epsilon2 ?

# Output

* `errormap` : array of the same dimensions as the analysis including the error variance at the point, relative to a constant background variance.
* `diagB` : array of the same dimensions as the analysis including the spatial structure of the background variance.

"""
function DIVAnd_diagapp(
    P,
    pmn,
    len,
    sv;
    wheretocalculate = fill(true, size(pmn[1])),
    Hobs = (),
    Rmatrix = (),
    Binv = false,
    iterB = 0,
)

    # To change later
    diagRobs = ()
    y = ()
    zy = ()
    diagB = ()
    tutuub = ()

    if Binv
        if Rmatrix == ()
            @warn "No R provided, assuming 1"

            diagRobs = 1.0 .* ones(size(Hobs)[1])
        else
            if isa(Rmatrix, Number)
                diagRobs = Rmatrix .* ones(size(Hobs)[1])
            else
                if ndims(Rmatrix) == 1
                    diagRobs = Rmatrix
                else
                    diagRobs = diag(Rmatrix)
                end

            end
        end
    end



    # Hardwired parameter to control the accuracy. Increase finesse to be closer to "exact" field
    finesse = 1.0
    # Box size for an error calculation and associated half sizes for sums around the center point
    mystrides = zeros(Int32, nfields(pmn))
    halfstrides = zeros(Int32, nfields(pmn))
    # Calculate resolution and grid points needed around a point for the summing
    # A way to relax the following very strong requirement is to use quantiles making for example sure that 90 percent of points are dealt correctly with if L*pmn is strongly variable
    #########################################################
    # Also advection would need a a reassement of the strides
    #########################################################
    # Also maybe include a possible restriction on where the error is to be calculated: eg if done in 3D with moving window one would
    # restrict error calculation to the center part of the domain in z ?
    # Idea optional input boolean matrix wheretocalculate. Start with first point found (and complete the eij) and mask those calculated.
    # proceed until completion of all points: DONE
    # Still to do: match the strides to the desired calculation points if those are already subsampling !!!!
    # That would allow a subsampled parallel calculation. If you throw in the data points in case data coverage is low one should have a nice error calculation !!
    # So take out the stride calculation or rather add a test if provided. Adding the calculation at SOME data points can be done in outer function. To exploit parallel aspect there, one would need to  find a way to deal with data in parallel (quadtree again ?)
    # Proble? What if P is not factorized ? Not feasible anymore
    #########################################################
    ranges = DIVAnd_Lpmnrange(pmn, len)
    for k = 1:nfields(pmn)
        mystrides[k] = Int(round(6.0 * finesse * ranges[k][2]))
        halfstrides[k] = Int(ceil(mystrides[k] / 2))
    end




    # idea: calculate error not with ei but several ei summed up. If points distant enough, error is just the sum AROUND each points
    # if P=C*C'
    # single points z=inv(C)*ei and error in point i = z'*z
    # in 2D needs to be adapted to cover the domain: pack unpack and sum on real domain ?
    # Allocate arrays once
    eij = zeros(Int, size(pmn[1]))
    diagerror = zeros(Float64, size(pmn[1])) .* NaN
    tutuu = zeros(Float64, size(pmn[1]))
    tutu = statevector_pack(sv, (eij,))
    #@show size(tutu)
    z = zeros(Float64, size(P)[1])
    zs = zeros(Float64, size(P)[1])

    if Binv
        y = zeros(Float64, size(P)[1])
        zy = zeros(Float64, size(P)[1])
        diagB = zeros(Float64, size(pmn[1])) .* NaN
        tutuub = zeros(Float64, size(pmn[1]))
    end
    # Get the permutations to apply
    inversep = [findall(x -> x == i, P.factors.p)[1] for i = 1:size(P)[1]]



    mystep = CartesianIndices(tuple(halfstrides...))[end]
    mystride = CartesianIndices(tuple(mystrides...))[end]
    IFI = CartesianIndices(diagerror)[1]
    ILA = CartesianIndices(diagerror)[end]

    # Loop over small box
    # HERE IS THE PLACE WHERE WE DEFINE THE STARTING AND ENDING PLACE IN CASE YOU WANT TO RESTRICT
    uu = findfirst(wheretocalculate)
    while uu != nothing
        I = uu
        #for I in CartesianIndices(tuple(mystrides...))

        eij[:] .= 0
        # and cover the domain with points at the mystride distances
        eij[[I[j]:mystride[j]:ILA[j] for j = 1:ndims(diagerror)]...] .= 1
        #@show size(eij)
        # Go to statevector
        tutu[:] = statevector_pack(sv, (eij,))
        # Get square root part of error
        z[:] = P.factors.PtL \ tutu

        if Binv
            y[:] = z[:]
            for iter = 1:iterB
                # Forces diagonal part of R

                y[:] =
                    z[:] +
                    P.factors.PtL \ (Hobs' * (diagRobs .\ (Hobs * (P.factors.UP \ y[:]))))
            end

        end
        # idea to get inB:
        # y=z
        # iterate
        # y=z+inv(L) H inv(R) Ht inv(U)  y[inversep]
        # then diagnose as for z to get B components from z'y
        # probably still some permutations to check ??
        #
        # Now squared and on the original locations
        zs[:] = z[inversep] .^ 2
        if Binv
            zy[:] = z[inversep] .* y[inversep]
        end
        # Projected back to real space
        tutuu[:], = statevector_unpack(sv, zs)
        if Binv
            tutuub[:], = statevector_unpack(sv, zy)
        end
        # Now each point of eij sums up the contribution aournd its box center
        for IG in findall(x -> x == 1, eij)
            ####################################
            # NEED TO TAKE INTO ACCOUNT MODDIM
            ####################################
            diagerror[IG] = sum(tutuu[max(IFI, IG - mystep):min(IG + mystep, ILA)])
            if Binv
                diagB[IG] = sum(tutuub[max(IFI, IG - mystep):min(IG + mystep, ILA)])
            end

        end
        wheretocalculate[eij.==1] .= false
        uu = findfirst(wheretocalculate)

    end




    if Binv
        return diagerror, diagB
    else
        return diagerror
    end
end

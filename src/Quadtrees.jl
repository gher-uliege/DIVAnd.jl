module Quadtrees

if VERSION >= v"0.7.0-beta.0"
    using Dates
else
    using Compat: @debug
end
using Compat
import Base.length

"""
quadtree (of the higher-dimensional equivalent)
T the type of the coordinates
TA the type of the attributes
N number of dimensions
"""
mutable struct QT{T,TA,N}
    children :: Vector{QT{T,TA,N}}  # vector of child nodes (empty if node is a leaf)
    # list of coordinates (only non-empty if node is a leaf)
    # points[:,i] coordinates of the i-th point
    points :: Array{T,2}
    min :: Vector{T}                # minimum of bounding box
    max :: Vector{T}                # maximim of bounding box
    attribs :: Vector{TA}           # additional attributes (only non-empty if node is a leaf)
end

"""create empty quadtree
"""
QT(TA::DataType,min::Vector{T}, max::Vector{T}) where T =
    QT(QT{T,TA,size(min,1)}[],Matrix{T}(undef,size(min,1),0),min,max,TA[])

"""create a quadtree
"""
QT(points::AbstractArray{T,2},attribs::AbstractVector{TA}) where {T,TA} =
    if VERSION >= v"0.7.0-beta.0"
        QT(QT{T,TA,size(points,1)}[],
           points,
           minimum(points,dims = 2)[:],
           maximum(points,dims = 2)[:],
           attribs)
    else
        QT(QT{T,TA,size(points,1)}[],
           points,
           minimum(points,2)[:],
           maximum(points,2)[:],
           attribs)
    end

function QT(points::AbstractArray{T,2}, min::Vector{T}, max::Vector{T}, attribs::AbstractVector{TA}) where {T,TA}

    if length(attribs) != size(points,2)
        error("QT inconsistent size $(length(attribs)) versus $(size(points,2))")
    end

    return QT(QT{T,TA,size(points,1)}[],points,min,max,attribs)
end

"""
             x1
  +----------+
  |          |
  |   +      |
  |   y      |
  +----------+
 x0

"""
@inline function inside(x0,x1,y)
    insd = true

    @inbounds for i = 1:length(y)
        insd = insd & (x0[i] <= y[i] <= x1[i])
    end
    return insd
end

"""
Test if the n-th bit in a is set. The least significant bit is n = 1.
"""
bitget(a,n) = Bool((a & (1 << (n-1))) >> (n-1))


"""
Test if the rectanges defined by x0,x1 and y0,y1 intersects/overlap
             x1
  +----------+
  |          |
  |   +----------+ y1
  |   |      |   |
  +----------+   |
 x0   |          |
      |          |
      +----------+
     y0
"""
function intersect(x0,x1,y0,y1)
    n = size(x0,1)
#    if (n != length(x1)) || (n != length(y0)) || (n != length(y1))
#        throw(ArgumentError("all arguments of intersect must have the same length"))
#    end

    # https://stackoverflow.com/a/306332/3801401
    cond = true
    for i = 1:length(x0)
        cond = cond && (x0[i] <= y1[i]) && (x1[i] >= y0[i])
    end
    return cond
end






"""
Test if the rectangle defined by x0,x1 is included in rectangle y0,y1
             x1
  +------------+
  |            |
  |   +--+ y1  |
  |   |  |     |
  |   +--+     |
  | y0         |
  +------------+
 x0

"""
@inline function include(x0,x1,y0,y1)
    return inside(x0,x1,y0) && inside(x0,x1,y1)
end



"""
number of points per node
it is always zero for non-leaf nodes
"""
Base.length(qt::QT) = size(qt.points,2)

isleaf(qt) = length(qt.children) == 0

inside(qt::QT,y) = inside(qt.min,qt.max,y)
Base.intersect(qt::QT,y0,y1) = intersect(qt.min,qt.max,y0,y1)
Base.ndims(qt::QT{T,TA,N}) where {T,TA,N} = N

function count(qt::QT)
    if isleaf(qt)
        return length(qt)
    else
        c = 0
        for child in qt.children
            c = c + count(child)
        end

        return c
    end
end

"""
sucess = add!(qt,x,attrib,max_cap = 10)
Add point `x` with the attribute `attrib` to the quadtree `qt`.
`sucess` is true if `x` is within the bounds of the quadtree node `qt` (otherwise
false and the point has not been added)
"""
function add!(qt::QT{T,TA,N},x,attrib,max_cap = 10) where {T,TA,N}
    if length(attrib) != size(x,2)
        @show length(attrib), size(x,2)
        error("QT inconsistent size")
    end

    if !inside(qt,x)
        return false
    else
        if isleaf(qt)
            qt.points = hcat(qt.points,x)

            push!(qt.attribs,attrib)

            if length(qt.attribs) != size(qt.points,2)
                @show length(qt.attribs), size(qt.points,2)
                error("QT inconsistent size")
            end

            # split if necessary
            rsplit!(qt, max_cap)

            return true
        else
            # try to add to all children and returns on first sucessful
            for child in qt.children
                if add!(child,x,attrib,max_cap)
                    return true
                end
            end

            # bounds of child
            cmin = Vector{T}(undef,N)
            cmax = Vector{T}(undef,N)
            xcenter = (qt.max + qt.min)/2

            # create missing child
            @inbounds for i = 1:2^N
                for j = 1:N
                    # all corners of a hypercube
                    if bitget(i-1, j)
                        cmin[j] = qt.min[j]
                        cmax[j] = xcenter[j]
                    else
                        cmin[j] = xcenter[j]
                        cmax[j] = qt.max[j]
                    end
                end

                if all(cmin .< x .<= cmax)
                    child = QT(TA,cmin,cmax)
                    add!(child,x,attrib)
                    push!(qt.children,child)
                    return true
                end
            end

            # should never happen
            error("could not add $(x)")
            return false
        end
    end
end

"""
split a single node
"""
function split!(qt::QT{T,TA,N}) where {T,TA,N}
    # N: number of dimenions

    if isleaf(qt)
        xcenter = (qt.max + qt.min)/2

        nchildren = 2^N
        qt.children = Vector{QT{T,TA,N}}(undef,nchildren)

        # bounds of child
        cmin = Vector{T}(undef,N)
        cmax = Vector{T}(undef,N)
        sel = trues(size(qt.points,2))
        nchildreneff = 0

        @inbounds for i = 1:nchildren
            fill!(sel,true)

            for j = 1:N
                # all corners of a hypercube
                if bitget(i-1, j)
                    sel = sel .& (qt.points[j,:] .<= xcenter[j])
                    cmin[j] = qt.min[j]
                    cmax[j] = xcenter[j]
                else
                    sel = sel .& (qt.points[j,:] .> xcenter[j])
                    cmin[j] = xcenter[j]
                    cmax[j] = qt.max[j]
                end
            end

            child = QT(qt.points[:,sel],copy(cmin),copy(cmax),qt.attribs[sel])
            #=
            points_sel = qt.points[:,sel]
            child = QT(points_sel,
                       minimum(points_sel, dims = 1)[1,:],
                       maximum(points_sel, dims = 1)[1,:],
                       qt.attribs[sel])
            =#
            # add only children with data
            if length(child) > 0
                 nchildreneff = nchildreneff+1
                 qt.children[nchildreneff] = child
            end
        end

        # trim leaves with no data
        resize!(qt.children,nchildreneff)

        # remove points from node
        qt.points = Matrix{T}(undef,N,0)
        qt.attribs = Vector{T}(undef,0)
    end
end

"""
recursive split
"""
function rsplit!(qt::QT{T,TA,N}, max_cap = 10, min_size = zeros(N)) where {T,TA,N}



    if isleaf(qt)
        if length(qt) < max_cap
            # no enougth points, stop recursion
            return
        end

        # check if the minimum size is reached
        min_size_reached = true

        for i = 1:N
            min_size_reached = min_size_reached && ((qt.max[i] - qt.min[i])  < min_size[i])
        end

        if min_size_reached
            # small enought, stop recursion
            return
        end

        # check of all are equal
        allequal = true

        @inbounds for i = 2:size(qt.points,2)
            allequal = allequal & ((@view qt.points[:,i]) == (@view qt.points[:,1]))
        end

        if allequal
            # all points are equal, stop recursion
            return
        end

        split!(qt)
    end

    for child in qt.children
        rsplit!(child,max_cap,min_size)
    end

end


"""
    attribs = within(qt,min,max)

Search all the points within a bounding box defined by the vectors `min` and `max`.
"""
function within(qt::QT{T,TA,N}, min, max) where {T,TA,N}
    nattrib = within_count(qt, min, max)
    attribs = Vector{TA}(undef,nattrib)

    within_buffer!(qt,min,max,attribs)
    return attribs
end


function within_count(qt::QT{T,TA,N}, min, max) where {T,TA,N}
    nattrib = 0

    if !Base.intersect(qt, min, max)
        # nothing to do
        return nattrib
    end

    if isleaf(qt)
        #@show "checking"
        @inbounds for i = 1:length(qt)
            if inside(min,max,@view qt.points[:,i])
                nattrib += 1
            end
        end
        return nattrib
    end

    for child in qt.children
        nattrib += within_count(child, min, max)
    end
    return nattrib
end


function within_buffer!(qt::QT{T,TA,N}, min, max, attribs, nattribs = 0) where {T,TA,N}
    #@show Base.intersect(qt, min, max), min, max, qt.min, qt.max
    if !Base.intersect(qt, min, max)
        # nothing to do
        return nattribs
    end

    @debug "within_buffer $(qt.min) - $(qt.max)"

    if isleaf(qt)
        @debug "leaf $(qt.min) - $(qt.max)"
        # check if node is entirely inside search area min-max
        if include(min,max,qt.min,qt.max)
            # add all
            if nattribs+length(qt) > length(attribs)
                # buffer too small
                return -1
            end

            @inbounds for i = 1:length(qt)
                nattribs += 1
                attribs[nattribs] = qt.attribs[i]
            end
        else
            #@show "check $qt"
            @inbounds for i = 1:length(qt)
                if inside(min,max,@view qt.points[:,i])
                    nattribs += 1
                    if nattribs > length(attribs)
                        # buffer too small
                        return -1
                    end
                    attribs[nattribs] = qt.attribs[i]
                end
            end
        end
        return nattribs
    end

    for child in qt.children
        nattribs = within_buffer!(child, min, max, attribs, nattribs)

        if nattribs == -1
            return -1
        end
    end

    return nattribs
end



function maxdepth(qt::QT{T,TA,N},d = 0) where {T,TA,N}
    if isleaf(qt)
        return d
    end

    dchild = 0
    for child in qt.children
        tmp = maxdepth(child,d)
        if tmp > dchild
            dchild = tmp
        end
    end
    return d + dchild + 1
end


# function qplot(qt::QT)
#     plot([qt.min[1], qt.max[1], qt.max[1], qt.min[1], qt.min[1]],
#          [qt.min[2], qt.min[2], qt.max[2], qt.max[2], qt.min[2]])
# end

# function rplot(qt::QT)
#     qplot(qt)
#     for child in qt.children
#         #@show child
#         rplot(child)
#     end
# end

function Base.show(io::IO,qt::QT; indent = "  ")
    if isleaf(qt)
        printstyled(io, indent,"Leaf $(length(qt))",color=:green)
    else
        printstyled(io, indent,"Node ",color=:blue)
    end
    print(io,"  from $(qt.min) to $(qt.max)\n")

    if !isleaf(qt)
        for child in qt.children
            show(io,child; indent = indent * "  ")
        end
    end
end


# duplicates


function dupset(duplicates)
    d = Vector{Set{Int}}()
    sizehint!(d,length(duplicates) ÷ 2)

    for i = 1:length(duplicates)
        if !(duplicates[i] in d)
            push!(d,duplicates[i])
        end
    end

    return d
end



function catx(x::Tuple)
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(undef,n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end
    return X
end


"""
    dupl = checkduplicates(x,value,delta,deltavalue)

Based on the coordinates `x` (a tuple of longitudes `lons`, latitudes `lats`, depths (`zs`)
and times (`times` vector of `DateTime`)), search for points which are in the same spatio-temporal bounding
 box of length `delta`. `delta` is a vector with 4 elements corresponding to
longitude, latitude, depth and time
(in days). `dupl` a vector of vectors containing the indices of the duplicates.
"""
function checkduplicates(x::Tuple,value,delta,deltavalue;
                         maxcap = 10_000,
                         label = collect(1:size(x[1],1)),
                         factor = 5
                         )
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(undef,n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end

    qt = Quadtrees.QT(X,label)
    Quadtrees.rsplit!(qt, maxcap, delta ./ factor)

    duplicates = Set{Int}[]

    xmin = zeros(n)
    xmax = zeros(n)

    index_buffer = zeros(Int,Nobs)

    @fastmath @inbounds for i = 1:Nobs
        for j = 1:n
            xmin[j] = X[j,i] - delta[j]
            xmax[j] = X[j,i] + delta[j]
        end

        nindex = Quadtrees.within_buffer!(qt,xmin,xmax,index_buffer)

        if nindex > 0
            index = @view index_buffer[1:nindex]

            # check for values
            vv = value[index]
            ii = sortperm(vv)

            istart = 1
            for i=1:length(vv)-1
                if vv[ii[i+1]] - vv[ii[i]] > deltavalue
                    #@show istart:i;

                    if i > istart
                        push!(duplicates,Set(index[ii[istart:i]]))
                    end

                    istart=i+1
                end
            end
            i = length(vv)
            if i > istart
                push!(duplicates,Set(index[ii[istart:i]]))
            end

            #push!(duplicates,Set(index))
            #@show index
        end
    end

    # collect(Set(...)) returns unique elements
    # collect.() transform the list of sets into a list of list
    return sort.(collect.(collect(Set(duplicates))))

end

"""
    dupl = checkduplicates(x1,value1,x2,v2,value2,delta,deltavalue)

Report duplicates of observations in data set (x2,v2) which are also in data set
(x1,v1). `x1` and `x2` are tuples of vectors with the coordinates, `v1` and `v2` are the
corresponding values.
"""
function checkduplicates(x1::Tuple,value1,
                         x2::Tuple,value2,
                         delta, deltavalue;
                         maxcap = 10_000,
                         label1 = collect(1:length(x1[1])),
                         factor = 5
                         )
    X1 = catx(x1)
    X2 = catx(x2)

    n = size(X1,1)
    Nobs1 = size(X1,2)
    Nobs2 = size(X2,2)

    qt = Quadtrees.QT(X1,label1)
    Quadtrees.rsplit!(qt, maxcap, delta ./ factor)
    #@show qt

    duplicates = Vector{Vector{Int}}(undef,Nobs2)

    xmin = zeros(n)
    xmax = zeros(n)

    index_buffer = zeros(Int,Nobs1)

    @fastmath @inbounds for i = 1:Nobs2
        for j = 1:n
            xmin[j] = X2[j,i] - delta[j]
            xmax[j] = X2[j,i] + delta[j]
        end

        nindex = Quadtrees.within_buffer!(qt,xmin,xmax,index_buffer)

        if nindex > 0
            index = @view index_buffer[1:nindex]
            # check for values
            vv = value1[index]
            duplicates[i] = sort(index[abs.(vv .- value2[i]) .< deltavalue])
        else
            duplicates[i] = Int[]
        end
    end

    return duplicates
end


export  QT, rsplit!, add!, show, ndims, count, checkduplicates

end # module

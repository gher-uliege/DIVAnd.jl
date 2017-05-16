# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
using Base.Test

# final grid
#gridsize = (101,101)
gridsize = (21,21,21)

n = length(gridsize)

# observations
xy = ntuple(i -> [0.],n)
f = [2.]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...

mask,pmn,xyi = divand_rectdom([linspace(-1,1,s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
len = ntuple(i -> 0.2,n)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;


@time fi,s = divandrun(mask,pmn,xyi,xy,f,len,epsilon2)
nothing


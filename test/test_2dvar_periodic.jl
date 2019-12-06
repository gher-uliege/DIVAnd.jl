# Testing DIVAnd in 2 dimensions with a periodic dimension

using Test
using DIVAnd

# grid of background field
mask, (pm, pn), (xi, yi) = DIVAnd_rectdom(
    2.0:4.0:358,
    -88.0:4.0:88.0)

# observations
# one at the periodic boundary (360,0) and one at the center (180,0)

x, y = ([360, 180.],[0., 0])
v = [1.,1.]

lenx = 5.;
leny = 5.;

epsilon2 = 0.0001;

#,err,s
va, s = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (x, y),
    v,
    (lenx, leny),
    epsilon2,
    moddim = [360,0],
)

imax = size(va,1)

@test circshift(va,(imax ÷ 2,0)) ≈ va

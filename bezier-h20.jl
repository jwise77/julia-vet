# 1D and 2D Bezier interpolation from Hennicker et al. (2020)
# Gives a better and more thorough description than the other papers.
#
# We only need the 1D and 2D cases for short characteristics because we
# need the quantities on the ray intersections with the cell edges 
# (2D domain) or cell faces (3D domain). No need for linear interpolation
# routines because of the control parameter omega.
#
# Also see bezier-hayek.jl
#
# =========================================================================
#
# Input:
# * f: 1D or 2D field
# * x: actual x-positions in units of cell widths
# * y: (optional) same but for y-positions
# * omega: Bezier control parameter between linear (=1) and quadratic (=0.5)
#
# Assumption: all cell widths are equal
#
# Output: interpolated quantity


#=
Julia example
=============

An example for code written in Julia.
=#

function sphere_vol(r)
    return 4/3 * pi * r^3
end

# %%
# This should work in notebook mode, at least
println(sphere_vol(3))

#%%
# Here's a subsection
# -------------------

phi = (1 + sqrt(5)) / 2

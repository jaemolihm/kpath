using StaticArrays
using Brillouin
using PlotlyJS
using LinearAlgebra
import Spglib, Bravais

a = 1.0
c = 8.0
# Trigonal lattice, space group R-3m, 166
# Example: Bi2Se3 (https://materialsproject.org/materials/mp-541837/)
lattice_nonstandard = SVector(SVector(a*sqrt(3)/2, -a/2, c/3),
                              SVector(0.0, a, c/3),
                              SVector(-a*sqrt(3)/2, -a/2, c/3))
conv_lattice = SVector(SVector(a*sqrt(3), 0, 0),
                       SVector(-a*sqrt(3)/2, a*3/2, 0),
                       SVector(0, 0, c))
sgnum = 166
lattice_standard = SVector(Bravais.primitivize(Bravais.DirectBasis(collect(conv_lattice)), Bravais.centering(sgnum, 3)))

θ = 10 / 180 * pi
rotation = [[cos(θ) sin(θ) 0]; [-sin(θ) cos(θ) 0]; [0 0 1]]
lattice_rotated = Ref(rotation) .* lattice_standard


function irrfbz_path_for_cell(cell)
    # standardize cell
    dset = Spglib.get_dataset(cell)
    sgnum = dset.spacegroup_number
    std_lattice = Bravais.DirectBasis(collect(eachcol(dset.std_lattice)))

    # Calculate kpath for standard primitive cell
    kp = Brillouin.irrfbz_path(sgnum, std_lattice)

    # Convert to original lattice
    # cell.lattice = rotation * primitive_lattice * transformation
    rotation = dset.std_rotation_matrix
    transformation = inv(Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))) * dset.transformation_matrix'

    kp_cart = cartesianize(kp)
    recip_basis = Bravais.reciprocalbasis(Bravais.DirectBasis(collect(eachcol(Matrix(cell.lattice)))))
    for key in keys(kp.points)
        k_cart_rotated = rotation * kp_cart.points[key]
        kp.points[key] = Brillouin.latticize(k_cart_rotated, recip_basis)
    end
    # Need to create a new object because kp.basis is not mutable
    kp_new = KPath(kp.points, kp.paths, recip_basis, kp.setting)

    kp_new, kp
end

cell_to_DirectBasis(cell) = Bravais.DirectBasis(collect(eachcol(Matrix(cell.lattice))))


cell = Spglib.Cell(lattice_nonstandard, [[0, 0, 0]], [0])
kp_new, kp_old = irrfbz_path_for_cell(cell)
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(cell_to_DirectBasis(cell))), kp_old, Layout(title="non-standard cell, old kpath"))
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(cell_to_DirectBasis(cell))), kp_new, Layout(title="non-standard cell, new kpath"))

cell = Spglib.Cell(lattice_rotated, [[0, 0, 0]], [0])
kp_new, kp_old = irrfbz_path_for_cell(cell)
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(cell_to_DirectBasis(cell))), kp_old, Layout(title="rotated cell, old kpath"))
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(cell_to_DirectBasis(cell))), kp_new, Layout(title="rotated cell, new kpath"))

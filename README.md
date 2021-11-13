# kpath_notes
Note on ideas and scripts for k-point path generation for an arbitrary lattice.

## Conventions
* The lattice vectors are represented as a 3\*3 matrix, whose columns are the lattice vectors. For example, lattice[:, 1] is the Cartesian coordinates of the first lattice vector. This is the convention used in `DFTK.jl` and `Spglib.jl`.

### Conventions for standard cell
TODO: The relation between the conventions used in the following.
* spglib
* Bravais.jl
* SeeK-path
* Brillouin.jl
* CDML
* ITA
* SC

### Conventions for primitive cell
There are (at least) three different conventions.
* "CDML primitive" = used in [Bravais.jl](https://thchr.github.io/Crystalline.jl/dev/bravais/#Bravais.primitivebasismatrix) = used in spglib
* "SC (Setyawan and Curtarolo) standard primitive" = Table 2 of HPKOT paper
* "ITA primitive"

TODO: what are these?
* used in `SeeK-path`
* "crystallographic primitive" of HPKOT = Table 3 of HPKOT paper


## Standard conventional cell
Existing kpath algorithms require the cell to be standardized. The standard cell is may contain 1, 2, or 3 primitive cells. SeeK-path uses spglib for the standardization (except for the tricilnic cell: see Issues).

The standardisation consists of two steps: (copied from [SeeK-path](https://github.com/giovannipizzi/seekpath/compare/fix_transformation_matrix))

1. **rearrangement of the lattice vectors as a linear combination of them** (this, e.g., involves permutations, inversion of axes, changing a vector with a linear combination of them, ...). For instance, this will reshuffle the three orthogonal vectors of a tetragonal cell, so that the third vector is the one with different length. Or, for a supercell, it will rescale the cell to get the conventional cell. \[Note by JML: In other words, one multiplies an integer-valued matrix to the right of the lattice matrix.\]

2. **rotation of Cartesian axes**: vectors, their lengths and relative angles are now chosen. The triplet of vectors can now be rotated in Cartesian space, e.g., to have the first vector along the ``x`` axis, the second along ``y`` etc. Note that, in this step, if the cell is refined, cell lengths and relative angles can be slightly adjusted (e.g. the length of the three vectors is set to be the same for almost-cubic systems, and the angles to be exactly 90 degrees even if they weren't so). \[Note by JML: this is called ["idealization"](https://spglib.github.io/spglib/definition.html#def-idealize-cell) in spglib.\]

The first step does not change the system (the point lattice) described by the lattice vectors. The second step can change the system.
But, note that the conventional lattice itself may differ from the primitive one, even before step 2, because of the rescaling of the cell. The conventional system (before step 2) is identical to the original system described by the primitive lattice only when the additional atoms in the cell (e.g. body-centered, face-centered, ...) are considered.

### In `Spglib.jl`,
* The standard lattice can be obtained from `dataset.std_lattice`.
* The transformation matrix for step 1 is `inv(dataset.transformation_matrix)'`. This is an integer-valued matrix.
* The rotation matrix for step 2 is `dset.std_rotation_matrix'`.

Hence, one finds
```julia
dset = Spglib.get_dataset(cell)
@assert cell.lattice ≈ dset.std_rotation_matrix * dset.std_lattice * dset.transformation_matrix' # true
@assert dset.std_lattice ≈ dset.std_rotation_matrix' * cell.lattice * inv(dset.transformation_matrix)' # true
```

For details, see spglib documentation on  [get-symmetry-dataset](https://spglib.github.io/spglib/python-spglib.html#get-symmetry-dataset), [transformation_matrix](https://spglib.github.io/spglib/dataset.html#dataset-origin-shift-and-transformation), and [std_rotation_matrix](https://spglib.github.io/spglib/dataset.html#std-rotation-matrix).

## Standard primitive cell
`SeeK-path` and `Brillouin.jl` gives the k path for the standard primitive cell.
The conversion from the standard conventional cell to the standard primitive cell does not rotate the physical system described by the lattice (i.e. here, there is no step 2 which was needed when obtaining the [standard conventional cell](#standard-conventional-cell)).

The explicit transformation matrix can be obtained using `Bravais.jl`:
```julia
conv_to_prim_matrix = Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))
std_prim_lattice = std_conv_lattice * conv_to_prim_matrix
```

The conversion from an standard conventional lattice to a standard primitive lattice can be done using `Bravais.jl`, or using `Spglib.jl`.

```julia
# 1. Using Bravais
# input lattice (cell) -> standard conventional (std_lattice) -> standard primitive (bravais_prim_lattice)
dset = Spglib.get_dataset(cell)
sgnum = dset.spacegroup_number
conv_lattice = Bravais.DirectBasis(collect(eachcol(dset.std_lattice)))
bravais_prim_lattice = Bravais.primitivize(std_lattice, Bravais.centering(sgnum, 3))

# 1-2. Another way: explicitly multiply the rotation matrix
conv_to_prim_matrix = Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))
@assert hcat(conv_lattice...) * conv_to_prim_matrix ≈ hcat(bravais_prim_lattice...)

# 2. Using Spglib
# input lattice (cell) -> standard primitive (spglib_prim_lattice) at once
spglib_prim_cell = Spglib.standardize_cell(cell, to_primitive=true)
spglib_prim_lattice = Spglib.basis_vectors(spglib_prim_cell)

# Check the two are identical
@assert all(spglib_prim_lattice .≈ bravais_prim_lattice)
```

## Conversion between input lattice and standard primitive lattice
Finally, we can use the two conversions (input cell <-> standard conventional <-> standard primitive) to convert the primitive lattice back to our input lattice.
```julia
# cell (input) <- standard conventional
dset = Spglib.get_dataset(cell)
sgnum = dset.spacegroup_number
@assert cell.lattice ≈ dset.std_rotation_matrix * dset.std_lattice * dset.transformation_matrix'

# standard conventional <- standard primitive
conv_to_prim_matrix = Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))
std_prim_lattice = hcat(Bravais.primitivize(std_lattice, Bravais.centering(sgnum, 3))...)
@assert dset.std_lattice ≈ std_prim_lattice * inv(conv_to_prim_matrix)

# cell (input) <- standard primitive
@assert cell.lattice ≈ dset.std_rotation_matrix * std_prim_lattice * inv(conv_to_prim_matrix) * dset.transformation_matrix'
```

## Issues
- [ ] SeeK-path uses a different standardization for triclinic cell (https://github.com/giovannipizzi/seekpath/tree/fix_transformation_matrix#transformation-matrix)
- [ ] What happens for an un-distorted supercell?
- [ ] How to test? (see [test in `SeeK-path`](https://github.com/giovannipizzi/seekpath/blob/7bb5a3c400dcfd2b1e8f17c636e482f776845ced/seekpath/test_paths_hpkot.py))
- [ ] The conversion between "crystallographic" conventional to primitive cell (Table 3 of HPKOT), and the conversion between "SC standard" conventional to primitive cell (Table 2 of HPKOT) is different (for hR and oA). So, is it okay to just use Spglib.jl for standardization?
- [ ] Is "standard conventional" the same for Spglib, Brillouin.jl, and SeeK-Path?
- [ ] Finish the primitive/standard cell conventions


## Links (references)
* https://github.com/giovannipizzi/seekpath/compare/fix_transformation_matrix
* HPKOT paper: Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Band structure diagram paths based on crystallography*, Comp. Mat. Sci. 128, 140 (2017), https://www.sciencedirect.com/science/article/pii/S0927025616305110?via%3Dihub

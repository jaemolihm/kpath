# kpath_notes
Note on ideas and scripts for k-point path generation for an arbitrary lattice.

## Convention
* The lattice vectors are represented as a 3\*3 matrix, whose columns are the lattice vectors. For example, lattice[:, 1] is the Cartesian coordinates of the first lattice vector. This is the convention used in `DFTK.jl` and `Spglib.jl`.

## Standardization of the cell
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

For details, spglib documentation on  [get-symmetry-dataset](https://spglib.github.io/spglib/python-spglib.html#get-symmetry-dataset), [transformation_matrix](https://spglib.github.io/spglib/dataset.html#dataset-origin-shift-and-transformation), and [std_rotation_matrix](https://spglib.github.io/spglib/dataset.html#std-rotation-matrix).



## Issues
- [ ] SeeK-path uses a different standardization for triclinic cell (https://github.com/giovannipizzi/seekpath/tree/fix_transformation_matrix#transformation-matrix)


## Links (references)
* https://github.com/giovannipizzi/seekpath/compare/fix_transformation_matrix

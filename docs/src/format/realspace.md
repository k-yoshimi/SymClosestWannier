# Real-space data exchange (draft)

**Status: draft for discussion.** Points marked **[open]** are not decided yet. This page is not in the table of contents until it is agreed.

This page proposes how the real-space operators of a Wannier/closest-Wannier (CW) tight-binding model are passed between programs:

- producers: `pw2cw` (symCW), wannier90 (which is planned to contain the CW construction);
- consumers: MultiPie `ModelAnalyzer` (`wannier` and `symcw` modes), wannier-berri, `postcw`.

The conventions themselves (definition of $O_{ab}(\boldsymbol{R})$, Fourier transform, `tb_gauge`) are fixed in [Phase conventions of real-space operators](convention.md).

## Scope of the first version

The first version is deliberately narrow:

- **orthonormal basis only**: the Wannier functions, or the orthonormalized CW orbitals (`Hr`). The non-orthogonal CW basis (`Hr_nonortho`, `Sr`) is outside this version; the position conversions below do not hold for it;
- **operators**: the Hamiltonian $H$ (required) and the position $\boldsymbol{r}$ (optional);
- **lattices**: the producer's cell and the consumer's primitive cell are related by an integer unimodular matrix in the same Cartesian frame (see [Lattice](#lattice-and-reference-cell)).

Spin, $BB$, $CC$, $SR$, $SHR$, $SH$ and symmetrized operators are left to later versions (see [Later versions](#later-versions)).

## Principles

1. **Only real-space data are exchanged.** k-space data on the DFT mesh (`Ek`, `Uk`, `mmn`, …) stay with the producer. The DFT mesh enters only through the lattice vectors of the real-space data and their weights.
2. **wannier90 file layouts are reused** (`_hr.dat`, `_tb.dat`, `_wsvec.dat`, with `seedname.win` and `seedname.nnkp` for the basis), so that the same files serve wannier90, symCW and wannier-berri.
3. **Stored convention: wannier90**, with no orbital positions in the Fourier phase. A consumer that works in the tb gauge converts after reading. The Hamiltonian needs no conversion.
4. **Units and indices as in wannier90:** Å, eV, lattice vectors as rows, integer $\boldsymbol{R}$ in units of the lattice vectors of the file, orbital indices 1-based in files.

## Current support

| | symCW (`pw2cw`) | wannier90 (develop) | MultiPie (develop) | wannier-berri |
| --- | --- | --- | --- | --- |
| `_hr.dat` | writes `seedname_hr.dat.cw` (`write_hr`) | writes | reads `seedname_hr.dat` | reads |
| `_tb.dat` ($H$ and $\boldsymbol{r}$) | writes `seedname_tb.dat.cw` (`write_tb`, needs `mmn`, not with `tb_gauge = true`) | writes | — | reads |
| `_wsvec.dat` | does not use `use_ws_distance`; writes none | writes with `use_ws_distance` | reads if present | does not read it with `_tb.dat`/`_hr.dat` |
| "ndegen applied" form | — | `write_ndegen_applied` | see [Weights](#lattice-vectors-and-weights) | — |

`postcw` reads a Hamiltonian with `hr_input` only if its lattice vectors and `ndegen` are identical to those of the CW model. MultiPie's alternative path (`read_KS`, reading `Ek`, `Uk`, `mmn`, … directly) is not implemented (`NotImplementedError`).

## Lattice vectors and weights

Two representations exist, both in the wannier90 layout:

| | folded (default) | ndegen applied (wannier90 `write_ndegen_applied = true`) |
| --- | --- | --- |
| $\boldsymbol{R}$ list | Wigner–Seitz supercell of `mp_grid` | the list expanded by the `use_ws_distance` images $\boldsymbol{T}$ |
| `ndegen` block | Wigner–Seitz degeneracies $N_{\boldsymbol{R}}$ | all 1 |
| stored value | not divided | divided by $N_{\boldsymbol{R}}$ and by the number of images $n_{\boldsymbol{T}}$ |
| `_wsvec.dat` | applied if `use_ws_distance` | **not** applied (wannier90 still writes it, with an $\boldsymbol{R}$ list that does not match) |

The k-space operator is, for the folded form,

$$
O_{ab}(\boldsymbol{k}) = \sum_{\boldsymbol{R}} \frac{1}{N_{\boldsymbol{R}}} \frac{1}{n_{ab}(\boldsymbol{R})}\sum_{\boldsymbol{T}\in \mathcal{T}_{ab}(\boldsymbol{R})} e^{i\boldsymbol{k}\cdot(\boldsymbol{R}+\boldsymbol{T})} O_{ab}(\boldsymbol{R}),
$$

with $\mathcal{T}_{ab}(\boldsymbol{R}) = \{\boldsymbol{0}\}$ without `use_ws_distance`, and for the "ndegen applied" form a plain sum $\sum_{\boldsymbol{R}} e^{i\boldsymbol{k}\cdot\boldsymbol{R}} O_{ab}(\boldsymbol{R})$.

Rules:

- A file with all `ndegen` equal to 1 may be either form, so the reader cannot tell them apart from the file alone. The representation has to be stated (**[open]**: in a small metadata file, or as consumer input).
- `_wsvec.dat` is applied to the folded form only. MultiPie (develop) applies it whenever the file is present, so it cannot read "ndegen applied" files written with `use_ws_distance`.
- wannier90 writes correct position matrix elements for `transl_inv_full = true` with `use_ws_distance = true` only in the "ndegen applied" form.
- wannier90 `_r.dat` has no `ndegen` block. In the folded form it is not self-contained, so the first version takes $\boldsymbol{r}$ from `_tb.dat`.
- The consumer converts both forms to unit-weight coefficients on reading (expanding the folded form with $N_{\boldsymbol{R}}$ and the images), and all operators of one calculation must use the same image assignment.

Agreed with the MultiPie side: the exchange uses the "ndegen applied" form, one value per bond. MultiPie builds its bonds from symmetry, so a regular $\boldsymbol{R}$ grid is not needed. symCW is to write this form as well.

## Requirements on files

A reader may reject a file that does not satisfy these:

- every $\boldsymbol{R}$ is listed with all `num_wann`² elements;
- $-\boldsymbol{R}$ is listed for every $\boldsymbol{R}$, with $N_{-\boldsymbol{R}} = N_{\boldsymbol{R}}$;
- $H(-\boldsymbol{R}) = H(\boldsymbol{R})^\dagger$ and $r_\alpha(-\boldsymbol{R}) = r_\alpha(\boldsymbol{R})^\dagger$ for each Cartesian component $\alpha$. MultiPie accepts $|O_{ab}(\boldsymbol{R}) - O_{ba}(-\boldsymbol{R})^*| \le 10^{-8} + 10^{-8}|O_{ab}(\boldsymbol{R})|$. The values in a file are rounded: wannier90 writes 6 decimals, symCW 8 significant digits. A file is hermitian to this tolerance only if the producer writes conjugate pairs consistently. In the `_hr.dat.cw` files of the symCW examples and tutorials, the pairs agree exactly;
- in this profile, $\boldsymbol{R}$ and `ndegen` of the folded form are the Wigner–Seitz supercell of `mp_grid` (`seedname.win`). This is the interpolation convention of wannier90. Other $\boldsymbol{R}$ lists that reproduce the same values on the DFT mesh define different values between mesh points.

## Basis information

The wannier90 files do not say which orbital each index is. In the first version this information is taken from the producer's `seedname.win` and `seedname.nnkp`:

- lattice vectors, atoms and `mp_grid` (`.win`);
- for each Wannier function in file order (`.nnkp`, projections): the centre, $l$, $m_r$, $r$, the local $z$ and $x$ axes, and for spinors the spin and the quantization axis.

In the first version, the Wannier basis must equal the model basis up to a permutation and lattice offsets. Bases that need any other transformation (rotation of orbitals or spin, mixing of orbitals) are rejected. Concretely:

- the projection data must be present in `seedname.nnkp`. With `auto_projections` there are none, and the input is rejected;
- local axes $z = (0,0,1)$, $x = (1,0,0)$, and for spinors the quantization axis $(0,0,1)$. MultiPie reads but does not use these axes;
- one radial channel per site and angular momentum. MultiPie drops the radial index $r$;
- each Wannier function is the orbital of its projection. This holds for the CW construction, which is a projection. Maximally localized Wannier functions can mix orbitals, and the matching projection labels alone do not guarantee this; the producer has to assert it.

An explicit correspondence (`ket_wannier` in MultiPie) only selects the permutation. It does not replace the projection data, which are still needed for the origin shift $t$ and the offsets $L_a$.

## Lattice and reference cell

Let $A_{\rm win}$ be the lattice vectors of `seedname.win` and $A_{\rm model}$ the primitive vectors of the model. The first version requires:

- `num_wann` equals the dimension of the model;
- $A_{\rm win} = U A_{\rm ref}$ with an integer unimodular $U$, in the same Cartesian frame, where $A_{\rm ref} \approx A_{\rm model}$ (the model's lattice constants may differ from the DFT ones);
- otherwise the input is rejected (agreed with the MultiPie side). MultiPie (develop) currently warns and continues with $U = 1$.

Cartesian quantities (positions, $\boldsymbol{r}$) refer to the producer's lattice. The consumer therefore uses $A_{\rm ref} = U^{-1} A_{\rm win}$, not the nominal $A_{\rm model}$, for everything that enters the calculation: the conversions below, the Cartesian orbital positions, the reciprocal lattice and the Cartesian k derivatives. MultiPie (develop) currently builds these from $A_{\rm model}$. Using $A_{\rm model}$ with data from $A_{\rm win}$ would deform the model, which this version does not define.

MultiPie accepts deviations of $A_{\rm win} A_{\rm model}^{-1}$ from integers up to 0.1. **[open]** Is this the intended tolerance, or should it be tighter?

The Wannier calculation is done in a primitive cell (agreed with the MultiPie side): MultiPie builds its $H(\boldsymbol{R})$ in the primitive cell, and a primitive cell is cheaper. A conventional DFT cell for a centred lattice (non-unimodular relation) is not supported.

## Mapping onto a MultiPie model

MultiPie (develop) relates the Wannier functions to the model kets with `map_wannier_to_model`:

- `w2m`: model ket of each Wannier function;
- $U$: as above, $x_{\rm model} = x_{\rm win} U + t$ (fractional, row vectors);
- $t$: origin shift;
- $L_a$: lattice offset of Wannier function $a$, $c_a U + t = r_{m} + L_a$ ($c$: projection centre, $r$: model ket position).

Wannier function $a$ in cell $\boldsymbol{X}$ is model ket $m =$ `w2m`$[a]$ in cell $\boldsymbol{X}U + L_a$. The consumer proceeds as follows:

1. **Read and normalize**: unit-weight coefficients, as described in [Weights](#lattice-vectors-and-weights).
2. **Shift the orbital origins** of operators that contain the position. Apply the conversion of [convention.md](convention.md) with the fractional positions $L_a - t$ in the lattice $A_{\rm ref}$, i.e. `to_tb_gauge(ops, irvec, A_ref, atoms_frac=L - t)`. For $\boldsymbol{r}$ this is $r_{aa}(\boldsymbol{0}) \to r_{aa}(\boldsymbol{0}) - (L_a - t)A_{\rm ref}$.
3. **Re-index**: $O'_{mn}(\boldsymbol{R}U + L_b - L_a) = O_{ab}(\boldsymbol{R})$ for all operators (`convert_hr_to_model` does this for $H$). Entries that land on the same key are added, and the union of the new $\boldsymbol{R}$ is zero-filled. The old `ndegen` cannot be carried over, which is why step 1 removes it.
4. **Consumer gauge**: for the tb gauge (MultiPie default `tb_gauge: True`), apply `to_tb_gauge` again with the fractional positions of the model kets, in model ket order, that the consumer uses in its Fourier phase. These need not be the Wannier centres.
5. **Fourier transform** with unit weights.

$H$ is unchanged in steps 2 and 4, and Cartesian components are unchanged under $U$ because the frames are the same.

The algebra of steps 2 and 3 was checked for $U = 1$ with lattice offsets and an origin shift, using the reference model of `tests/test_tb_gauge.py`. The mapped $\boldsymbol{r}$ (and $BB$, $CC$) agree with those of the translated orbitals within the finite-difference accuracy ($10^{-7}$). This checks the formulas only, not an implementation of the pipeline.

An implementation of the pipeline is accepted only when tests of the following pass. They compare $H$ and $\boldsymbol{r}$ at k points off the mesh with an independent reference: the Bloch frame of the producer, transformed to the model basis.

- the actual mapping functions with a nontrivial `w2m`, $U \neq 1$, offsets and an origin shift;
- $A_{\rm ref} \neq A_{\rm model}$ (lattice constants different from the model's);
- folded and "ndegen applied" inputs with nontrivial `ndegen` and `_wsvec.dat`;
- file round trips at file precision;
- the final tb-gauge conversion.

## Later versions

- **Operator sets.** $BB$ and $CC$ are only meaningful together with the Hamiltonian they were built from, $CC$ also needs $BB$, $SR$ needs the spin and $SHR$ needs $SH$. A file set must name its producer calculation, and a symmetrized Hamiltonian (`Hr_sym`) combined with unsymmetrized $BB$/$CC$ must be labelled as an approximation (see [convention.md](convention.md)).
- **Layout for operators without a wannier90 file.** Proposal: the `_hr.dat` header (comment, `num_wann`, `nrpts`, `ndegen`), followed by lines `n1 n2 n3 m n` and the real and imaginary parts of each component, with the components in the order of the operator's indices (for $SR$, $SHR$: spin component first). Operator definitions, units and spin normalization have to be fixed first.
- **Existing symCW files outside the first version**:
  - `_r.dat.cw` and `.s.cw` (spin): `_hr.dat`-style header with an `ndegen` block, unlike wannier90 `_r.dat`;
  - `_sr.dat.cw`: non-orthogonal overlap;
  - `_hr_sym.dat.cw`: a 7-line description block before the `_hr.dat` layout, and SAMB ket order.

  Changing them needs a migration policy: keep the old files, or distinguish the formats explicitly.

## Decisions

Agreed with the MultiPie side (2026-10-10):

1. Operators are passed as real-space files written by symCW or wannier90, and MultiPie maps them onto its basis order. MultiPie does not read `Ek`, `Uk`, `mmn`, … (`read_KS`).
2. The Wannier calculation uses a primitive cell. A lattice that cannot be related to the model's primitive lattice is an error.
3. The exchange uses one value per bond ("ndegen applied" form).
4. symCW writes the symmetrized Hamiltonian (`_hr_sym.dat.cw`) in the plain wannier90 layout, in the basis order of MultiPie.

## Open items

1. The lattice tolerance (MultiPie: 0.1 on the entries of $U$).
2. How the representation (folded or "ndegen applied"), the producer calculation and the basis correspondence are recorded: a small metadata file, or consumer input.
3. MultiPie: do not apply `_wsvec.dat` to "ndegen applied" files; reject lattice mismatch instead of falling back to $U = 1$.
4. Operators needed by MultiPie in later versions, and their layout.

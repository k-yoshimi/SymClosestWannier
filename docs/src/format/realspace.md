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

Spin, $BB$, $CC$, $SR$, $SHR$, $SH$ and symmetrized operators follow in the second version: [Operator sets](#operator-sets), [Matrices needed by each calculation](#matrices-needed-by-each-calculation) and [Files for operators without a wannier90 layout](#files-for-operators-without-a-wannier90-layout).

## Principles

1. **Only real-space data are exchanged.** k-space data on the DFT mesh (`Ek`, `Uk`, `mmn`, …) stay with the producer. The DFT mesh enters only through the lattice vectors of the real-space data and their weights.
2. **wannier90 file layouts are reused** (`_hr.dat`, `_tb.dat`, `_wsvec.dat`, with `seedname.win` and `seedname.nnkp` for the basis), so that the same files serve wannier90, symCW and wannier-berri.
3. **Stored convention: wannier90**, with no orbital positions in the Fourier phase. A consumer that works in the tb gauge converts after reading. The Hamiltonian needs no conversion.
4. **Units and indices as in wannier90:** Å, eV, lattice vectors as rows, integer $\boldsymbol{R}$ in units of the lattice vectors of the file, orbital indices 1-based in files.

## Current support

| | symCW (`pw2cw`) | wannier90 (develop) | MultiPie (develop) | wannier-berri |
| --- | --- | --- | --- | --- |
| `_hr.dat` | writes `seedname_hr.dat.cw` (`write_hr`) | writes | reads `seedname_hr.dat` | reads |
| `_tb.dat` ($H$ and $\boldsymbol{r}$) | writes `seedname_tb.dat.cw` (`write_tb`, not with `tb_gauge = true`; $\boldsymbol{r}$ from `mmn`, or from the projection centres with `tb_position`) | writes | — | reads |
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

- A file with all `ndegen` equal to 1 may be either form, so the reader cannot tell them apart from the file alone. In the exchange the representation is given by the manifest (see [Normalization contract](#normalization-contract)); for legacy files it is given by the consumer.
- `_wsvec.dat` is applied to the folded form only. MultiPie (develop) applies it whenever the file is present, so it cannot read "ndegen applied" files written with `use_ws_distance`.
- wannier90 writes correct position matrix elements for `transl_inv_full = true` with `use_ws_distance = true` only in the "ndegen applied" form.
- wannier90 `_r.dat` has no `ndegen` block. In the folded form it is not self-contained, so the first version takes $\boldsymbol{r}$ from `_tb.dat`.
- The consumer converts both forms to unit-weight coefficients on reading (expanding the folded form with $N_{\boldsymbol{R}}$ and the images), and all operators of one calculation must use the same image assignment.

Agreed with the MultiPie side: the exchange uses the "ndegen applied" form, one value per bond. MultiPie builds its bonds from symmetry, so a regular $\boldsymbol{R}$ grid is not needed. symCW is to write this form as well.

### Normalization contract

To avoid dividing twice (or not at all), the stored value and the side that normalizes are fixed:

- **Stored value in the exchange**: the "ndegen applied" coefficient $\tilde O_{ab}(\boldsymbol{R})$, defined by
  $O_{ab}(\boldsymbol{k}) = \sum_{\boldsymbol{R}} e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\,\tilde O_{ab}(\boldsymbol{R})$ exactly, with every `ndegen` equal to 1. For the folded form, $\tilde O_{ab}(\boldsymbol{R}+\boldsymbol{T}) = O_{ab}(\boldsymbol{R}) / (N_{\boldsymbol{R}}\, n_{ab}(\boldsymbol{R}))$ summed over the images.
- **The writer normalizes, once.** symCW's arrays in memory and in `seedname.hdf5` (`get_oper_R`, `Hr`, …) are folded and not divided (see [convention.md](convention.md)); its writer for the exchange divides them when writing. A reader of the exchange never divides.
- **The image assignment is fixed by the writer.** The $\boldsymbol{R}+\boldsymbol{T}$ list in the file is final: a reader does not apply `_wsvec.dat` or its own minimal-distance replica selection (MDRS) to it. The assignment depends on the centres used and on how ties are decided: on the wannier90 Fe example with a 2×2×2 mesh, image pairs whose distances differ by 1e-5 Å (the tie tolerance) are assigned differently by wannier90 and by wannier-berri, which rounds the centre differences first, and the anomalous Hall conductivity changes from 368 to −1127 S/cm. Writing the assignment explicitly removes this ambiguity.
- **Legacy folded files** (wannier90 default `_hr.dat`/`_tb.dat`, symCW `_hr.dat.cw`) are imported by a separate legacy reader, which normalizes once and then hands over the same unit-weight representation as the exchange reader. The representation of such a file is given by its metadata or by explicit consumer input; a file whose `ndegen` are all 1 and whose representation is given by neither is rejected.
- **Invariant after reading**: every consumer works with unit-weight coefficients and unit `ndegen`. An adapter that divides by `ndegen` (as `symcw_wberri.to_wannierberri_system` does) is called with unit `ndegen`.

## Requirements on files

A reader may reject a file that does not satisfy these:

- every $\boldsymbol{R}$ is listed with all `num_wann`² elements;
- $-\boldsymbol{R}$ is listed for every $\boldsymbol{R}$, with $N_{-\boldsymbol{R}} = N_{\boldsymbol{R}}$;
- $\boldsymbol{R} = 0$ is listed, with zeros if all on-site values vanish (the tb-gauge conversion of symCW, `to_tb_gauge`, needs it);
- $H(-\boldsymbol{R}) = H(\boldsymbol{R})^\dagger$ and $r_\alpha(-\boldsymbol{R}) = r_\alpha(\boldsymbol{R})^\dagger$ for each Cartesian component $\alpha$. MultiPie accepts $|O_{ab}(\boldsymbol{R}) - O_{ba}(-\boldsymbol{R})^*| \le 10^{-8} + 10^{-8}|O_{ab}(\boldsymbol{R})|$. The values in a file are rounded: wannier90 writes 6 decimals, symCW 9 significant digits (`15.8E`). A file is hermitian to this tolerance only if the producer writes conjugate pairs consistently. In the `_hr.dat.cw` files of the symCW examples and tutorials, the pairs agree exactly;
- in this profile, $\boldsymbol{R}$ and `ndegen` of the folded form are the Wigner–Seitz supercell of `mp_grid` (`seedname.win`). This is the interpolation convention of wannier90. Other $\boldsymbol{R}$ lists that reproduce the same values on the DFT mesh define different values between mesh points.

## Basis information

The wannier90 files do not say which orbital each index is. In the first version this information is taken from the producer's `seedname.win` and `seedname.nnkp`:

- lattice vectors, atoms and `mp_grid` (`.win`);
- for each Wannier function in file order (`.nnkp`, projections): the centre, $l$, $m_r$, $r$, the local $z$ and $x$ axes, and for spinors the spin and the quantization axis.

In the first version, the Wannier basis must equal the model basis up to a permutation and lattice offsets. Bases that need any other transformation (rotation of orbitals or spin, mixing of orbitals) are rejected. Concretely:

- the projection data must be present in `seedname.nnkp`. With `auto_projections` there are none, and the input is rejected;
- local axes $z = (0,0,1)$, $x = (1,0,0)$, and for spinors the quantization axis $(0,0,1)$. MultiPie reads but does not use these axes;
- one radial channel per site and angular momentum. MultiPie drops the radial index $r$;
- each Wannier function is the orbital of its projection. This holds for the CW construction, which is a projection. Maximally localized Wannier functions can mix orbitals, and the matching projection labels alone do not guarantee this; the producer has to assert it with `basis.wannier_is_projection` in the [manifest](#manifest-of-an-export). The exchange always has a manifest, in the first version too. Files without one are legacy files (see [Normalization contract](#normalization-contract)): they must be in the projection order (`projection_index` the identity), and the consumer gives the assertion, which includes this order (the files carry no orbital labels, so the order cannot be checked), as explicit input together with the representation. A reader rejects the input if the assertion is absent or false.

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

1. **Read**: unit-weight coefficients, from the exchange reader or, for legacy files, from the legacy reader (see [Normalization contract](#normalization-contract)).
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
- a Hamiltonian with no on-site values (the $\boldsymbol{R} = 0$ block written as zeros);
- legacy files: accepted with the basis assertion given (which includes that the order is the projection order), rejected without it or with a declared order other than the projection order;
- file round trips at file precision;
- the final tb-gauge conversion;
- the manifest checks: a wrong digest, a wrong `export_id`, a missing or false `wannier_is_projection`, a `projection_index` that is not a permutation, a mismatched embedded Hamiltonian of `_tb.dat`, and an export order different from the projection order (`projection_index` not the identity), read back to the same matrices in the model basis.

## Operator sets

Operators other than $H$ are only meaningful together with the calculation they come from. An **export** is the group of files written by one producer run for one consumer; a manifest (below) lists them. All matrices of an export share:

- **one canonical basis order**. The writer permutes every matrix to it before writing: the model ket order of MultiPie when the export is for MultiPie or contains the symmetrized Hamiltonian, otherwise the projection order (`ket_amn`). A set never mixes orders. Because `seedname.nnkp` keeps the projection order, the manifest records for each export index $i$ the index of its projection record, `basis.projection_index[i]` (1-based, identity for the projection order). A reader takes the centre, axes, radial channel and spin of export index $i$ from that record, applies the permutation exactly once (it is not applied to the matrices again), and then relates export indices to model kets with `w2m` and computes the offsets $L_a$ as in [Mapping](#mapping-onto-a-multipie-model);
- the lattice ($A_{\rm win}$), the "ndegen applied" representation and one image assignment;
- the phase convention (wannier90, no orbital positions in the phase), the units and the energy reference.

**Hamiltonian variants.** An export names each Hamiltonian it contains or depends on, with its complete modification record:

- `H_dft`: the Hamiltonian of the Wannier or CW construction (`Hr`);
- `H_sym`: the symmetrized Hamiltonian (`Hr_sym`), if any. In symCW it has the same lattice, $\boldsymbol{R}$ list and `ndegen` as `Hr`, with symmetrized values, in the SAMB ket order; bringing it to the canonical order is a permutation;
- variants with `scissors_shift` or with `shc_bandshift` (first shifted band and energy shift), if used. `first_band` is the 1-based index of the first shifted band among the bands of `seedname.eig` after `exclude_bands`, as `shc_bandshift_firstband` of wannier90; all bands from it upward are shifted. (symCW currently reads `shc_bandshift_firstband` as a float and uses it as a 0-based slice, `get_SHC_R` and `get_response.py`; this has to be fixed before such variants are exported.)

Rules:

1. **Construction Hamiltonian.** $BB$ and $CC$ name the variant they were built from (`built_from`); $SH$ names its variant, and $SHR$ names the same variant and the $SH$ it pairs with. The conversions of [convention.md](convention.md) (tb gauge, origin shift in [Mapping](#mapping-onto-a-multipie-model)) use the construction Hamiltonian of each operator. symCW builds $SH$ and $SHR$ with `shc_bandshift` applied and rejects `scissors_shift` for $BB$ and $CC$. Because the conversions combine operators, an export is closed under them: if $CC$ is present, so is $BB$ with the same `built_from`; if $SR$ is present, so is $SS$; if $SHR$ is present, so is the $SH$ it pairs with, with the same `built_from`; and every variant named in `built_from` is listed in `hamiltonians`. The reader checks this before any mapping or conversion.
2. **Response Hamiltonian.** The manifest names the variant used for the responses (`response_hamiltonian`). If it differs from a construction Hamiltonian (for example `H_sym` with $BB$, $CC$ from `H_dft`), quantities that combine them (orbital magnetization, gyrotropic responses) are approximate, and the consumer reports this.
3. **Embedded Hamiltonians.** The Hamiltonian contained in `_tb.dat` is declared as one of the variants, and the reader checks that it equals that variant (within the number format).
4. **External Hamiltonians** (`hr_input` of `postcw`) belong to an export only if they are written into it as a variant, with a digest; otherwise the consumer rejects the combination.
5. **No mixing.** Files are accepted only through the manifest of their export; a file whose digest does not match is rejected.

### Three kinds of centres

Three sets of positions appear and are recorded separately:

| | meaning | where recorded |
| --- | --- | --- |
| projection centres | positions of the basis orbitals (`seedname.nnkp`), used for the basis correspondence and for the tb gauge of symCW | `seedname.nnkp`, copied into the manifest |
| image centres | centres used by the writer to assign the images $\boldsymbol{T}$ (`use_ws_distance`; for wannier90 the centres of `seedname.chk`) | manifest (`images`) |
| phase centres of a response | positions in the Fourier phase of the consumer's calculation (zero = wannier90 convention, projection centres, Wannier centres, model ket positions) | chosen by the consumer and reported with its results |

Most responses do not depend on the phase centres. Two do, and for them the consumer has to state the choice:

- Qiao's spin Hall conductivity when $SH$ is not $S H$ in the Wannier space (disentangled models). With the same operators, wannier90 and symCW (`tb_gauge = false`) correspond to zero phase centres, and symCW with `tb_gauge = true` to the projection centres (checked against wannier-berri to 1e-8).
- The shift current with a finite `sc_eta` without the correction of Lihm, PRB 103, 247101 (2021).

For point-orbital models (no `mmn`) the phase centres are the orbital positions by construction (see below).

### Manifest of an export

**[open]** Proposal: a JSON file `seedname_operators.json`, written by the producer next to the files. It is mandatory for the exchange; `seedname.win` and `seedname.nnkp` remain the authoritative basis description (local axes, radial channels, spin axes), are listed with digests, and a disagreement between them and the manifest is an error.

```json
{
  "format": "realspace-operators",
  "format_version": 1,
  "export_id": "<random identifier written into every file header>",
  "producer": {"code": "SymClosestWannier", "version": "2.0.7", "seedname": "Fe",
               "source_digest": "<sha256 of seedname.hdf5>"},
  "representation": "ndegen_applied",
  "convention": "wannier90",
  "units": {"length": "angstrom", "energy": "eV"},
  "lattice": {"unit_cell_cart": [[1.435, 1.435, 1.435], [-1.435, 1.435, 1.435], [-1.435, -1.435, 1.435]],
              "mp_grid": [2, 2, 2]},
  "images": {"use_ws_distance": false, "centres_frac": null, "tolerance": null,
             "assignment_digest": "<sha256 of the sorted list of (R, m, n) with nonzero weight>"},
  "basis": {"num_wann": 18, "order": "ket_amn", "kets": ["..."], "projection_index": [1, 2, "..."],
            "wannier_is_projection": true, "projection_centres_frac": [["..."]],
            "spinors": true, "files": {"Fe.win": "<sha256>", "Fe.nnkp": "<sha256>"}},
  "energy": {"fermi_energy": 12.6279},
  "hamiltonians": {
    "H_dft": {"file": "Fe_hr.dat", "sha256": "<...>", "modifications": {}},
    "H_shc": {"file": "Fe_hr_shc.dat", "sha256": "<...>",
              "modifications": {"shc_bandshift": {"first_band": 5, "energy_shift": 0.3}}}
  },
  "response_hamiltonian": "H_dft",
  "operators": {
    "AA":  {"file": "Fe_tb.dat", "sha256": "<...>", "embedded_hamiltonian": "H_dft", "source": "mmn"},
    "BB":  {"file": "Fe_bb.dat", "sha256": "<...>", "built_from": "H_dft"},
    "CC":  {"file": "Fe_cc.dat", "sha256": "<...>", "built_from": "H_dft"},
    "SS":  {"file": "Fe_spin.dat", "sha256": "<...>"},
    "SH":  {"file": "Fe_spin_h.dat", "sha256": "<...>", "built_from": "H_shc"},
    "SR":  {"file": "Fe_spin_r.dat", "sha256": "<...>"},
    "SHR": {"file": "Fe_spin_hr.dat", "sha256": "<...>", "built_from": "H_shc", "pairs_with": "SH"}
  }
}
```

- `format_version` is the version of this schema. Which operators an export contains is given by `operators`; a reader rejects an unknown `format_version`. A listed operator that the requested calculation does not need may be unsupported by the reader and is then not used, but its file still passes the structural checks below (digest, `export_id`, $\boldsymbol{R}$ list); a file whose layout the reader cannot parse for these checks is an error; an operator it needs that is missing or unsupported is an error. All mandatory checks below apply in either case.
- `operators.AA.source`: `mmn` (finite differences of the overlaps) or `projection_centres` (`tb_position`, a point-orbital approximation). With `projection_centres` the consumer computes with internal terms only (see [Matrices](#matrices-needed-by-each-calculation)), with the projection centres, mapped to the consumer's cell representatives, as phase centres.
- `basis.projection_index` is a permutation of 1, ..., `num_wann`; a reader rejects anything else. `basis.wannier_is_projection` is `true` only if export function $i$ is the orbital of projection record `projection_index[i]` (see [Basis information](#basis-information)).
- `export_id` is written in the header comment of every operator and Hamiltonian file of the export. `seedname.win` and `seedname.nnkp` are checked by their digests only.
- `images.assignment_digest` is an opaque provenance token of the producer, not a description from which the assignment can be recomputed.

Mandatory checks of a reader: `format_version`; digests of all listed files; the closure of the operator set under the conversions (rule 1 of [Operator sets](#operator-sets)); `projection_index` and `wannier_is_projection`; the `export_id` in each operator and Hamiltonian file; lattice and basis against `seedname.win` and `seedname.nnkp` through `projection_index`; the embedded Hamiltonian of `_tb.dat`; and that every file has the same $\boldsymbol{R}$ list. The common $\boldsymbol{R}$ list only proves a common grid. That the per-pair image assignment is the same is guaranteed by provenance: all files are written by one writer pass with one assignment.

## Matrices needed by each calculation

The table lists the real-space matrices that the wannier-berri (26.7) calculators need, from `system/needed_data.py` and the formulas they use. "Internal terms" is wannier-berri's `force_internal_terms_only`, which keeps only $H$ and $SS$: it gives the tight-binding (point-orbital) approximation, in which the Wannier functions are taken as point orbitals at the centres of the Fourier phase.

| quantity | wannier-berri calculators | matrices | internal terms only |
| --- | --- | --- | --- |
| density of states, Fermi surface, band velocity and inverse mass | `DOS`, `CumDOS`, `JDOS`, tabulators `Energy`, `Velocity`, `InvMass` | $H$ | $H$ |
| Ohmic, classical Hall and nonlinear Drude conductivities | `Ohmic_*`, `Hall_classic_*`, `NLDrude_FermiSea`, `_FermiSurf`, `_Fermider2` | $H$ | $H$ |
| anomalous Hall conductivity, Berry curvature (dipole), eMChA | `AHC`, `BerryDipole_*`, `OmegaOmega`, `eMChA_FermiSurf`, `NLDrude_Zeeman_orb_Omega`, tabulator `BerryCurvature` | $H$, $AA$ | $H$ |
| optical conductivity, shift and injection currents | `OpticalConductivity`, `ShiftCurrent`, `InjectionCurrent` | $H$, $AA$ | $H$ |
| spin, spin gyrotropic magnetic effect, Zeeman (spin) Drude | `Spin`, `GME_spin_*`, `NLDrude_Zeeman_spin`, tabulator `Spin` | $H$, $SS$ | $H$, $SS$ |
| Zeeman (spin) anomalous Hall | `AHC_Zeeman_spin` | $H$, $AA$, $SS$ | $H$, $SS$ |
| orbital magnetization, orbital GME, Zeeman (orbital) | `Morb`, `GME_orb_*`, `AHC_Zeeman_orb`, `NLDrude_Zeeman_orb`, tabulator `OrbitalMoment` | $H$, $AA$, $BB$, $CC$ | $H$ |
| spin Hall conductivity, Qiao | static `SHC` and tabulator `SpinBerry` with `kwargs_formula={"spin_current_type": "qiao"}`; dynamic `SHC` with `SHC_type="qiao"` | $H$, $AA$, $SS$, $SR$, $SH$, $SHR$ | not available |
| spin Hall conductivity, Ryoo | same, `"ryoo"` | $H$, $AA$, $SS$, $SA$, $SHA$ (`NeededData` also loads $SH$) | not available |
| spin Hall conductivity, simple | same, `"simple"` | $H$, $AA$, $SS$ | $H$, $SS$ |
| quantum metric | `QuantumMetric_FermiSea`, `QuantumMetric_Vel_DQ` | $H$, $AA$, $FF$ | $H$ |
| orbital-moment-related (OSD; loading requirements of `NeededData(OSD=True)`, not checked per formula) | — | $H$, $AA$, $BB$, $CC$, $OO$, $GG$ (normally replaced by $FF$) | $H$ |

What symCW provides:

| matrix | symCW | file today |
| --- | --- | --- |
| $H$ | yes | `_hr.dat.cw`, `_tb.dat.cw` |
| $AA$ | from `seedname.mmn`, or the projection centres (`tb_position`) | `_tb.dat.cw`, `_r.dat.cw` |
| $BB$ | from `mmn` (and the eigenvalues) | none |
| $CC$ | from `mmn` and `uHu` | none |
| $SS$ | from `seedname.spn` | `.s.cw` (own layout) |
| $SR$, $SH$, $SHR$ | from `mmn` and `spn` | none |
| $SA$, $SHA$, $OO$, $GG$, $FF$ | not implemented (need `siu`, `shu`, `uIu`) | — |

$SA$, $SHA$, $OO$, $GG$ and $FF$ are outside this exchange. Calculations that need them use wannier90's files with wannier-berri's own reader (which assigns the images itself), not real-space files from symCW; they are not planned for symCW.

wannier90's defaults broaden the energy denominators of the spin Hall conductivity (`kubo_adpt_smr = true`) and drop states above `kubo_eigval_max`; wannier-berri does neither. Comparisons set these explicitly.

For point-orbital models (no `mmn`), the consumer uses the internal terms with the projection centres as phase centres, instead of synthetic $BB$, $CC$ built from the centres. The spin Hall conductivity is then available only in the simple form.

## Files for operators without a wannier90 layout

**[open]** Proposal for $BB$, $CC$, $SS$, $SH$, $SR$, $SHR$. (The symmetrized Hamiltonian in MultiPie order uses the plain `_hr.dat` layout under a new name, see [Compatibility](#compatibility-with-existing-symcw-files).)

| operator | definition (wannier90 convention) | components | unit | file |
| --- | --- | --- | --- | --- |
| $BB_a$ | $\langle\phi_m(\boldsymbol{0})|H(\boldsymbol{r}-\boldsymbol{R})_a|\phi_n(\boldsymbol{R})\rangle$ | 3: $a = x, y, z$ | eV Å | `seedname_bb.dat` |
| $CC_{ab}$ | $\langle\phi_m(\boldsymbol{0})|r_a H (\boldsymbol{r}-\boldsymbol{R})_b|\phi_n(\boldsymbol{R})\rangle$ | 9: $(a, b)$, $b$ fastest | eV Å² | `seedname_cc.dat` |
| $SS_s$ | $\langle\phi_m(\boldsymbol{0})|\sigma_s|\phi_n(\boldsymbol{R})\rangle$, Pauli matrices (eigenvalues ±1) | 3: $s = x, y, z$ | 1 | `seedname_spin.dat` |
| $SH_s$ | $\langle\phi_m(\boldsymbol{0})|\sigma_s H|\phi_n(\boldsymbol{R})\rangle$ | 3 | eV | `seedname_spin_h.dat` |
| $SR_{sa}$ | $\langle\phi_m(\boldsymbol{0})|\sigma_s(\boldsymbol{r}-\boldsymbol{R})_a|\phi_n(\boldsymbol{R})\rangle$ | 9: $(s, a)$, $a$ fastest | Å | `seedname_spin_r.dat` |
| $SHR_{sa}$ | $\langle\phi_m(\boldsymbol{0})|\sigma_s H(\boldsymbol{r}-\boldsymbol{R})_a|\phi_n(\boldsymbol{R})\rangle$ | 9 | eV Å | `seedname_spin_hr.dat` |

- $\boldsymbol{r}$ in these definitions is in the same finite-difference sense as $AA$: symCW, wannier90 (`get_SHC_R`) and wannier-berri compute $SR$, $SHR$ as $i\sum_{\boldsymbol{b}} w_b\,\boldsymbol{b}\,(\dots)$, the factor $i$ being that of $A = i\langle u|\nabla_{\boldsymbol{k}} u\rangle$. Independently constructed finite differences need not agree exactly between producers (wannier-berri, for example, puts the centre phase into the finite difference of $AA$); a fixed image assignment removes only the interpolation ambiguity.
- The component order $(s, a)$ is that of symCW and wannier90. Conversions are done by the consumer's adapter, on the full tensors: the origin shift and the tb gauge of [convention.md](convention.md) on full $CC$, and only then, for wannier-berri, the transposition of $SR$, $SHR$ to its $(a, s)$ order and the pseudovector $i(CC_{ab} - CC_{ba})$ (as `symcw_wberri.to_wannierberri_system` does).
- Hermiticity: $SS(-\boldsymbol{R}) = SS(\boldsymbol{R})^\dagger$ and $CC_{ab}(-\boldsymbol{R}) = CC_{ba}(\boldsymbol{R})^\dagger$. $BB$, $SH$, $SR$, $SHR$ have no such relation, and neither a writer nor a reader symmetrizes them; $CC$ is not antisymmetrized either.
- Layout: the `_hr.dat` header (comment line with the operator name, its definition and the `export_id`, `num_wann`, `nrpts`, the `ndegen` block), then one line per $(\boldsymbol{R}, m, n)$ in the order of `_hr.dat` ($m$ fastest, then $n$, then $\boldsymbol{R}$): `n1 n2 n3 m n` followed by the real and imaginary parts of each component in the order of the table.
- Number format **[open]**: at least 10 significant digits (wannier90 writes 6 decimals, symCW 9 significant digits). For the relations above to hold at the reader's tolerance, the writer first checks them on the unrounded arrays (an error there is a defect of the producer, not hidden by the export), then writes the conjugate partners from one array and checks the written file by reading it back.
- Size: after the image expansion the file lists all `num_wann`² pairs for every $\boldsymbol{R}$ in the union of the images, most of them zero; with 9 components this grows quickly. **[open]** A sparse layout (only nonzero pairs) may be needed; a benchmark on a realistic model (file size and the peak memory of reading and mapping, where the orbital offsets enlarge the $\boldsymbol{R}$ list further) is to be done before freezing the layout.
- The file names avoid symCW's `_sr.dat.cw` (the overlap matrix); names differing only in case are avoided because file systems may ignore case.

## Compatibility with existing symCW files

Existing files keep their names and layouts; new layouts get new names:

- `_r.dat.cw` and `.s.cw` (spin): `_hr.dat`-style header with an `ndegen` block, unlike wannier90 `_r.dat`. Unchanged.
- `_sr.dat.cw`: non-orthogonal overlap. Unchanged.
- `_hr_sym.dat.cw`: a 7-line description block before the `_hr.dat` layout, folded, SAMB ket order. Unchanged. The symmetrized Hamiltonian of decision 4 (plain wannier90 layout, MultiPie order, "ndegen applied") is written under a new name **[open]**, proposed `seedname_hr_sym.dat` with the manifest, so that existing readers of `_hr_sym.dat.cw` keep working. Whether the old files are dropped later is to be decided with the users.

## Acceptance tests

In addition to those of [Mapping](#mapping-onto-a-multipie-model), an implementation of the second version is accepted when tests of the following pass:

- write–read round trips of every operator at file precision, with nontrivial `ndegen`, image assignments and centres, and the component order of $CC$, $SR$, $SHR$ checked with unequal components;
- conversion (tb gauge, origin shift) with a construction Hamiltonian different from the response Hamiltonian;
- rejection of an operator set not closed under the conversions: $CC$ without $BB$, $BB$ and $CC$ built from different variants, $SR$ without $SS$, $SHR$ without its $SH$, a `built_from` variant not in `hamiltonians`;
- `shc_bandshift` variants: with excluded bands, `first_band` = 1 and an interior `first_band`, the shifted eigenvalues are exactly those of bands `first_band`, …, `num_bands` of `seedname.eig` after `exclude_bands`;
- rejection of: a file not in the manifest or with a wrong digest, a stale manifest, a mismatched embedded Hamiltonian, mixed exports, and a reader that would apply `_wsvec.dat` to an "ndegen applied" export;
- responses from the read files equal to those of the producer's own operators (as `symcw-wberri` does for `postcw`), including the spin Hall conductivity with stated phase centres;
- internal-terms-only calculations for a point-orbital model, and rejection of Qiao and Ryoo spin Hall conductivities in that mode.

## Decisions

Agreed with the MultiPie side (2026-10-10):

1. Operators are passed as real-space files written by symCW or wannier90, and MultiPie maps them onto its basis order. MultiPie does not read `Ek`, `Uk`, `mmn`, … (`read_KS`).
2. The Wannier calculation uses a primitive cell. A lattice that cannot be related to the model's primitive lattice is an error.
3. The exchange uses one value per bond ("ndegen applied" form).
4. symCW writes the symmetrized Hamiltonian in the plain wannier90 layout, in the basis order of MultiPie (file name: see [Compatibility](#compatibility-with-existing-symcw-files)).

## Open items

1. The lattice tolerance (MultiPie: 0.1 on the entries of $U$).
2. The manifest (`seedname_operators.json`): name, fields and checks.
3. MultiPie: do not apply `_wsvec.dat` to "ndegen applied" files; reject lattice mismatch instead of falling back to $U = 1$.
4. File names (including the new name of the plain symmetrized Hamiltonian), layout (dense or sparse) and number format of the operators without a wannier90 layout.
5. Which of these operators MultiPie needs, and when.

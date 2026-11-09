"""
generate_permeability.py
========================
Reads an existing OpenFOAM volScalarField file (constant/K) and treats each
cell's current value as the *local* geometric mean of a log-normal
distribution.  A single random perturbation is drawn for every cell, centred
on that cell's own mean, producing a spatially heterogeneous field that honours
the original spatial structure.

Cells where K = 0 are treated as impermeable (e.g. solid inclusions) and are
written back unchanged as 0.

Supported existing internalField formats
-----------------------------------------
    internalField   uniform 1e-05;          ← same mean broadcast to all cells
    internalField   nonuniform List<scalar> ← per-cell local mean
                    4000
                    ( ... )
                    ;

If the target file does not yet exist, FALLBACK_MEAN_K is broadcast to every
cell as the uniform local mean and a complete new file is written.

Usage
-----
    python generate_permeability.py

Outputs
-------
    constant/K          – existing file with internalField replaced in-place
                          (or a new file if constant/K does not exist)
    K_field.npy         – raw NumPy array (for inspection / post-processing)
    K_histogram.png     – histogram of non-zero generated values (log-scale)

Dependencies
------------
    numpy, matplotlib
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend (safe for servers)
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. USER-CONFIGURABLE PARAMETERS
# ---------------------------------------------------------------------------

# Mesh dimensions — must match the blockMeshDict / polyMesh cell count
NX, NY, NZ = 300, 200, 1          # cells in x, y, z

# Log-normal spread parameter (std-dev of ln K, dimensionless)
# Typical vadose-zone values: sandy loam ≈ 1.5, silt loam ≈ 1.2
SIG_LN = 1.5

# Fallback geometric mean used ONLY when constant/K does not yet exist
# (K in m/s — e.g. 1e-5 for sandy loam, 1e-6 for silt loam)
FALLBACK_MEAN_K = 1e-5

# Minimum permissible K for active (non-zero) cells. Any sampled value below
# this threshold is clamped up to K_MIN. Does NOT affect cells that are
# already 0 — those remain untouched (see "K = 0" handling in section 2).
K_MIN = 9.72e-7   # m/s

# Random seed (set to None for a different field every run)
SEED = 42

# Paths
OF_OUTPUT_DIR = "constant"
OF_FIELD_NAME = "K"
NPY_OUTPUT    = "K_field.npy"
HIST_OUTPUT   = "K_histogram.png"

# ---------------------------------------------------------------------------
# 1b. READ LOCAL MEANS FROM EXISTING FILE
# ---------------------------------------------------------------------------

of_path  = os.path.join(OF_OUTPUT_DIR, OF_FIELD_NAME)
original = None   # will hold the file text if the file exists

def _parse_values_from_field(text: str, nx: int, ny: int, nz: int) -> np.ndarray | None:
    """
    Extract per-cell K values from an OpenFOAM internalField entry.

    • uniform <value>  → broadcasts the scalar to every cell, returning a
                         flat array of length nx*ny*nz.
    • nonuniform list  → returns the parsed value array as-is.

    Returns None if no internalField block is found.
    """
    n_cells = nx * ny * nz

    # --- uniform scalar: broadcast to all cells ---
    uniform_match = re.search(
        r"internalField\s+uniform\s+([0-9eE+\-\.]+)\s*;",
        text,
    )
    if uniform_match:
        scalar = float(uniform_match.group(1))
        print(f"  internalField uniform {scalar:.3e}  →  broadcasting to {n_cells} cells")
        return np.full(n_cells, scalar)

    # --- nonuniform list: parse every value ---
    nonuniform_match = re.search(
        r"internalField\s+nonuniform\s+List<\w+>\s*"
        r"(\d+)\s*\(\s*(.*?)\s*\)\s*;",
        text,
        re.DOTALL,
    )
    if nonuniform_match:
        declared_n = int(nonuniform_match.group(1))
        values = np.fromstring(nonuniform_match.group(2), sep="\n")
        print(f"  internalField nonuniform  →  parsed {values.size} values "
              f"(declared {declared_n})")
        if values.size != n_cells:
            raise ValueError(
                f"Cell count mismatch: file has {values.size} values but "
                f"NX×NY×NZ = {nx}×{ny}×{nz} = {n_cells}. "
                f"Update NX/NY/NZ to match your mesh."
            )
        return values

    return None

if os.path.exists(of_path):
    with open(of_path, "r") as fh:
        original = fh.read()

    local_means = _parse_values_from_field(original, NX, NY, NZ)
    if local_means is not None:
        n_nonzero = int(np.count_nonzero(local_means))
        n_zero    = len(local_means) - n_nonzero
        print(f"Read existing constant/K  →  {n_nonzero} active cells, "
              f"{n_zero} zero/impermeable cells")
    else:
        local_means = np.full(NX * NY * NZ, FALLBACK_MEAN_K)
        print(f"Warning: constant/K found but no internalField parsed. "
              f"Falling back to uniform mean = {FALLBACK_MEAN_K:.3e} m/s")
else:
    local_means = np.full(NX * NY * NZ, FALLBACK_MEAN_K)
    print(f"constant/K not found — using fallback mean K = {FALLBACK_MEAN_K:.3e} m/s")

# ---------------------------------------------------------------------------
# 2. GENERATE THE LOG-NORMAL FIELD (cell-local means)
# ---------------------------------------------------------------------------

rng     = np.random.default_rng(SEED)
n_cells = NX * NY * NZ

# Identify zero cells (impermeable / solid inclusions) before sampling
zero_mask    = (local_means == 0.0)         # True where K should stay 0
active_mask  = ~zero_mask                   # True where we sample
n_zero       = int(zero_mask.sum())
n_active     = int(active_mask.sum())

# Allocate output array; zero cells remain 0 throughout
K_flat = np.zeros(n_cells)

if n_active > 0:
    # Per-cell log-space means: mu_ln_i = ln(K_local_i)
    # SIG_LN is a global spread applied uniformly to every active cell
    mu_ln_local = np.log(local_means[active_mask])   # shape: (n_active,)
    ln_K_active = rng.normal(loc=mu_ln_local, scale=SIG_LN)
    sampled     = np.exp(ln_K_active)

    # Enforce minimum threshold: any sampled value below K_MIN is raised
    # to K_MIN. Cells held at 0 (zero_mask) are unaffected — handled separately.
    n_clamped = int(np.count_nonzero(sampled < K_MIN))
    if n_clamped > 0:
        print(f"  Clamped {n_clamped} cell(s) below K_MIN = {K_MIN:.3e}  →  raised to K_MIN")
    sampled = np.maximum(sampled, K_MIN)

    K_flat[active_mask] = sampled

if n_zero > 0:
    print(f"  Skipped {n_zero} zero/impermeable cells — written back as 0.0")

# Reshape for .npy save (preserves spatial ordering)
K = K_flat.reshape((NX, NY, NZ), order="C")

# Reported statistics (active cells only)
active_vals = K_flat[active_mask]
if n_active > 0:
    geo_mean   = float(np.exp(np.mean(np.log(active_vals))))
    arith_mean = float(active_vals.mean())
    print(f"Generated {n_cells} cells  ({NX}×{NY}×{NZ})  —  {n_active} active")
    print(f"  Geometric mean K  (active): {geo_mean:.3e}  m/s")
    print(f"  Arithmetic mean K (active): {arith_mean:.3e}  m/s")
    print(f"  Min K (active)            : {active_vals.min():.3e}  m/s")
    print(f"  Max K (active)            : {active_vals.max():.3e}  m/s")
else:
    geo_mean = arith_mean = 0.0
    print("Warning: all cells are zero — no active cells to sample.")

# Save the raw array
np.save(NPY_OUTPUT, K)
print(f"Saved raw field → {NPY_OUTPUT}")

# ---------------------------------------------------------------------------
# 3. DIAGNOSTIC HISTOGRAM
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7, 4))
if n_active > 0:
    ax.hist(active_vals, bins=60, color="#1D9E75", edgecolor="white", linewidth=0.3)
    ax.axvline(geo_mean,   color="#D85A30", linewidth=1.5,
               label=f"Geometric mean = {geo_mean:.2e}")
    ax.axvline(arith_mean, color="#378ADD", linewidth=1.5, linestyle="--",
               label=f"Arithmetic mean = {arith_mean:.2e}")
    ax.legend(fontsize=9)
else:
    ax.text(0.5, 0.5, "No active cells", transform=ax.transAxes,
            ha="center", va="center", fontsize=12)
ax.set_xscale("log")
ax.set_xlabel("Permeability K  (m/s)")
ax.set_ylabel("Cell count")
ax.set_title(f"Log-normal K distribution  ({n_active} active, {n_zero} zero cells)")
fig.tight_layout()
fig.savefig(HIST_OUTPUT, dpi=150)
plt.close(fig)
print(f"Saved histogram    → {HIST_OUTPUT}")

# ---------------------------------------------------------------------------
# 4. PATCH internalField IN AN EXISTING OpenFOAM K FILE
# ---------------------------------------------------------------------------

os.makedirs(OF_OUTPUT_DIR, exist_ok=True)
# of_path and original are already set in section 1b

# Build the replacement internalField block
values_block = "\n".join(f"    {v:.6e}" for v in K_flat)
new_internal = (
    f"internalField   nonuniform List<scalar>\n"
    f"{n_cells}\n"
    f"(\n"
    f"{values_block}\n"
    f")\n"
    f";"
)

# Metadata comment to inject / replace just above internalField
meta_comment = (
    f"// Heterogeneous hydraulic conductivity (log-normal, cell-local means)\n"
    f"//   sig_ln = {SIG_LN:.4f}  (spread, std-dev of ln K)\n"
    f"//   K_min  = {K_MIN:.3e}  (floor applied to active cells)\n"
    f"//   seed   = {SEED}\n"
    f"//   cells  = {n_cells}  ({NX} x {NY} x {NZ})  —  {n_active} active, {n_zero} zero\n"
    f"//"
)

if original is not None:
    # ------------------------------------------------------------------
    # PATCH MODE: replace the existing internalField block in-place.
    # Two dedicated patterns are tried in order:
    #   1. nonuniform List<scalar> N ( ... ) ;
    #   2. uniform <value> ;
    # Whichever matches first is replaced.  A trailing append is only
    # used as a last resort when neither pattern matches (e.g. a skeleton
    # file with no internalField line at all).
    # ------------------------------------------------------------------

    # Pattern 1: nonuniform list (multi-line, values inside parens)
    pat_nonuniform = re.compile(
        r"internalField\s+nonuniform\s+List<\w+>\s*"
        r"\d+\s*\(\s*.*?\)\s*;",
        re.DOTALL,
    )

    # Pattern 2: uniform scalar (single line)
    pat_uniform = re.compile(
        r"internalField\s+uniform\s+[0-9eE+\-\.]+\s*;",
    )

    match = pat_nonuniform.search(original) or pat_uniform.search(original)

    if match:
        replacement = meta_comment + "\n\n" + new_internal
        patched = original[: match.start()] + replacement + original[match.end():]
        with open(of_path, "w") as fh:
            fh.write(patched)
        print(f"Replaced internalField in existing file -> {of_path}")
    else:
        # Last resort: no internalField found — insert before final separator
        print("Warning: no existing internalField block found — appending.")
        separator = "// ************************************************************************* //"
        if separator in original:
            insert_pos = original.rfind(separator)
            patched = (
                original[:insert_pos]
                + "\n"
                + meta_comment + "\n\n"
                + new_internal + "\n\n"
                + original[insert_pos:]
            )
        else:
            patched = original.rstrip() + "\n\n" + meta_comment + "\n\n" + new_internal + "\n"
        with open(of_path, "w") as fh:
            fh.write(patched)
        print(f"Appended internalField to existing file -> {of_path}")

else:
    # ------------------------------------------------------------------
    # FALLBACK: file does not exist — write a complete new file
    # ------------------------------------------------------------------
    of_content = f"""\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {OF_FIELD_NAME};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
{meta_comment}

dimensions      [0 1 -1 0 0 0 0];   // m/s

boundaryField
{{
    top
    {{
        type            zeroGradient;
    }}
    bottom
    {{
        type            zeroGradient;
    }}
    left
    {{
        type            zeroGradient;
    }}
    right
    {{
        type            zeroGradient;
    }}
    front
    {{
        type            zeroGradient;
    }}
    back
    {{
        type            zeroGradient;
    }}
}}

{new_internal}

// ************************************************************************* //
"""
    with open(of_path, "w") as fh:
        fh.write(of_content)
    print(f"No existing file found — wrote new OpenFOAM field → {of_path}")

print("Done.")

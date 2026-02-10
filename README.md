# TunnelWiz

TunnelWiz is a small Python package for **computing biomolecular tunnels** from structural data and **visualizing** the tunnel axis and tunnel-lining atoms.

**As of 10.02.2026 this package is not yet fully functional and published.**

It contains two modules:

- `builder` — tunnel geometry computation (axis + atom selection + coordinate transforms)
- `visual` — interactive visualization and plotting utilities

---

## Dependencies

Core scientific stack:

- NumPy
- SciPy
- MDAnalysis

Visualization and analysis utilities:

- Pandas
- Plotly
- NGLView
- RDKit *(currently used to guess charges from a PDB; this may change in future versions)*

---

## API Overview

### `get_axis(universe, endpoints, endpoint_type)`
Compute a tunnel axis and a local reference frame along that axis.

**Parameters**
- `universe` (`MDAnalysis.Universe`): structure loaded with MDAnalysis
- `endpoints` (list): two endpoints defining the tunnel. Supported formats depend on `endpoint_type`:
  - residue IDs (`resid`)
  - atom IDs (`atomid`)
  - Cartesian coordinates (`[x, y, z]`)
- `endpoint_type` (str): how to interpret `endpoints` (e.g. `"resid"`, `"atomid"`, `"coords"`)

**Returns**
- `axis` (`np.ndarray`): sampled axis points along the tunnel, shape `(N, 3)`
- `t` (`np.ndarray`): tangent vectors, shape `(N, 3)`
- `n` (`np.ndarray`): normal vectors, shape `(N, 3)`
- `b` (`np.ndarray`): binormal vectors, shape `(N, 3)`

---

### `get_atoms(universe, axes)`
Select atoms near the tunnel axis (i.e., tunnel-lining atoms).

**Parameters**
- `universe` (`MDAnalysis.Universe`)
- `axes`: output of `get_axis()` (axis + frame vectors)

**Returns**
- `atoms` (`MDAnalysis.core.groups.AtomGroup` or `MDAnalysis.Universe` depending on implementation):
  atoms selected around the tunnel axis

---

### `get_projection(points, axes)`
Project Cartesian coordinates into the tunnel coordinate system and return cylindrical coordinates.

**Parameters**
- `points` (`np.ndarray`): Cartesian coordinates, shape `(M, 3)`
- `axes`: output of `get_axis()`

**Returns**
- `z` (`np.ndarray`): position along the tunnel axis
- `theta` (`np.ndarray`): angular coordinate around the tunnel axis
- `r` (`np.ndarray`): radial distance from the tunnel axis

---

### `construct_tunnel(universe, endpoints, endpoint_type="resid")`
Convenience function that runs the full computation pipeline.

**What it does**
1. Compute the tunnel axis (`get_axis`)
2. Select tunnel-lining atoms (`get_atoms`)
3. Project those atoms into cylindrical coordinates (`get_projection`)

**Returns**
- `atoms`: tunnel-lining atoms
- `cyl`: cylindrical coordinates of those atoms (z, theta, r)

---

## Visualization & Plotting

### `show_tunnel(universe, axes)`
Visualize the tunnel axis inside the structure using NGLView.

**Parameters**
- `universe` (`MDAnalysis.Universe`)
- `axes`: output of `get_axis()`

**Returns**
- an `nglview.NGLWidget` (interactive viewer)

---

### `write_df(universe, cyl, pdb_path)`
Create a Pandas DataFrame for downstream plotting.

**Parameters**
- `universe`: atoms for which cylindrical coordinates were computed
- `cyl`: cylindrical coordinates from `get_projection()` / `construct_tunnel()`
- `pdb_path` (str or path): PDB file used to infer charges via RDKit *(temporary approach)*

**Returns**
- `df` (`pandas.DataFrame`): structured table for plotting and aggregation

---

### `show_scatter(df, color_by)`
Interactive scatter plot of tunnel-lining atoms in cylindrical space.

**Parameters**
- `df`: DataFrame created by `write_df()`
- `color_by` (str): column name used for coloring points, e.g.:
  - `"Chain"`, `"Charge"`, `"R"`, `"Type"`, `"ResName"`

---

### `show_heatmap(df, value)`
Heatmap of a chosen property in tunnel coordinates.

**Parameters**
- `df`: DataFrame created by `write_df()`
- `value` (str): currently supported: `"Charge"`

---

## Typical Pipeline

1. Load a structure using MDAnalysis:
   ```python
   import MDAnalysis as mda
   u = mda.Universe("structure.pdb")
2. Compute tunnel geometry and tunnel-lining atoms
    ```python
    atoms, cyl = construct_tunnel(u, endpoints=[10, 250], endpoint_type="resid")
3. Build a DataFrame for plotting
    ```python
    df = write_df(atoms, cyl, pdb_path="structure.pdb")
4. Visualize and plot
     ```python
    show_scatter(df, color_by="Chain")
    show_heatmap(df, value="Charge")

## Notes
Charge computation via RDKit + PDB is currently a placeholder and may change in future release. 

Selection of the tunnel-lining atoms is done by spherical selection. In future release this is intended to be improved.

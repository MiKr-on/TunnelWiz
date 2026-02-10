import pandas as pd
import plotly.express as px
import numpy as np
import MDAnalysis as mda


def show_tunnel(universe: mda.Universe, axes: np.array):
    """
    This function creates a quick overview of the tunnel axis within the biomolecules' structure.
    
    :param universe: biomolecule loaded with mda, or it's parts
    :param axes: array containing sampled axis points and tangent, normal, binormal vectors for each sample point
    """
    import nglview as nv
    view = nv.show_mdanalysis(universe)
    view.clear_representations()
    view.add_licorice()
    for i in range(len(axes[0])):
        view.shape.add_arrow(axes[0][i],axes[0][i]+axes[1][i]*5,[0,0,1],0.3)
        view.shape.add_arrow(axes[0][i],axes[0][i]+axes[2][i]*5,[0,1,0],0.3)
        view.shape.add_arrow(axes[0][i],axes[0][i]+axes[3][i]*5,[1,0,0],0.3)
    view


def make_df(universe: mda.Universe, cyl: np.array, pdb_path: np.array):
    """
    Create a dataframe with characteristics of atoms from the given universe, from those on the inner surface of the the tunnel

    :param universe: biomolecule loaded with mda
    :param 
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Guess atomic charges
    mol = Chem.MolFromPDBFile(pdb_path)
    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.ComputeGasteigerCharges(mol)
    charges = np.asarray([float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()])
    tunnel_charges = charges[universe.indices]

    Z     = cyl[0]
    Theta = cyl[1]
    R     = cyl[2]
    sizes = np.interp(R, (R.min(), R.max()), (1000, 100))

    df = pd.DataFrame({
        "Z": Z,
        "Theta": Theta,
        "ThetaDeg": Theta*180/3.14,
        "R": R,
        "Rinv": sizes,
        "Type": universe.types,
        "Chain": universe.chainIDs,
        "Resid": universe.resids,
        "ResName": universe.resnames,
        "Charge": tunnel_charges
        # hydrophobicity to be added
    })

    return df


def show_heatmap(df: pd.DataFrame, type: str="Charge", xbin_size: int=3, ybin_size: int=1):
    """
    Create a heatmap to visualize the inner surface of the computed tunnel, this functions only works for numerical description i.e. charges, hydrophobicity
    
    :param df: dataframe containg at least columns named Z, Theta, Charge
    :param type: characteristic to be plotted out, must match the column names in df
    """
    from scipy.stats import binned_statistic_2d

    x_bins = [i*xbin_size for i in range( int(df["Z"].max()/xbin_size)+2)]
    y_bins = [-3.5+i*ybin_size for i in range(8/ybin_size)]
    y_label = [3.5-i*ybin_size for i in range(8/ybin_size)]
    stat, x_edge, y_edge, binnum = binned_statistic_2d(df["Z"],df["Theta"],df[type],bins=[x_bins,y_bins])
    new_df = pd.DataFrame(stat).T
    new_df = new_df.set_axis(x_edge[1:],axis="columns")
    new_df = new_df.set_axis(y_label[1:],axis="index")

    fig = px.imshow(new_df,width=800,height=300)
    fig.update_layout(plot_bgcolor="white")
    fig.update_yaxes(
        title="Theta (Rad)",
        tickmode="array",
        tickvals=new_df.index,
        ticktext=y_label,
        linecolor="black",
        showgrid=False,
        mirror=True,
        linewidth=2
    )
    fig.update_xaxes(
        title="Z (Å)",
        tickmode="array",
        tickvals=new_df.columns,
        ticktext=new_df.columns,
        linecolor="black",
        showgrid=False,
        mirror=True,
        linewidth=2
        )
    fig.show()


def show_scatter(df: pd.DataFrame, type: str="Chain"):
    """
    Docstring pro show_scatter
    
    :param df: Popis
    :param type: Popis
    """

    fig = px.scatter(df, x="Z", y="Theta", color=f"{type}", size="Rinv", hover_data="Resid",width=800,height=300,opacity=1)
    fig.update_layout(plot_bgcolor="white")
    fig.update_yaxes(
        range=[-3.7,3.7],
        title="Theta (rad)",
        tickmode="array",
        tickvals=[0,3.14,-3.14],
        ticktext=["0","π","-π"],
        linecolor="black",
        showgrid=False,
        mirror=True,
        linewidth=2)
    fig.update_xaxes(
        title_text="Z (Å)",
        linecolor="black",
        showgrid=False,
        mirror=True,
        linewidth=2
        )
    fig.add_hline(y=0,line_width=2,line_color="gray",line_dash="dash")
    fig.show()
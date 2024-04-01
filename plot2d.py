from numpy import sin, cos, log10, diff
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import AxesGrid
from h5py import File
from typer import Typer

app = Typer()
u_range = dict(vmin=+0.0, vmax=10.0)
e_range = dict(vmin=-5.0, vmax=+0.0)
d_range = dict(vmin=-3.5, vmax=-0.5)
s_range = dict(vmin=-3.0, vmax=+2.2)


def plot_one_panel_frame(fig, filename):
    """
    Plot a one-panel image frame from a products file
    """
    print(f"load {filename}")
    with File(filename, "r") as h5f:
        r_faces = h5f["face_positions_i"][...]
        q_faces = h5f["face_positions_j"][...]
        u = h5f["radial_gamma_beta"][...]
        d = h5f["comoving_mass_density"][...]
        p = h5f["gas_pressure"][...]
        s = h5f["radial_momentum"][...] / diff(-cos(q_faces))
    e = p / d * 3.0
    x = r_faces * sin(q_faces)
    z = r_faces * cos(q_faces)
    grid = AxesGrid(
        fig,
        111,
        nrows_ncols=(1, 1),
        axes_pad=0.6,
        label_mode="L",
        share_all=False,
        cbar_location="top",
        cbar_mode="each",
        cbar_size="8%",
        cbar_pad="0%",
    )
    grid[0].set_aspect("equal")
    c0 = grid[0].pcolormesh(x, z, u, cmap="viridis", edgecolors="none")
    grid.cbar_axes[0].colorbar(c0)
    grid.cbar_axes[0].axis["top"].set_label(r"$\gamma \beta_r$")
    grid[0].set_xlabel(r"$x / r_{\rm shell}$")
    grid[0].set_ylabel(r"$z / r_{\rm shell}$")
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)


def plot_two_panel_frame(fig, filename):
    """
    Plot a two-panel image frame from a products file
    """
    print(f"load {filename}")
    with File(filename, "r") as h5f:
        r_faces = h5f["face_positions_i"][...]
        q_faces = h5f["face_positions_j"][...]
        u = h5f["radial_gamma_beta"][...]
        d = h5f["comoving_mass_density"][...]
        p = h5f["gas_pressure"][...]
    e = p / d * 3.0
    x = r_faces * sin(q_faces)
    z = r_faces * cos(q_faces)
    grid = AxesGrid(
        fig,
        111,
        nrows_ncols=(1, 2),
        axes_pad=0.10,
        label_mode="L",
        share_all=False,
        cbar_location="top",
        cbar_mode="edge",
        cbar_size="7%",
        cbar_pad="2%",
    )
    ax1 = grid[0]
    ax2 = grid[1]
    ax1.set_aspect("equal")
    ax2.set_aspect("equal")
    c1 = ax1.pcolormesh(x, z, u, cmap="viridis", **u_range)
    c2 = ax2.pcolormesh(x, z, log10(e), cmap="magma", **e_range)
    grid.cbar_axes[0].colorbar(c1)
    grid.cbar_axes[1].colorbar(c2)
    grid.cbar_axes[0].axis["top"].set_label(r"$\gamma \beta_r$")
    grid.cbar_axes[1].axis["top"].set_label(r"$e$")
    fig.subplots_adjust(left=0.05, right=0.95)


def plot_four_panel_frame(fig, filename):
    """
    Plot a four-panel image frame from a products file
    """
    print(f"load {filename}")
    with File(filename, "r") as h5f:
        r_faces = h5f["face_positions_i"][...]
        q_faces = h5f["face_positions_j"][...]
        u = h5f["radial_gamma_beta"][...]
        d = h5f["comoving_mass_density"][...]
        p = h5f["gas_pressure"][...]
        s = h5f["radial_momentum"][...] / diff(-cos(q_faces))
    e = p / d * 3.0
    x = r_faces * sin(q_faces)
    z = r_faces * cos(q_faces)
    grid = AxesGrid(
        fig,
        111,
        nrows_ncols=(2, 2),
        axes_pad=0.6,
        label_mode="L",
        share_all=False,
        cbar_location="top",
        cbar_mode="each",
        cbar_size="8%",
        cbar_pad="0%",
    )
    grid[0].set_aspect("equal")
    grid[1].set_aspect("equal")
    grid[2].set_aspect("equal")
    grid[3].set_aspect("equal")
    c0 = grid[0].pcolormesh(x, z, u, cmap="viridis")  # , **u_range)
    c1 = grid[1].pcolormesh(x, z, log10(e), cmap="magma")  # , **e_range)
    c2 = grid[2].pcolormesh(x, z, log10(d), cmap="plasma")  # , **d_range)
    c3 = grid[3].pcolormesh(x, z, s, cmap="viridis")  # , **s_range)
    grid.cbar_axes[0].colorbar(c0)
    grid.cbar_axes[1].colorbar(c1)
    grid.cbar_axes[2].colorbar(c2)
    grid.cbar_axes[3].colorbar(c3)
    grid.cbar_axes[0].axis["top"].set_label(r"$\gamma \beta_r$")
    grid.cbar_axes[1].axis["top"].set_label(r"$log_{10} e$")
    grid.cbar_axes[2].axis["top"].set_label(r"Comoving mass density $\rho$")
    grid.cbar_axes[3].axis["top"].set_label(r"Radial Momentum")
    grid[2].set_xlabel(r"$x / r_{\rm shell}$")
    grid[3].set_xlabel(r"$x / r_{\rm shell}$")
    grid[0].set_ylabel(r"$z / r_{\rm shell}$")
    grid[2].set_ylabel(r"$z / r_{\rm shell}$")
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)


@app.command()
def frame(filename: str, figsize: tuple[int, int] = (12, 10), panels: int = 4):
    fig = plt.figure(figsize=figsize)
    if panels == 1:
        plot_one_panel_frame(fig, filename)
    if panels == 4:
        plot_four_panel_frame(fig, filename)
    plt.show()


@app.command()
def movie(
    filenames: list[str],
    figsize: tuple[int, int] = (12, 10),
    title: str = "Movie",
    output: str = "out.mp4",
):
    FFMpegWriter = animation.writers["ffmpeg"]
    metadata = dict(title=title)
    writer = FFMpegWriter(fps=30, metadata=metadata)
    fig = plt.figure(figsize=figsize)
    with writer.saving(fig, output, dpi=200):
        for filename in filenames:
            plot_four_panel_frame(fig, filename)
            writer.grab_frame()
            fig.clf()
    print(f"write {output}")


if __name__ == "__main__":
    app()

from numpy import sin, cos, log10
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import AxesGrid
from h5py import File
from typer import Typer

app = Typer()


def plot_frame(fig, filename):
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
    c1 = ax1.pcolormesh(x, z, u, vmin=0.0, vmax=11.0, cmap="viridis")
    c2 = ax2.pcolormesh(x, z, e, cmap="magma", vmin=0.0, vmax=2.0)
    grid.cbar_axes[0].colorbar(c1)
    grid.cbar_axes[1].colorbar(c2)
    grid.cbar_axes[0].axis["top"].set_label(r"$\gamma \beta_r$")
    grid.cbar_axes[1].axis["top"].set_label(r"$e$")
    fig.subplots_adjust(left=0.05, right=0.95)


@app.command()
def frame(filename: str):
    fig = plt.figure(figsize=[10, 5])
    plot_frame(fig, filename)
    plt.show()


@app.command()
def movie(filenames: list[str], title: str = "Movie", output: str = "out.mp4"):
    FFMpegWriter = animation.writers["ffmpeg"]
    metadata = dict(title=title)
    writer = FFMpegWriter(fps=30, metadata=metadata)
    fig = plt.figure(figsize=[10, 5])
    with writer.saving(fig, output, dpi=200):
        for filename in filenames:
            plot_frame(fig, filename)
            writer.grab_frame()
            fig.clf()
    print(f"write {output}")


if __name__ == "__main__":
    app()

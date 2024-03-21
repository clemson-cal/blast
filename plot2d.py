from numpy import sin, cos
from matplotlib import pyplot as plt
from h5py import File
from typer import Typer

app = Typer()


@app.command()
def main(filename: str, field="radial_gamma_beta"):
    fig, ax1 = plt.subplots(figsize=[10, 10])
    with File(filename, "r") as h5f:
        r_faces = h5f["face_positions_i"][...]
        q_faces = h5f["face_positions_j"][...]
        y = h5f[field][...]
    x = r_faces * sin(q_faces)
    z = r_faces * cos(q_faces)
    c = ax1.pcolormesh(x, z, y)
    plt.colorbar(c)
    ax1.set_xlim(0.0, 0.3)
    ax1.set_ylim(0.0, 0.3)
    ax1.set_aspect("equal")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    app()

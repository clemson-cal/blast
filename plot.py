from pathlib import Path
import typer
from matplotlib import pyplot as plt
from h5py import File

app = typer.Typer()

@app.command()
def main(items: list[str]):
    prim = "comoving_mass_density"
    filenames = []

    for item in items:
        if item.endswith('.h5'):
            filenames.append(Path(item))
        else:
            prim = item

    fig, ax1 = plt.subplots()
    for filename in filenames:
        with File(filename, 'r') as h5f:
            xi, xo = h5f["__config__"]["domain"][...]
            x = h5f["cell_coordinate"][...]
            if prim in h5f:
                y = h5f[prim][...]
                ax1.plot(x, y, label=filename.name, color='k')
                # ax1.set_xlim(xi, xo)
                ax1.set_ylim(0.0)
                plt.show()
            else:
                print(f"Primitive quantity '{prim}' not found in file {filename}.")

if __name__ == "__main__":
    app()
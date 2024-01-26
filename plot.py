from pathlib import Path


def main(filenames: list[Path]):
    from matplotlib import pyplot as plt
    from h5py import File

    fig, ax1 = plt.subplots()
    for filename in filenames:
        h5f = File(filename)
        x = h5f["cell_coordinate"][...]
        d = h5f["comoving_mass_density"][...]
        ax1.plot(x, d, label=filename, color='k')
        ax1.set_xlim(1.0, 10.0)
        # ax1.set_ylim(0.0, 1.0)
    # ax1.legend()
    plt.show()


if __name__ == "__main__":
    from typer import Typer

    app = Typer(pretty_exceptions_enable=False)
    app.command()(main)

    try:
        app()
    except RuntimeError as e:
        print("Error:", e)

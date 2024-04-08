from argparse import ArgumentParser
from subprocess import call, DEVNULL
from itertools import chain
from h5py import File


def measure_energy(cfg, **kwargs):
    call(
        ["bin/blast2d_omp"] + list(f"{k}={v}" for k, v in dict(cfg, **kwargs).items()),
        stdout=DEVNULL,
    )
    h5f = File("chkpt.0000.h5")
    energy = h5f["state"]["cons"][...][:, :, 3]
    return energy.sum()


def main():
    parser = ArgumentParser()
    parser.add_argument("config_file")
    parser.add_argument("extra_configs", nargs="*")
    args = parser.parse_args()
    cfg = dict()
    with open(args.config_file) as f:
        for line in chain(f, args.extra_configs):
            name, var = line.partition("=")[::2]
            cfg[name.strip()] = var.strip()
    e0 = measure_energy(cfg, tfinal=0.0, shell_u=0.0)
    e1 = measure_energy(cfg, tfinal=0.0)
    print(f"The shell energy is {e1 - e0 : .4e}")


if __name__ == "__main__":
    main()

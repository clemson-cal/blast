from re import compile
from glob import glob
from subprocess import run
from numpy import arange, sin
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
from h5py import File


def invoke(fig, cfg):
    exe = run(["bin/blast_cpu", cfg], capture_output=True)
    r = compile("\[(\d+)\].*Mzps=([0-9]*\.*[0-9]*)")
    n = list()
    M = list()
    for line in str(exe.stdout, "utf-8").split("\n"):
        if m := r.match(line):
            n.append(int(m.groups()[0]))
            M.append(float(m.groups()[1]))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    h5f = File(sorted(glob("prods.????.h5"))[-1])
    x = h5f["cell_coordinate"][...]
    d = h5f["comoving_mass_density"][...]
    u = h5f["gamma_beta"][...]
    p = h5f["gas_pressure"][...]

    ax1.plot(x, d)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"Comoving mass density $\rho$")

    ax2.plot(x, p)
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"Gas pressure $p$")

    ax3.plot(n, M)
    ax3.set_xlabel(r"Iteration number")
    ax3.set_ylabel(r"Zone updates per second [$10^6$]")


def main():
    rc(group="text", usetex=True)
    with PdfPages("report.pdf") as pdf:
        for cfg in ["sod", "bmk"]:
            fig = plt.figure(figsize=[8, 6])
            invoke(fig, f"setups/{cfg}.cfg")
            fig.tight_layout()
            pdf.savefig(fig)
            fig.clf()

        d = pdf.infodict()
        d["Title"] = "Blast Code Diagnostics Report"
        d["Author"] = "Jonathan Zrake"
        d["Subject"] = str()
        d["Keywords"] = str()


if __name__ == "__main__":
    main()

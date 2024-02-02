import re, pathlib, platform
from pprint import pformat
from tempfile import TemporaryDirectory
from glob import glob
from subprocess import run, check_output
from numpy import arange, sin
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
from h5py import File


def invoke(fig, cfg, outdir):
    exe = run(["bin/blast_cpu", cfg, f"outdir={outdir}", "fold=1"], capture_output=True)
    r = re.compile("\[(\d+)\].*Mzps=([0-9]*\.*[0-9]*)")
    n = list()
    M = list()
    for line in str(exe.stdout, "utf-8").split("\n"):
        if m := r.match(line):
            n.append(int(m.groups()[0]))
            M.append(float(m.groups()[1]))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(323)
    ax3 = fig.add_subplot(325)
    ax4 = fig.add_subplot(222)
    ax5 = fig.add_subplot(224)
    h5f = File(sorted(glob(f"{outdir}/prods.????.h5"))[-1])
    x = h5f["cell_coordinate"][...]
    d = h5f["comoving_mass_density"][...]
    u = h5f["gamma_beta"][...]
    p = h5f["gas_pressure"][...]

    sha = check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    config = {k: h5f["__config__"][k][...] for k in h5f["__config__"]}
    description = f"""
        Code version: {sha}

        Platform:
        machine = {platform.machine()}
        platform = {platform.platform()}
        system = {platform.system()}
        processor = {platform.processor()}

        Run configuration:
        method = {str(config['method'], 'utf-8')}
        rk = {config['rk']}
        bc = {chr(config['bc'][0])}{chr(config['bc'][1])}
        dx = {config['dx']}
    """

    ax1.plot(x, d)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"Comoving mass density $\rho$")
    ax2.plot(x, u)
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"Gamma-beta $u$")
    ax3.plot(x, p)
    ax3.set_xlabel(r"$x$")
    ax3.set_ylabel(r"Gas pressure $p$")
    ax4.plot(n, M)
    ax4.set_xlabel(r"Iteration number")
    ax4.set_ylabel(r"Zone updates per second [$10^6$]")
    ax5.axes.xaxis.set_visible(False)
    ax5.axes.yaxis.set_visible(False)
    ax5.text(0.05, 0.0, description)


def main():
    rc(group="text", usetex=True)
    with PdfPages("report.pdf") as pdf:
        for cfg in glob("setups/*"):
            fig = plt.figure(figsize=[8, 6])
            with TemporaryDirectory() as outdir:
                print(cfg)
                invoke(fig, cfg, outdir)
            fig.suptitle(f"Setup: {pathlib.Path(cfg).stem}")
            fig.subplots_adjust(left=0.06, right=0.98, bottom=0.08)
            pdf.savefig(fig)
            fig.clf()

        d = pdf.infodict()
        d["Title"] = "Blast Code Diagnostics Report"
        d["Author"] = "Jonathan Zrake"
        d["Subject"] = str()
        d["Keywords"] = str()


if __name__ == "__main__":
    main()

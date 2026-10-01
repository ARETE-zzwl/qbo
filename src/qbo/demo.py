"""Run the input builders and boundary diagnostics on small synthetic datasets."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys

from netCDF4 import Dataset
import numpy as np

from .diagnostics import diagnose
from .io import sha256
from .targets import PROFILE_HPA


def synthetic_inputs(directory):
    directory.mkdir()
    pressure = PROFILE_HPA[::-1]
    lines = ["Synthetic demonstration data", "P (hPa): " + " ".join(map(str, pressure))]
    for year in range(1980, 2025):
        for month in range(1, 13):
            t = (year - 1980) * 12 + month - 1
            wind = 24 * np.cos(2 * np.pi * t / 28 - 1.8 * np.log(pressure / 70))
            lines.append(f"{year}{month:02d} " + " ".join(f"{u:.8f}" for u in wind))
    (directory / "wind.txt").write_text("\n".join(lines) + "\n", encoding="ascii")
    lat, lon = np.arange(-10., 11., 5.), np.arange(0., 360., 45.)
    for variable in ("sst", "sic"):
        path = directory / f"{variable}.nc"
        with Dataset(path, "w") as ds:
            for name, size in (("time", 528), ("latitude", len(lat)), ("longitude", len(lon))):
                ds.createDimension(name, size)
            time = ds.createVariable("time", "f8", ("time",))
            time.units, time.calendar = "days since 1981-01-01", "360_day"
            time[:] = np.arange(528) * 30
            ds.createVariable("latitude", "f8", ("latitude",))[:] = lat
            ds.createVariable("longitude", "f8", ("longitude",))[:] = lon
            field = ds.createVariable(variable, "f4", ("time", "latitude", "longitude"))
            if variable == "sst":
                values = 27 + .02 * (np.arange(528) // 12) + np.sin(2 * np.pi * np.arange(528) / 12)
                field[:] = np.broadcast_to(values[:, None, None], field.shape)
            else:
                field[:] = 0
            ds.source = "Synthetic demonstration data"
        with gzip.GzipFile(filename=str(path) + ".gz", mode="wb", mtime=0) as handle:
            handle.write(path.read_bytes())
    with Dataset(directory / "baseline.nc", "w") as ds:
        for name, size in (("time", 12), ("lat", len(lat)), ("lon", len(lon))):
            ds.createDimension(name, size)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        time = ds.createVariable("time", "f8", ("time",))
        time.units = "days since 0001-01-01"
        time[:] = np.arange(12) * 30
        ds.createVariable("date", "i4", ("time",))[:] = 10101 + np.arange(12) * 100
        ds.createVariable("datesec", "i4", ("time",))[:] = 0
        for name, value, units in (("SST_cpl", 27., "deg_C"), ("ice_cov", 0., "1")):
            field = ds.createVariable(name, "f4", ("time", "lat", "lon"))
            field.long_name, field.units = name, units
            field[:] = value
    with Dataset(directory / "snapshot.nc", "w") as ds:
        ds.input_kind = "synthetic"
        for name, size in (("diagnostic", 3), ("level", 5), ("lat", len(lat)), ("lon", len(lon))):
            ds.createDimension(name, size)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("layer_bottom_hpa", "f8", ("level",))[:] = [30, 50, 61.52, 73.75, 90]
        base = 85 + 12 * np.cos(np.deg2rad(lon))[None, :] + .3 * lat[:, None]
        ds.createVariable("tropopause_hpa", "f8", ("diagnostic", "lat", "lon"))[:] = np.stack([base, base + 4, base + 8])
        ds.createVariable("found", "i4", ("diagnostic", "lat", "lon"))[:] = 1
        ds.source = "Synthetic demonstration data"


def plot_overview(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with Dataset(output / "diagnostics/weights.nc") as ds:
        column = np.asarray(ds["column_weight"][:])
        zonal = np.asarray(ds["zonal_weight"][:])
        lat, lon = np.asarray(ds["lat"][:]), np.asarray(ds["lon"][:])
        bottom = np.asarray(ds["layer_bottom_hpa"][:, 0, 0])
    profiles = []
    for phase in ("W", "E"):
        with Dataset(output / f"targets/qbo_target_{phase}.nc") as ds:
            profiles.append(np.asarray(ds["qbo"][6]))
    np.savetxt(output / "march_profiles.csv", np.column_stack([PROFILE_HPA, *profiles]),
               delimiter=",", header="pressure_hpa,W_mps,E_mps", comments="")
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .6, "lines.linewidth": 1.2,
                         "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 78 / 25.4),
                                 gridspec_kw={"width_ratios": [1, 1.25, 1]}, layout="constrained")
        for values, label, color in zip(profiles, ("W target", "E target"), ("#D55E00", "#0072B2")):
            axes[0].plot(values, PROFILE_HPA, "o-", markersize=2.5, color=color, label=label)
        axes[0].axvline(0, color=".7", lw=.6)
        axes[0].set(xlabel="March zonal wind (m s$^{-1}$)", ylabel="Pressure (hPa)",
                    ylim=(105, 5), xlim=(-27, 27))
        axes[0].legend(frameon=False, fontsize=6.5, loc="lower left")
        heat = axes[1].pcolormesh(lon, lat, column[2], vmin=0, vmax=1, cmap="cividis", shading="nearest", rasterized=True)
        axes[1].set(xlabel="Longitude (°E)", ylabel="Latitude (°N)", xticks=[0, 90, 180, 270], yticks=[-10, 0, 10])
        bar = fig.colorbar(heat, ax=axes[1], orientation="horizontal", pad=.17, fraction=.08)
        bar.set_label("Column weight at 61.52 hPa")
        area = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], column.shape[1:])
        for values, label, color in ((column, "Column", "#009E73"), (zonal, "Zonal minimum", "#0072B2")):
            means = [np.average(v, weights=area) for v in values]
            axes[2].plot(means, bottom, "o-", markersize=3, label=label, color=color)
        axes[2].set(xlabel="Area-mean weight", ylabel="Layer-bottom pressure (hPa)",
                    xlim=(-.04, 1.04), ylim=(100, 20), xticks=[0, .5, 1])
        axes[2].legend(frameon=False, fontsize=6.5, loc="lower right")
        for ax, title in zip(axes, ("a  Target profiles", "b  Column support", "c  Zonal constraint")):
            ax.set_title(title, loc="left", fontsize=8, fontweight="bold", pad=7)
        fig.suptitle("Synthetic QBO demonstration", fontsize=9, fontweight="bold")
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(output / f"overview.{suffix}", dpi=300)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    synthetic_inputs(inputs)
    commands = [
        ("targets", "--input", inputs / "wind.txt", "--output", output / "targets"),
        ("backgrounds", "--sst-gz", inputs / "sst.nc.gz", "--ice-gz", inputs / "sic.nc.gz",
         "--baseline", inputs / "baseline.nc", "--output", output / "backgrounds"),
    ]
    for module, *arguments in commands:
        subprocess.run([sys.executable, "-m", f"qbo.{module}", *map(str, arguments)], check=True, timeout=60)
    # Keep the provenance of these demonstrator files explicit when copied elsewhere.
    for path in [*output.glob("targets/*.nc"), *output.glob("backgrounds/*.nc")]:
        with Dataset(path, "a") as ds:
            ds.source = "Synthetic demonstration data generated by qbo-demo"
            ds.title = "Synthetic demonstration: " + ds.title
    for path in (output / "targets/target_metadata.json", output / "backgrounds/background_metadata.json"):
        metadata = json.loads(path.read_text())
        metadata["input_kind"] = "synthetic"
        if "target_files" in metadata:
            metadata["target_files"] = {f"{phase}QBO": sha256(output / f"targets/qbo_target_{phase}.nc") for phase in ("W", "E")}
        else:
            for item in metadata["outputs"].values():
                item["sha256"] = sha256(Path(item["path"]))
        path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    subprocess.run([sys.executable, "-m", "qbo.validate", str(output / "targets"),
                    "--output", str(output / "targets/qc.json")], check=True, timeout=60)
    diagnose(inputs / "snapshot.nc", output / "diagnostics")
    plot_overview(output)
    (output / "README.md").write_text(
        "# Synthetic QBO demonstration\n\n"
        "Generated locally with `qbo-demo`. All source fields are synthetic.\n\n"
        "![Targets and boundary weights](overview.png)\n\n"
        "`targets/` contains seven-level W/E composites and file QC. "
        "`backgrounds/` contains the two monthly climatologies. "
        "`diagnostics/weights.nc`, `diagnostics/summary.csv` and `march_profiles.csv` "
        "contain the plotted data. The figure is also saved as editable SVG and PDF.\n\n"
        "The zonal minimum is uniform around each latitude ring and never exceeds a column weight. "
        "The weights describe the operator's support; realized winds require model output.\n",
        encoding="utf-8",
    )
    print(f"Demo saved to {output}")


if __name__ == "__main__":
    main()

"""multi-run / folder-run submit through AaronTools.job_control, which needs jinja2 and a job template."""
import os
import stat

from hfrpkg.conformers import ConformerSettings
from hfrpkg.run_single import run_jobs
from hfrpkg.runner import run_reaction


def test_folder_run_writes_and_submits_slurm_jobs(tmp_path, monkeypatch):
    import AaronTools.job_control as jc

    lib, bin_ = tmp_path / "lib", tmp_path / "bin"
    lib.mkdir()
    bin_.mkdir()
    (lib / "Gaussian_template.txt").write_text(
        "#SBATCH --job-name={{ name }} --cpus-per-task={{ processors }} --mem={{ memory }}G --time={{ walltime }}:00:00\n")
    sbatch = bin_ / "sbatch"
    sbatch.write_text("#!/bin/sh\necho \"$1\" >> submitted.txt\n")
    sbatch.chmod(sbatch.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setattr(jc, "QUEUE_TYPE", "SLURM")
    monkeypatch.setattr(jc, "AARONLIB", str(lib))

    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCC", outfolder=str(out), extension="g",
                 conformers=ConformerSettings(enabled=False))
    run_jobs(str(out), walltime=24, processors=12, memory=24)

    coms = sorted(p.stem for p in out.glob("*.com"))
    assert sorted(os.path.basename(p) for p in (out / "submitted.txt").read_text().split()) == sorted(c + ".job" for c in coms)
    job = (out / f"{coms[0]}.job").read_text()
    assert f"--job-name={coms[0]}" in job and "--cpus-per-task=12" in job and "--mem=24G" in job and "--time=24:00:00" in job

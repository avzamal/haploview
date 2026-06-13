"""Parity test against the original Haploview jar.

Runs the original Haploview 4.1 and this tool on the same committed sample data
and asserts the LD table and block membership are identical. Skips cleanly when
Java or the jar are unavailable (e.g. CI without the downloaded jar).
"""

import os
import shutil
import subprocess

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA = os.path.join(ROOT, "data")
JAR = os.path.join(ROOT, "third_party", "Haploview4.1.jar")

pytestmark = pytest.mark.skipif(
    shutil.which("java") is None or not os.path.exists(JAR),
    reason="Java or Haploview4.1.jar not available")


def _parse_ld(path):
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 7:
                continue
            out[(t[0], t[1])] = tuple(float(x) for x in t[2:7])
    return out


def _hv_blocks(path, info_names):
    blocks = []
    if not os.path.exists(path):
        return blocks
    for line in open(path):
        if line.startswith("BLOCK"):
            nums = [t.rstrip("!") for t in line.split("MARKERS:")[1].split()]
            blocks.append([info_names[int(n) - 1] for n in nums])
    return blocks


def test_parity_on_sample(tmp_path):
    from haploview.pipeline import run_analysis, write_outputs

    hv_out = str(tmp_path / "hv")
    subprocess.run(
        ["java", "-Djava.awt.headless=true", "-jar", JAR, "-n",
         "-pedfile", os.path.join(DATA, "sample.ped"),
         "-info", os.path.join(DATA, "sample.info"),
         "-dprime", "-blockoutput", "GAB", "-maxdistance", "500",
         "-out", hv_out],
        check=True, cwd=str(tmp_path), capture_output=True)

    result = run_analysis(os.path.join(DATA, "sample.vcf"), method="gabriel",
                          max_distance_kb=500)
    mine = str(tmp_path / "mine")
    write_outputs(result, mine, "gabriel")

    hv_ld = _parse_ld(hv_out + ".LD")
    my_ld = _parse_ld(mine + ".LD")
    assert set(hv_ld) == set(my_ld)
    for key in hv_ld:
        for i in range(5):
            assert abs(hv_ld[key][i] - my_ld[key][i]) < 1e-9, (key, i)

    info_names = [l.split()[0] for l in open(os.path.join(DATA, "sample.info"))]
    hv_blocks = _hv_blocks(hv_out + ".GABRIELblocks", info_names)
    my_blocks = [[result.dataset.markers[m].name for m in b]
                 for b in result.blocks]
    assert hv_blocks == my_blocks

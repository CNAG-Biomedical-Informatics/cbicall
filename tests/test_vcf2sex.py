import gzip
import os
import subprocess

import pytest


REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
VCF2SEX_SCRIPTS = [
    os.path.join(REPO_ROOT, "workflows", "bash", "gatk-3.5", "vcf2sex.sh"),
    os.path.join(REPO_ROOT, "workflows", "bash", "gatk-4.6", "vcf2sex.sh"),
    os.path.join(REPO_ROOT, "workflows", "snakemake", "gatk-4.6", "vcf2sex.sh"),
]


def _write_vcf(path, records):
    text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            *records,
            "",
        ]
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _run_vcf2sex(script, tmp_path, records):
    vcf = tmp_path / "sample.vcf.gz"
    _write_vcf(vcf, records)

    env_file = tmp_path / "env.sh"
    env_file.write_text("ARCH=x86_64\n", encoding="utf-8")

    env = os.environ.copy()
    env["CBICALL_ENV_FILE"] = str(env_file)

    return subprocess.run(
        [script, str(vcf)],
        check=False,
        text=True,
        capture_output=True,
        env=env,
    )


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_reports_unknown_when_x_y_depths_are_absent(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "22\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t0/1:1,1:2",
            "22\t200\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t0/1:2,2:4",
        ],
    )

    assert result.returncode == 0
    assert "Illegal division by zero" not in result.stderr
    assert result.stdout.splitlines() == [
        "MEAN DEPTH FOR AUTOSOMES=3.00",
        "MEAN DEPTH FOR X=NA",
        "MEAN DEPTH FOR Y=NA",
        "THRESHOLD=NA",
        "SEX=UNKNOWN",
    ]


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_uses_dp_format_position(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:60:30,30",
            "X\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20",
            "Y\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:1:1,0",
        ],
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "MEAN DEPTH FOR AUTOSOMES=60.00",
        "MEAN DEPTH FOR X=40.00",
        "MEAN DEPTH FOR Y=1.00",
        "THRESHOLD=20",
        "SEX=FEMALE",
    ]

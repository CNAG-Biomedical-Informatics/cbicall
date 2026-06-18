import gzip
import os
import subprocess

import pytest


REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
METHOD_COMMENT = "# METHOD: sex inferred for QC from VCF-derived depth proxies. If X/autosome ratio >= 0.75, classify FEMALE; otherwise classify FEMALE only when X-Y depth difference >= THRESHOLD."

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
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=3.00",
        "MEAN DEPTH FOR X=NA",
        "MEAN DEPTH FOR Y=NA",
        "X_AUTOSOME_RATIO=NA",
        "X_MINUS_Y_DEPTH=NA",
        "THRESHOLD=NA",
        "DECISION=Insufficient autosome, X, or Y depth values",
        "SEX=UNKNOWN",
    ]


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_reports_likely_female_when_y_depth_is_absent(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:100",
            "2\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:102",
            "X\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:63",
        ],
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=101.00",
        "MEAN DEPTH FOR X=63.00",
        "MEAN DEPTH FOR Y=NA",
        "X_AUTOSOME_RATIO=0.62",
        "X_MINUS_Y_DEPTH=NA",
        "THRESHOLD=NA",
        "DECISION=X/autosome ratio >= 0.55 and no usable Y records",
        "SEX=FEMALE_LIKELY",
    ]


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_keeps_unknown_when_y_depth_is_absent_and_x_ratio_is_low(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:100",
            "2\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:102",
            "X\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:50",
        ],
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=101.00",
        "MEAN DEPTH FOR X=50.00",
        "MEAN DEPTH FOR Y=NA",
        "X_AUTOSOME_RATIO=0.50",
        "X_MINUS_Y_DEPTH=NA",
        "THRESHOLD=NA",
        "DECISION=Insufficient autosome, X, or Y depth values",
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
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=60.00",
        "MEAN DEPTH FOR X=40.00",
        "MEAN DEPTH FOR Y=1.00",
        "X_AUTOSOME_RATIO=0.67",
        "X_MINUS_Y_DEPTH=39",
        "THRESHOLD=20",
        "DECISION=X-Y depth difference >= threshold",
        "SEX=FEMALE",
    ]


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_uses_x_autosome_ratio_when_y_variant_depth_is_noisy(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:31",
            "2\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:31",
            "X\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:27",
            "Y\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:44",
        ],
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=31.00",
        "MEAN DEPTH FOR X=27.00",
        "MEAN DEPTH FOR Y=44.00",
        "X_AUTOSOME_RATIO=0.87",
        "X_MINUS_Y_DEPTH=-17",
        "THRESHOLD=8",
        "DECISION=X/autosome ratio >= 0.75",
        "SEX=FEMALE",
    ]


@pytest.mark.parametrize("script", VCF2SEX_SCRIPTS)
def test_vcf2sex_keeps_male_when_x_autosome_ratio_is_low(script, tmp_path):
    result = _run_vcf2sex(
        script,
        tmp_path,
        [
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:31",
            "2\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:31",
            "X\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15",
            "Y\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:28",
        ],
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        METHOD_COMMENT,
        "MEAN DEPTH FOR AUTOSOMES=31.00",
        "MEAN DEPTH FOR X=15.00",
        "MEAN DEPTH FOR Y=28.00",
        "X_AUTOSOME_RATIO=0.48",
        "X_MINUS_Y_DEPTH=-13",
        "THRESHOLD=8",
        "DECISION=X/autosome ratio < 0.75 and X-Y depth difference < threshold",
        "SEX=MALE",
    ]

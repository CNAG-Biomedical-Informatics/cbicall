import os, gdown

# Dictionary mapping the desired output filenames to their Google Drive file IDs
files = {
    "data.tar.gz.md5":      "1KPQdljTv-jDk8kodRgNGcajby97chva_",
    "data.tar.gz.part-00":  "1xNX0vH5Ya4XqhQEylk_ZZNVtpEeCsZD9",
    "data.tar.gz.part-01":  "11sJTFojJkgUyrKMROAQihP4K9z8515Kt",
    "data.tar.gz.part-02":  "1jkUMU0iyDMeZ_jkad_24yk1YcTcYghM2",
    "data.tar.gz.part-03":  "1U2iZ1kzQ5vwbj3FuTcMJ2Vl2Oxegg0C3",
    "data.tar.gz.part-04":  "1zgCjShM0T6fOAMEirbQdoLT03V7-1RxR",
    "data.tar.gz.part-05":  "1faW5N8NECugJddJyoTfpn5rtkb1tNygl",
}

def download_if_missing(filename, file_id):
    if os.path.exists(filename):
        print(f"{filename} already exists. Skipping download.")
    else:
        url = f'https://drive.google.com/uc?export=download&id={file_id}'
        print(f"Downloading {filename}...")
        gdown.download(url, filename, quiet=False)

# Check and download files if not present
for filename, file_id in files.items():
    download_if_missing(filename, file_id)

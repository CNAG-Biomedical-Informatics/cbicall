import os, gdown

# Dictionary mapping the desired output filenames to their Google Drive file IDs
files = {
    'data.tar.gz.md5': '1gCnKxcM1NIJk2LnUdHqYQ1zRChpUVoaY',
    'data.tar.gz.part-00': '1TznaIvxCrcvuh3Zfce4HYv2E_hlhNAZE',
    'data.tar.gz.part-01': '16ufn1qhTkBVXMJe66S4WTXAms61i2xBy',
    'data.tar.gz.part-02': '1WWxX72m-1vOKe8DUVkUYhKnyMUZ7IaUy',
    'data.tar.gz.part-03': '1BpFuiNqygOi8Ir2CqVDxDIXWXnq11-ln',
    'data.tar.gz.part-04': '1tAuti2nLUM6QPEkP_UfDuVXyZOrLjE33'
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

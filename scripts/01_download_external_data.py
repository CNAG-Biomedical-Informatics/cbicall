import os, gdown

# Dictionary mapping the desired output filenames to their Google Drive file IDs
files = {
    'data.tar.gz.md5': '1m8Xge4bMPFWkuwoaLvww9DYj6DmFG8_B',
    'data.tar.gz.part-00': '1zg03tDnQD-arH8fK7Jic1_aQIDyTP0XP',
    'data.tar.gz.part-01': '1_fPH0_Fes9gDYnRI682oO_9Sg6S17NiZ',
    'data.tar.gz.part-02': '1wOxB25iuCK02USBTP5RK2gctWl07zLH8',
    'data.tar.gz.part-03': '19afh4MOr8oGMk0N8XIJJj3GB8bNUHFml'
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

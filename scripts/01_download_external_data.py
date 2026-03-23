import os, gdown

# Dictionary mapping the desired output filenames to their Google Drive file IDs
files = {
    "data.tar.gz.md5": "1MmeppkZ1xR9ODsBLuDWUh-2ob98Prxq6",
    "data.tar.gz.part-00": "12HQw_duxmlcASh7Yw8yL2-f-hyjqMim9",
    "data.tar.gz.part-01": "1Ejrn1oQ2WO3TSTmWSfoQ1IjYBOXivbRm",
    "data.tar.gz.part-02": "18CC1y2t5heV66-jOHCFIHqgj9GbT3QgF",
    "data.tar.gz.part-03": "1mq5t0U8GFk6ZdaMqBzWHD73LfYJ626vi",
    "data.tar.gz.part-04": "1v8rQA-qtgYgnPcOV8PzO-otw9dE9ePqt",
    "data.tar.gz.part-05": "1knI_DEdrlushYj5lOclnRNDdvSEmXQ5m",
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

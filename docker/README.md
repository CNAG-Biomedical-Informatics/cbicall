# Containerized Installation

## Downloading Required Databases and Software

Begin by downloading the required databases and software. Save the data outside the container; this preserves it across container restarts and lets you update the software without downloading the data again.

Install dependencies for Python 3:

```
pip3 install gdown
```

Finally, navigate to a directory where you want the databases stored and execute:

```
wget https://raw.githubusercontent.com/mrueda/cbicall/refs/heads/main/scripts/01_download_external_data.py
python3 ./01_download_external_data.py
```

Note: Google Drive can be a tad restrictive with the download. If you get an error, please use the error URL link in a browser and you should be able to retrieve it there.

Once downloaded, perform a checksum to make sure the files were not corrupted:

```
md5sum -c data.tar.gz.md5
```

Now let's reassemble the split files into the original tar archive:

```
cat data.tar.gz.part-?? > data.tar.gz
```

Clean up split files to save space (when you think you are ready!):

```
rm data.tar.gz.part-??
```

Extract the tar archive:

```
tar -xzvf data.tar.gz
```

Finally, in the `cbicall` repo:

Change `DATADIR` variable in `workflows/bash/*/parameters.sh` and `workflows/snakemake/*/config.yaml` so that it matches the location of your downloaded data.


### Method 2: Installing from Docker Hub (fast)

Pull the latest Docker image from [Docker Hub](https://hub.docker.com/r/manuelrueda/cbicall):

```bash
docker pull manuelrueda/cbicall:latest
docker image tag manuelrueda/cbicall:latest cnag/cbicall:latest
```

### Method 3: Installing from Dockerfile (slow)

Download the `Dockerfile` from [GitHub](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/Dockerfile):

```bash
wget https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/main/docker/Dockerfile
```

Then build the container:

- **For Docker version 19.03 and above (supports buildx):**

  ```bash
  docker buildx build -t cnag/cbicall:latest .
  ```

- **For Docker versions older than 19.03 (no buildx support):**

  ```bash
  docker build -t cnag/cbicall:latest .
  ```

## Running and Interacting with the Container

```bash
# Please update '/absolute/path/to/cbicall-data' with your actual local data path
docker run -tid --volume /absolute/path/to/cbicall-data:/cbicall-data -e USERNAME=root --name cbicall cnag/cbicall:latest
```

To connect to the container:

```bash
docker exec -ti cbicall bash
```

Finally, inside the `cbicall` repo:

Change `DATADIR` variable in `workflows/bash/parameters.sh` and `workflows/snakemake/config.yaml` to `/cbicall-data`.

## Performing unit test

Inside the container

```bash
cd examples/input
./run_test.sh
```

## System requirements

- OS/ARCH supported: **linux/amd64** and **linux/arm64**.
- Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well (untested).
- 16GB of RAM
- \>= 1 core (ideally i7 or Xeon).
- At least 100GB HDD.

## Common errors: Symptoms and treatment

  * Dockerfile:

          * DNS errors

            - Error: Temporary failure resolving 'foo'

              Solution: https://askubuntu.com/questions/91543/apt-get-update-fails-to-fetch-files-temporary-failure-resolving-error

# Non-containerized installation

### Method 1: Download from GitHub

First, we need to install a few system components (note you might have them already in your system):

```bash
sudo apt install gcc make git cpanminus libperl-dev
```

Use `git clone` to get the latest (stable) version:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
```

If you only new to update to the lastest version do:

```bash
git pull
```

We use `cpanm` to install the CPAN modules. We'll install the dependencies at `~/perl5`:

```bash
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
cpanm --notest --installdeps .
bin/cbicall
```
Testing the deployment:

```bash
prove
```

To ensure Perl recognizes your local modules every time you start a new terminal, run:

```bash
echo 'eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)' >> ~/.bashrc
```

### Install required external software

Install dependencies for Python 3:

```
pip3 install -r requirements.txt
```

Finally, navigate to a directory where you want the databases stored and execute:

```
python3 $path_to_cbicall/scripts/01_download_external_data.py  # Replace $path_to_cbicall with your CBICall installation path.
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

Ok, finally we are going to install `Java 8` in case you don't have it already:

```
sudo apt install openjdk-8-jdk # In some systems you might need Java 17 -> openjdk-17-jre
```

## Performing unit test

Once you are in the root directory of the repo:

```bash
cd examples/input
./run_test.sh
```

## System requirements

- OS/ARCH supported: **linux/amd64** and **linux/arm64**.
- Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well (untested).
- Perl 5 (>= 5.36 core; installed by default in many Linux distributions). Check the version with `perl -v`
- Java 8
- 16GB of RAM
- \>= 1 core (ideally i7 or Xeon).
- At least 100GB HDD.

## Platform Compatibility
This distribution is written in pure Perl and is intended to run on any platform supported by Perl 5. It has been tested on Debian Linux and macOS. Please report any issues.

## Common errors: Symptoms and treatment

* Perl errors:
    - Foo

      Solution: 

      `Bar`

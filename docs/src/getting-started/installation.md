# Installation

microcosm runs with [Nextflow](https://www.nextflow.io/) and gets every tool from containers or Conda: you install Nextflow and one software runtime, then the databases you need ([Databases](databases.md)).

### Software

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) **version 26.04.0 or later** (required since microcosm 2.0.0, which uses the Nextflow strict syntax)
- [Java](http://www.oracle.com/technetwork/java/javase/downloads/index.html) (version 16 to 22, required by Nextflow)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and/or [Mamba](https://github.com/conda-forge/miniforge) for package and environment management.
- Optional: A container runtime engine like [Docker](https://www.docker.com/) or [Apptainer](https://apptainer.org/docs/admin/main/installation.html).

### Hardware

- A POSIX-compatible operating system (e.g., Linux, macOS) or Windows with [WSL](https://learn.microsoft.com/en-us/windows/wsl/). We strongly recommend using a Unix-based OS for best compatibility with all tools.
- At least 32 GB of RAM for MEGAHIT, or 128 GB if you plan to use metaSPAdes for assembly.
- At least 256 GB of disk space.

  > ℹ️ **Storage Estimates** (based on shallow shotgun metagenomics):
  >
  > - Nextflow installation: ~150 MB
  > - Databases: ~120 GB
  > - Environment folder: ~20 GB
  > - Intermediate/cache files: 20–50 GB (varies with sample size and quality)
  > - Output per sample: ~9 GB

* Optional: Some steps in assembly-based analysis (for example, `COMEBin` in the binning step) will run faster if your system is equipped with a GPU.


## Install the pipeline

### 1. Clone the Pipeline

```bash
git clone --single-branch --branch main https://github.com/tnmquann/microcosm.git
```
- Or download the [latest release](https://github.com/tnmquann/microcosm/releases/latest).
- Or use this command to download and unzip in one step:
```bash
nextflow pull tnmquann/microcosm
# Use this command to check pipeline info (if necessary)
nextflow info tnmquann/microcosm
```

### 2. Install Nextflow

**a.** Check Java installation:
```bash
java -version
```
If missing, follow [Java install instructions](https://www.nextflow.io/docs/latest/install.html).

**b.** Download Nextflow:
```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
```

**c.** (Optional) Move Nextflow to your `$PATH`:
```bash
sudo mv nextflow /usr/local/bin
```

**d.** Test Nextflow:
```bash
nextflow info
```

### 3. Install Conda/Mamba

- **Mamba (Recommended):** Install via [Miniforge](https://github.com/conda-forge/miniforge).
- **Conda:** [Install guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### 4. Install Container runtime (Optional)

- **Apptainer:** [Install guide](https://apptainer.org/docs/admin/main/installation.html). Users are recommended to use the unprivileged version to avoid unexpected errors.
- **Docker:** [Install guide](https://docs.docker.com/engine/install/ubuntu/)

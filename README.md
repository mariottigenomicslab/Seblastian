# Seblastian

Seblastian predicts eukaryotic selenoprotein genes by searching genomic sequences for SECIS elements and then analysing the upstream regions with BLAST and Exonerate.

The workflow was described in:

> Mariotti M, Lobanov AV, Guigó R, Gladyshev VN. SECISearch3 and Seblastian: new tools for prediction of SECIS elements and selenoproteins. *Nucleic Acids Research*. 2013;41(15):e149. [https://doi.org/10.1093/nar/gkt550](https://doi.org/10.1093/nar/gkt550)

Seblastian is legacy software written for Python 2.7. Linux x86-64 is the best-supported platform, particularly for the bundled SECISearch and COVE executables.

## Installation with Conda

Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or Mamba, then clone this repository:

```bash
git clone https://github.com/mariottigenomicslab/Seblastian.git
cd Seblastian
```

Create and activate an environment containing Python 2.7 and the external command-line dependencies:

```bash
conda create -n seblastian \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority \
  python=2.7 \
  perl \
  gawk \
  blast-legacy=2.2.26 \
  exonerate=2.4.0 \
  infernal=1.0.2 \
  viennarna=2.4.18

conda activate seblastian
```

ViennaRNA 2.4.18 is used here because it is the latest Bioconda build compatible with Python 2.7. Make the bundled scripts discoverable and create a writable temporary directory:

```bash
export SEBLASTIAN_HOME="$(pwd)"
export PATH="$SEBLASTIAN_HOME:$SEBLASTIAN_HOME/bin:$PATH"
chmod +x "$SEBLASTIAN_HOME/blaster_parser.g"
mkdir -p "$SEBLASTIAN_HOME/tmp"
```

Display the command-line help:

```bash
python "$SEBLASTIAN_HOME/Seblastian.py" -h
```

The environment variable and `PATH` change apply to the current shell. Set them again after opening a new terminal.

## Running Seblastian

### SECIS search only

```bash
mkdir -p results

python "$SEBLASTIAN_HOME/Seblastian.py" \
  -t target.fa \
  -o results/seblastian \
  -SS \
  -temp "$SEBLASTIAN_HOME/tmp" \
  -infernal_cm "$SEBLASTIAN_HOME/SECIS_infernal.cm" \
  -infernal_stk "$SEBLASTIAN_HOME/SECIS_infernal.stk" \
  -covels_cm "$SEBLASTIAN_HOME/SECIS_covels.cm" \
  -bin_folder "$SEBLASTIAN_HOME/bin"
```

### Full selenoprotein prediction

```bash
mkdir -p results

python "$SEBLASTIAN_HOME/Seblastian.py" \
  -t target.fa \
  -o results/seblastian \
  -d "$SEBLASTIAN_HOME/uniref50.only_selenoproteins.fa" \
  -temp "$SEBLASTIAN_HOME/tmp" \
  -infernal_cm "$SEBLASTIAN_HOME/SECIS_infernal.cm" \
  -infernal_stk "$SEBLASTIAN_HOME/SECIS_infernal.stk" \
  -covels_cm "$SEBLASTIAN_HOME/SECIS_covels.cm" \
  -bin_folder "$SEBLASTIAN_HOME/bin"
```

Replace `target.fa` with a nucleotide FASTA file. Run `python "$SEBLASTIAN_HOME/Seblastian.py" -h` to see all available routines and options.

## Alternative: Docker

A prebuilt image is available from [Docker Hub](https://hub.docker.com/r/maxtico/seblastian):

```bash
docker pull maxtico/seblastian:latest
```

Display the help:

```bash
docker run --rm maxtico/seblastian:latest \
  python /Seblastian/Seblastian.py -h
```

Run SECIS prediction on `target.fa` in the current directory and write the results back to that directory:

```bash
docker run --rm \
  --mount type=bind,source="$(pwd)",target=/data \
  maxtico/seblastian:latest \
  python /Seblastian/Seblastian.py \
  -t /data/target.fa \
  -o /data/seblastian \
  -SS
```

Run the full pipeline with the protein database included in the image:

```bash
docker run --rm \
  --mount type=bind,source="$(pwd)",target=/data \
  maxtico/seblastian:latest \
  python /Seblastian/Seblastian.py \
  -t /data/target.fa \
  -o /data/seblastian
```

## Models and protein database

The repository currently includes the resources used by the program:

- `SECIS_infernal.cm`, `SECIS_infernal.1.0.2.cm`, and `SECIS_covels.cm` are covariance models.
- `SECIS_infernal.stk` is the alignment associated with the Infernal model.
- `uniref50.only_selenoproteins.fa` is the bundled protein database.
- The `.phr`, `.pin`, and `.psq` files are the prebuilt legacy BLAST indexes for that protein database.

The bundled protein database and its indexes occupy less than 1 MB, so they are kept in this repository. A different protein FASTA database can be supplied with `-d`; Seblastian will create legacy BLAST indexes when they are absent.

If a future database is too large for GitHub, publish a versioned archive in a data repository such as Zenodo and document its DOI, checksum, and expected filename here.

## Contributing

Changes can be proposed through GitHub issues and pull requests. When changing the prediction workflow, include the command and input used to verify the result.

## License

See [LICENSE](LICENSE).

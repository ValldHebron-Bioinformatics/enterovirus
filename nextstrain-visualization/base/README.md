# Template

Template for basic EV builds based on the [template repository](https://github.com/hodcroftlab/template_nextstrain) from Emma Hodcroft lab.


## Prerequisites
Ensure you have the following installed:
- Python=3.8 or higher
- Micromamba or Conda
- Snakemake=7: [Instalation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Nextstrain CLI: [Instalation guide](https://docs.nextstrain.org/en/latest/install.html)

## Structure

You should have the following folders inside your [ EV ] folder:
- config
- data

### config 
This folder contains the following files:
- `auspice_config.json`: Configuration file. Replace in the first line < EV > with your virus
- `colors.tsv`: Color settings.
- `geo_regions.tsv`: Association between country and regions.
- `lat_longs.tsv`: Coord for maps.
- `reference_sequence.gb`: Reference sequence for your EV. Check that in the proteins the codification is CDS.

### data
This folder contains the following files:
- `metadata.tsv`: Contains metadata used for filtering and coloring.
- `sequences.fasta`: Contains the sequences of the virus available in the database.

## Snakemake build
nextstrain build --cpus 8 [EV]
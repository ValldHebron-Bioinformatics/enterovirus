# Template

Template for basic EV builds based on the [template repository](https://github.com/hodcroftlab/template_nextstrain) from Emma Hodcroft lab.


## ⚙️ Prerequisites
Ensure you have the following installed:
- Python=3.8 or higher
- Snakemake=7: [Instalation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Augur=30.0.0: [Instalation guide](https://docs.nextstrain.org/projects/augur/en/stable/installation/installation.html), [Git repo](https://github.com/nextstrain/augur.git)
- Dependencies needed by augur, listed in setup.py file from its [git repo](https://github.com/nextstrain/augur.git)
- IQ-TREE=2.0.7
- MAFFT=7.505 

## 🛠 Structure

You should have the following inside your < EV > folder:
- 📁 config
- 📁 data
- 📄 snakefile

### 📁 config 
This folder contains the following files:
- `auspice_config.json`: Configuration file. Replace in the first line < EV > with your virus
- `colors.tsv`: Color settings.
- `geo_regions.tsv`: Association between country and regions.
- `lat_longs.tsv`: Coord for maps.
- `reference_sequence.gb`: Reference sequence for your EV. Check that in the proteins the codification is CDS.

### 📁 data
This folder contains the following files:
- `metadata.tsv`: Contains metadata used for filtering and coloring.
- `sequences.fasta`: Contains the sequences of the virus available in the database.

### 📄 snakefile
This file will create the json for visualization. You have to:
- Rename JSON file with the name of the new EV. You have to replace `<ev-name>.json` in lines 3 and 209 with your EV, for example `ev-d68.json`.
- You can change the length for sequence filtering (set to 6000 nt, to keep whole-genomes)

## 🐍 Snakemake build
To get the json file, go to the folder of your EV and run:
   ```sh
   snakemake --cores X
   ```
Replacing X by the number of cores you want to use.

This will create two additional folders in the directory:
- 📁 results: Contains several files used for building the json
- 📁 auspice: Contains the json file named `<ev-name>.json`

To remove data from the resulting folders (auspice, results), run:
   ```sh
   snakemake --cores X clean
   ```
Replacing X by the number of cores you want to use.

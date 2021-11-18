## PRECOGx webApp
![action main](https://github.com/raimondilab/precogx/actions/workflows/main.yml/badge.svg)

![Legacy version](https://github.com/gurdeep330/precog)

A Flask-based webApp to predict GPCR/G-protein couplings.

## File organization
```data```: Files/folders requried to run the app

```templates```: HTML templates to load the home and output page

```static```: JS files to load Navigation bar as well as sequence and structure panels

```.github/workflows```: Workflows

```test.py```: Python script (caled in the workflow) to fetch latest SIFT Pfam/Chain files and annotates each PDB-ID with GPCR and G-protein chains

## How to use?
```
git clone <repo>
pip install Flask
python3 run.py
```

Open the URL on your browser ```http://129.206.245.88:5000```

It currently runs using a sample output file ```static/OL820/out.tsv```

## Features enabled
1. Load FASTA (in Sequence panel) and 3D strucrture (Structure panel) of the first entry
2. Update Seqeuence and Structure panels by clicking on <b>GPCR</b>
3. Update PDB list on every click (sort based on sequence identity)

## Workflows
1. The current workflow ```.github/workflows/main.yml``` fetches the latest SIFT mappings
2. Extracts all GPCR/G-protein complexes and save them in ```data/pdb_list.txt```
3. Fetches PDB FASTA files (```data/fasta/```) of the complexes and concatenates them into ```data/fasta/all_pdbs.fasta```
4. Set up Miniconda and insall BLAST
5. Create a blastdb of ```data/fasta/all_pdbs.fasta```, which is required to find the best structures (Structure panel)

Use ![Issues](https://github.com/gurdeep330/precogx/issues) to suggest more features or report problems.

More updates to follow soon....

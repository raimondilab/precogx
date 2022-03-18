## PRECOGx webApp
![action main](https://github.com/raimondilab/precogx/actions/workflows/main.yml/badge.svg)

A Flask-based webApp to visualize PRECOGx predictions. Source code and analysisof PRECOGx can be found ![here](https://github.com/raimondilab/transformers/tree/main/predictor)

PRECOG (legacy version) can be found ![here](https://github.com/gurdeep330/precog)

## File organization
```data```: Files/folders requried to run the app (mostly loaded in the run.py script)

```templates```: HTML templates to load the home, output, result, etc pages

```static```: JS, css, .py files to implement the front-end

```.github/workflows```: Workflow to fetch the latest 3D complexes of GPCR/G-protein and generate mappings/annotations, which are used by various functions in the run.py script

## How to use?
```git clone <repo>
pip install Flask
python3 run.py
```
Open the URL on your browser ```http://129.206.245.88:5000```

# Output page

## Coupling probabilites
1. The predicted probabilites are displayed in the top left panel
2. Each cell is clickable, and every click will call for functions to update other panels

## PCA
1. For every layer, embeddings of 377 GPCRs were generated and their PCA was calculated
2. The PC1 and PC2 are showed in this panel
3. From the dropdown, the user can select the layer as well as the type of annotation

## Features enabled
1. Load FASTA (in Sequence panel) and 3D strucrture (Structure panel) of the first entry
2. Update Seqeuence and Structure panels by clicking on <b>GPCR</b>
3. Update PDB-ID list on every click (sort based on sequence identity)

## Automate workflow
1. The current workflow ```.github/workflows/main.yml``` fetches the latest SIFT mappings
2. Extracts all GPCR/G-protein complexes and save them in ```data/pdb_list.txt```
3. Fetches PDB FASTA files (```data/fasta/```) of the complexes and concatenates them into ```data/fasta/all_pdbs.fasta```
4. Set up Miniconda and insall BLAST
5. Create a blastdb of ```data/fasta/all_pdbs.fasta```, which is required to find the best structures (Structure panel)

Use ![Issues](https://github.com/gurdeep330/precogx/issues) to suggest more features or report problems.

More updates to follow soon....

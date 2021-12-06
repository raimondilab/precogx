### How to run PRECOGx?
1. Clone the ![repo](https://github.com/maticmarin/transformers) or if you have already cloned it once before, pull the latest commit:
```
git pull
```
2. To run the script, provide `input file` and `assay name` as positional arguments:
```
./run.py <input.file> <assay.name>
```
for ex:
```
./run.py sample/input.fasta shed
```
3. Use `--help` to know more
```
./run.py --help
```
4. This will then create a sub-directory with a **unique ID** as its name in `output/`. Predictions will be stored in **out.tsv**.
5. Only **FASTA formatted** files can be given provided as input for now (![example](https://github.com/maticmarin/transformers/blob/main/predictor/sample/input.fasta)). Even if it is a point mutation, a FASTA formatted (mutated) sequence ***must*** be provided. Input GPCRs must not contain `/` and `|`. You may use `_`.

A sample output can be found ![here](https://github.com/maticmarin/transformers/blob/main/predictor/output/9N4QF/out.tsv)

![Want to raise an issue?](https://github.com/maticmarin/transformers/issues)

### Additional information (optional)
1. Specifcy the location where you have saved the models in the variable **path_to_model** of `run.py`.
2. Write the names of the best models and the information if they include the seqeuenced-based (handcrafted) features (i.e. 1 or 0), respectively, in the matrix variable **data** of `run.py`.


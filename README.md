# AMPEvolve

### How to evaluate hill-climbing algorithms using DBAASP:

1. Run hill climbing on some base sequences:
```python .\src\hill_climbing\hill_climber.py --output <output_folder> --input_seqs .\inputs.fasta```
2. Put the resulting fasta file into the DBAASP we tool
3. Copy out results to a tsv file
4. Convert the tsv file to a dbaasp json file:
```python .\src\hill_climbing\tsv_to_json.py <tsv_file> <dbaasp_json_file>```
5. Plot the results: 
```python .src\hill_climbing\plot_hill_climbing.py <dbaasp_json_file_0> <dbaasp_json_file_1> ...```
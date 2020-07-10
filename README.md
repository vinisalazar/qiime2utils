### qiime2utils

Scripts to manipulate Qiime 2 data. Functionality so far includes:
* **filter_by_category:** Group data by category and filter most abundant organisms.
* **extract_and_blast:** Subset ASVs from a FASTA file and BLAST.

```
# Install
pip install qiime2utils

# Test installation
qiime2utils -h

# Example usage
qiime2utils filter_by_category -f feature_table.qza -t taxonomy.qza -m manifest.tsv -c 'Host_Family' -n 30
```

For any problems or suggestions, please open an issue at the [GitHub repository](https://github.com/vinisalazar/qiime2utils). Contributions are welcome!

**qiime2utils is in active development and comes without warranties of any kind (please see the [license file.](LICENSE))**

##### Test data
At the moment, I am unable to provide test data. It will come in future versions.
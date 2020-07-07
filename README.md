### qiime2utils

Scripts to manipulate Qiime 2 data. Functionality so far includes:
* Grouping data by metadata columns
* Filtering by *n* most abundant ASVs

```
# Install
pip install qiime2utils

# Test installation
qiime2utils -h

# Example usage
qiime2utils -f feature_table.qza -t taxonomy.qza -m manifest.tsv -c 'Host_Family' -n 30
```

For any problems or suggestions, please open an Issue at the [GitHub repository](https://github.com/vinisalazar/qiime2utils). Contributions are welcome!

**qiime2utils is in active development and comes without warranties of any kind (please see the [license file](LICENSE) for details.)**
# DeSigNate

The scientific manuscript of this work is currently under revision. Links to the publication will be added shortly.

If you want to access the tool without cloning the respository, we host a publicly available server which runs the web-interface of DeSigNate at https://designate.dbresearch.uni-salzburg.at/.

## Tool description

### Detecting Signature Nucleotides for Taxon Diagnoses
DeSigNate is an innovative tool for detecting diagnostic nucleotide positions for taxon diagnoses. The analysis is based on a novel representation of the gene sequence data, which enables a ranking of all positions according to their diagnostic relevance and a classification of the nucleotides. DeSigNate is also able to detect diagnostic combinations of nucleotides. The tool guides the user step-by-step through the analysis and presents the results without need to post-process the output data.

### Which nucleotides are suitable for taxon diagnoses?
In taxon diagnoses, only nucleotides that unambiguously distinguish the query from the reference group are of interest, i.e., they are uniform at homologous alignment positions in the query group. Two types of signature nucleotides are distinguished:
1. at binary positions the nucleotides of the reference group are uniform but different from the nucleotide in the query group
2. at asymmetric positions the nucleotides of the reference group are not uniform but different from the nucleotide in the query group.

## Usage

### Use DeSigNate via command-line

DeSigNate can be used as a fully functional command-line interface independent of the web-interface. Hence, the tool runs without installing Django and does not use Javascript libraries; however, Python 3 has to be installed.

To see the full list of available command-line parameters, execute the following command within the webapp directory:
```
python3 designate.py -h
```
Using the command-line interface, the query and the reference group are determined based on the preorder ID of an inner node in the phylogenetic tree. The phylogenetic tree can be printed with the following command:
```
python3 designate.py --phylotree <path_to_tree_file> --print_tree
```
In the following example, the query group contains all species of the subtree with the preorder ID 3 (resp. 5 for the reference group). With this configuration, individual and combined signature nucleotides (k-window: k = 5) are identified and the entropy is calculated.
```
python3 designate.py --alignment <path_to_alignment_file> --phylotree <path_to_tree_file> --reference_group_id 5 --query_group_id 3 --signature_nucleotides --two_pos 5 --shannon_entropy
```

### Use DeSigNate via web-interface

Next to the command-line interface, this repository contains all the data needed to run web-interface of this tool locally. In order to do that, Django (Version 2.2.4) must be installed. The local server can be started with the following command in the root directory of DeSigNate:
```
python3 manage.py runserver
```
Now, the tool is available at localhost port 8000 (usually 127.0.0.1:8000).
Be aware that Google Analytics will also track local usage of the web-interface.


## Used libraries
- **ChartJS:** used to display the entropy plot on the results page. (https://www.chartjs.org/)
- **BioPython:** used to process the alignment and tree data. (https://biopython.org/)
- **Phylogenetic tree library:** used to display the phylogenetic tree. (http://bl.ocks.org/jacksonhenry3/3083726)
- **D3:** used by the Phylogenetic tree library to display the phylogenetic tree. (https://d3js.org/)
- **newick.js:** used to parse the newick file. (https://github.com/jasondavies/newick.js)
- **Google Analytics:** used to track the number of users. (https://analytics.google.com/analytics/web/)
- **Cookie Consent:** used to display cookie message. (https://cookieconsent.osano.com/)
## Paper and Citation
If you use DeSigNate for your research, please cite the publication with DOI: >link.
The paper is available at >link.

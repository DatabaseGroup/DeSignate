# DeSigNate

This work is under revision. Links to the publication will be added shortly.

## Tool description

### Detecting Signature Nucleotides for Taxon Diagnoses
DeSigNate is an innovative tool for detecting diagnostic nucleotide positions for taxon diagnoses. The analysis is based on a novel representation of the gene sequence data, which enables a ranking of all positions according to their diagnostic relevance and a classification of the nucleotides. DeSigNate is also able to detect diagnostic combinations of nucleotides. The tool guides the user step-by-step through the analysis and presents the results without need to post-process the output data.

### Which nucleotides are suitable for taxon diagnoses?
In taxon diagnoses, only nucleotides that unambiguously distinguish the query from the reference group are of interest, i.e., they are uniform at homologous alignment positions in the query group. Two types of signature nucleotides are distinguished:
1. at binary positions the nucleotides of the reference group are uniform but different from the nucleotide in the query group
2. at asymmetric positions the nucleotides of the reference group are not uniform but different from the nucleotide in the query group.

## Usage

### Use DeSigNate via command-line

DeSigNate can be used as a fully functional standalone command-line program without the web-interface or django server.
To use it via command-line, Python3 has to be installed.  
To see the full list of available parameters, execute the following command in the webapp directory:

```
python3 designate.py -h
```

In the following example, the query group contains all species of the subtree with the preorder id 3 (resp. 5 for the reference group). With this configuration, individual and combined signature  nucleotides (k-window: k = 5) are identified and the entropy is calculated.

```
python3 designate.py --alignment <path_to_alignment_file> --phylotree <path_to_tree_file> --reference_group_id 5 --query_group_id 3 --signature_nucleotides --two_pos 5 --shannon_entropy
```

### Use DeSigNate via web-interface

We host a publicly available server which runs DeSigNate at  https://designate.dbresearch.uni-salzburg.at/.
However, the web-interface of DeSigNate can be run locally too.
The tool is implemented using Django (Version 2.2.4). If Django is installed at your machine and you cloned the repository, the local server can be started by using the following command in the root directory of DeSigNate:
```
python3 manage.py runserver
```
The tool will be available at localhost port 8000 (usually 127.0.0.1:8000).
However, google analytics will even track local usage of the web-interface.


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

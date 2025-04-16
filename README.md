# DeSignate


## Tool description

### Detecting Signature Characters for Taxon Diagnoses
DeSignate is an innovative tool for detecting diagnostic character positions for taxon diagnoses. The analysis is based on a novel representation of the gene sequence data, which enables a ranking of all positions according to their diagnostic relevance and a classification of the characters. DeSignate is also able to detect diagnostic combinations of characters. The tool guides the user step-by-step through the analysis and presents the results without need to post-process the output data.

### Which characters are suitable for taxon diagnoses?
In taxon diagnoses, only characters that unambiguously distinguish the query from the reference group are of interest, i.e., they are uniform at homologous alignment positions in the query group. Two types of signature characters are distinguished:
1. at binary positions, the characters of the reference group are uniform but different from the character in the query group.
2. at asymmetric positions, the characters of the reference group are not uniform but different from the character in the query group.


If you want to access the tool without cloning the repository, we host a publicly available server which runs the web-interface of DeSignate at https://designate.dbresearch.uni-salzburg.at/.


## Paper and Citation
The scientific manuscript of this tool was published at BMC Bioinformatics. The paper is open access and can be found at https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3498-6. Please cite this article as described on the bottom of the BMC Bioinformatics website.


## Usage

### Use DeSignate via command-line

DeSignate can be used as a fully functional command-line interface independent of the web-interface. Hence, the tool runs without installing Django and does not use Javascript libraries; however, Python 3 has to be installed.

To see the full list of available command-line parameters, execute the following command within the webapp directory:
```
python3 designate.py -h
```
Using the command-line interface, the query and the reference group are determined based on the preorder ID of an inner node in the phylogenetic tree. The phylogenetic tree can be printed with the following command:
```
python3 designate.py --phylotree <path_to_tree_file> --print_tree
```
In the following example, the query group contains all species of the subtree with the preorder ID 3 (resp. 5 for the reference group). With this configuration, individual and combined signature characters (k-window: k = 5) are identified and the entropy is calculated.
```
python3 designate.py --alignment <path_to_alignment_file> --phylotree <path_to_tree_file> --reference_group_id 5 --query_group_id 3 --signature_characters --k_window 5 --shannon_entropy
```

### Use DeSignate via web-interface

Next to the command-line interface, this repository contains all the data needed to run web-interface of this tool locally. In order to do that, Django (Version 2.2.4) must be installed. The required javascript libraries can be installed via the following command in the /webapp subdirectory:
```
npm install
```
After installing the javascript libraries, the local server can be started with the following command in the root directory of DeSignate:
```
python3 manage.py runserver
```
Now, the tool is available at localhost port 8000 (usually 127.0.0.1:8000).


## Used libraries
- **ChartJS:** used to display the entropy plot on the results page. (https://www.chartjs.org/)
- **BioPython:** used to process the alignment and tree data. (https://biopython.org/)
- **Phylogenetic tree library:** used to display the phylogenetic tree. (http://bl.ocks.org/jacksonhenry3/3083726)
- **D3:** used by the Phylogenetic tree library to display the phylogenetic tree. (https://d3js.org/)
- **newick.js:** used to parse the newick file. (https://github.com/jasondavies/newick.js)
- **Cookie Consent:** used to display cookie message. (https://cookieconsent.osano.com/)

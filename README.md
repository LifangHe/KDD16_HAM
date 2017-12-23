# A Python implementation of HAM algorithm

Please see the following paper for more details of HAM algorithm:
* He, L., Lu, C. T., Ma, J., Cao, J., Shen, L., & Philip, S. Y. Joint Community and Structural Hole Spanner Detection via Harmonic Modularity. Proceedings of the 22th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD), 2016.


## *** USAGE
To try the code, we provide a graph benchmark dataset “Karate” as an example, which includes adjacency matrix and ground-truth community (label) interface for graph files and code.

### karatern.txt - adjacency list
The first line is a summary of the graph (# Nodes || # Edges) (space-delimited)

### karatecrn.txt - ground-truth community
Each row is a COMMUNITY containing the indices of nodes (tab-delimited)

`—------—------—------—------NOTE—-----—-------—------—------`

Please change the following code according to your environment and applications.

datapath= os.getcwd() + '/data/'  <br />
adj_name = 'karatern'  &emsp; &emsp;   # Note: the data index starts from 0, not from 1  <br />
community_name = 'karatecrn'  <br />
topk=3     &emsp; &emsp;      # the number of selected SH spanner

*See comments inline of the code for more details.*

### In terminal
You can use the following routine to perform the code.

$ python HAM.py

*For SHII metric, please see `README` file in the `SHII_metric` folder.*

## Contact
Lifang He (Email: lifanghescut@gmail.com)

# A C++ implementation of SHII metric


**—------—------—------—------NOTE—-----—-------—------—------**

The computation of SHII metric is implemented based on the use of CELF++ code, which can be downloaded from http://www.cs.ubc.ca/~goyal/code-release.php. We followed the default setting of this code, and it can be tuned by changing the weight values of each edge. In the config file, you can change the parameter "edgeWeight" to control it. In our experiment, we set it to 0.1 on all datasets.

Specifically, SHII needs to set the seedRatio.

## *** USAGE
To try the code, we provide a graph benchmark dataset “Karate” as an example, which can be found in the “datasets” file. Note that the index of selected top-k SHS starts from 0, not from 1.

See comments inline of the code for more details.

### In terminal
You can use the following routine to perform the code.

$ make  <br />
$ ./InfluenceModels -c config_test.txt

### Configuration
In the config file, one can specify various parameters like input file for network, propagation model etc. These options can also be specified on the command line. If a parameter is present in both command line and config file, the command line has the preference. The code is written in MC.cc file. This module selects the seed set under LT or IC model using Monte Carlo Simulations. MC stands for Monte Carlo. 

Parameters needed:

1. propModel : Propagation Model. Should be IC or LT.
2. probGraphFile : File containing the edges of the network in Matrix Market format. Each line should contain at least 3 columns, separated by a space. First column is user 1 id, second column is user 2 id and the third column is the edge weight of user 1 on user 2. The first three lines of the file are ignored by the code, and can be used for some comments.
3. labelFile : File containing the labels of the each node. Each line represents the nodes within one community.
4. seedFile : File containing the node ID of the selected top-k structural holes spanners (SHSs).
5. datasetFile : Output file name, which will be automatically appended with the model name
6. randSeedIter : How many iterations to sample seeds (which contains a SHS, and its sampled neighbors). Default value is 100.
7. seedRatio : # of sampled seeds / # of nodes of the community that the given SHS belongs to. Default value is 0.05. For a small dataset, e.g., Karate dataset, one should set the value appropriately to guarantee that at least one neighbor of the SHS can be sampled. We set the value of seedRatio to be 0.1 for Karate dataset in our experiment.
8. edgeWeight: the weight to be assigned to each edge. Default value is 0.1.
9. mcruns : Number of monte carlo simulations to be used. Default value is 10000.
10. outdir : Output directory where the output should be written.

### Output
The output is created in the `<outdir>`. The outfile filename is of the form `datasetFile_<model>.txt`. <br />
Each row represents a seed node, with two columns  <br />
The first column represents the SHII metric: # of the influenced outsiders / # of the total influenced nodes  <br />
The second column represents a variant of the SHII metric: # of the influenced outsider / # of the influenced insiders

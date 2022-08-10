# MotifDiscovery-Algos
Integration of some of the best existing Motif discovery algorithms in time series.

# Introduction
This library integrates a number of time series motif discovery algorithms and incorporates pattern discovery, visualisation and motifs storage.

The mdalgos folder contains algos, i.e. algorithmic ensembles. These include STOMP, STUMP, GPU_STUMP, SCRIMP++, GRAMMARINTRODUCTION and MK algorithms.

The Test folder contains the data Euro_pfund.csv for testing and the trial code written in jupyter lab. Data source: [European Central Bank, Frankfurt am Main, Germany](https://sdw.ecb.europa.eu/browse.do?node=1495). Euro_Pfund.csv is the data that has been processed by the author.

Download all files as a zip file and install the zip file via pip.
But please install the necessary libraries mentioned in Requiements before installing this library. Otherwise it will not work.

The library can be imported at the time of use in the following way：

```python
from mdalgos.algos import *
```
The algorithm can be used in subsequent programs by entering the name of the algorithm, see **Example** for details.

Two main parameters are required when calling the corresponding algorithm.

One is the time series data, which must be univariate, i.e. similar to Euro_Pfund.csv in the Test folder, and in particular the time column index must be named "t".
The other is the storage path, in the following format: **path='/home/'**, with a "/" at the end. if you don't enter a path, it will default to the folder where the current application is located.

After observing the visualisation results, you can go to the specified storage location and see the motif data files, each with a motif data file name in the format: "Parameter name + Motif order + Algorithm name". Therefore, when making subsequent comparisons, it is possible to identify them directly by name.

When a comparison between different motifs is required, the similarity comparison algorithm can be invoked directly by importing the data from two motifs, and the scores allow direct confirmation of their similarity. This section can be used to compare, for example, the motifs of the same component from different manufacturers to distinguish the differences between their tuning of the component's motion characteristics.

# Background
The discovery of Motif in time series has large implications for the study of time series data. However, there are numerous existing algorithms with different principles and no one has integrated them.

By integrating the different algorithms, it will be easier to use different algorithms for motif discovery in time series and to compare the algorithms in order to choose the best one for motif discovery.

The algorithms chosen here are MK, STOMP, STUMP, SCRIMP++, Grammarintroduction. Based on raw sequence violence search, Matrirx profile and SAX.

# Advantages

This algorithm package has the following Advantages:

1. The call of different algorithm is simple and straightforward
2. Direct visualisation of results for easy comparison
3. Motif can be automatically saved as a csv file and exported

# limitations:
1. The input is restricted to univariate time series, so if the original data is a multivariate time series, the data needs to be extracted separately for motif discovery.
2. The algorithms used are not necessarily optimal, as some of them do not yet have good conversion code in python.

# Requirements
pandas

numpy

matplotlib

[RePair](https://github.com/axelroques/RePair)

[SAX](https://github.com/axelroques/SAX)

[Grammar Induction](https://github.com/axelroques/GrammarInduction)

[mkalgo](https://github.com/saifuddin778/mkalgo#mkalgo-mk-algorithm) 

[STUMPY](https://github.com/TDAmeritrade/stumpy)

[Matrix Profile](https://github.com/matrix-profile-foundation/matrixprofile)

# Example
## Data read
__ts__: 
1.One column is time (__must be named as 't'__)， One column for data. 
2.If all values in a column are identical, then the timeseries is not available.(Will cause a running error)

```python
ts = pd.read_csv('Euro_pfund.csv')
```

![Image text](https://github.com/7SiebenPunch/img-folder/blob/main/Testdata.png)
## Import parameters
__ts__: Time series data, in the format described in **Data read**, must meet the requirements, otherwise the algorithms in the library cannot be used.

__path__: The path to export the motif as a csv file for storage. 
          The format is: '/home/', if you don't enter a path, it will default to the folder where the current application is located.
          If you don't want to export the file, delete the csv export code in algo.py

__windowsize:__ The time window, in this case, can be thought of as the length of the motif.

__topk:__ The first k motifs to be extracted and visualised.

Note: The MK algorithm can only extract one motif, so the topk parameter is not required. In addition, the distance calculation method 'euclidean' or 'dtw' can be selected in the MK algorithm.
      
## MK(ts, windowsize = 150, metrix = 'euclidean')
The MK algorithm is based on:
Mueen, A., Keogh, E., Zhu, Q., Cash, S. and Westover, B. [Exact Discovery of Time Series Motif](http://alumni.cs.ucr.edu/~mueen/pdf/EM.pdf)
. SDM 2009

When using the MK algorithm, the distance calculation method can be chosen from "euclidean" and "dtw"


```python
from mdalgos.algos import *

mkal(ts, path='/home/')
```
## STOMP, SCRIMP++, STUMP, GPU_STUMP(ts, windowsize=200, topk=3)
These algorithms are based on：
1. Zhu Y, Zimmerman Z, Senobari NS, Yeh CC, Funning G, Mueen A, Brisk P, Keogh E. [Matrix profile ii: Exploiting a novel algorithm and gpus to break the one hundred million barrier for time series motifs and joins](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7837898). In2016 IEEE 16th international conference on data mining (ICDM) 2016 Dec 12 (pp. 739-748). IEEE.]
2. Zhu Y, Yeh CC, Zimmerman Z, Kamgar K, Keogh E. [Matrix profile XI: SCRIMP++: time series motif discovery at interactive speeds](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8594908). In2018 IEEE International Conference on Data Mining (ICDM) 2018 Nov 17 (pp. 837-846). IEEE.
3. Yeh CC, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, Silva DF, Mueen A, Keogh E. [Matrix profile I: all pairs similarity joins for time series: a unifying view that includes motifs, discords and shapelets](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7837992). In2016 IEEE 16th international conference on data mining (ICDM) 2016 Dec 12 (pp. 1317-1322). Ieee.

```python
stomp(ts, path='/home/')
scrimp_2plus(ts, path='/home/')
stump(ts, path='/home/')
gpu_stump(ts, path='/home/')
```
**Note:** For a multivariate time series, when using traversal to perform motif discovery on individual variables in turn
stomp and scrimp_2plus will run with an error if a column has exactly the same data (i.e. its data is a straight line).

## Grammarintroduction(ts, w=6, n=100, k=10, topk=3) 
Based on:
1. Senin, P., Lin, J., Wang, X., Oates, T., Gandhi, S., Boedihardjo, A.P., Chen, C., Frankenstein, S., Lerner, M., [GrammarViz 2.0: a tool for grammar-based pattern discovery in time series](http://www2.hawaii.edu/~senin/assets/papers/grammarviz2.pdf), ECML/PKDD Conference, 2014.
2. Larsson NJ, Moffat A. [Off-line dictionary-based compression](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=892708). Proceedings of the IEEE. 2000 Nov;88(11):1722-32.

__w:__ Subsequence character length

__n:__ subsequences length

__k:__ stride of k

```python
grammarintroduction(ts, path='/home/')
```

# License
[MIT LICENSE](https://github.com/7SiebenPunch/MotifDiscovery-Algos/blob/main/LICENSE)

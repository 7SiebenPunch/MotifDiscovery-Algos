# MotifDiscovery-Algos
Integration of some of the best existing Motif discovery algorithms in time series.

# Introduction


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

[STUMPY](https://github.com/TDAmeritrade/stumpy）

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
__ts__: ts

__path__: The path to export the motif as a csv file for storage. 
          The format is: '/home/', if you don't enter a path, it will default to the folder where the current application is located.
          If you don't want to export the file, delete the csv export code in algo.py

__windowsize:__ The time window, in this case, can be thought of as the length of the motif.

__topk:__ The first k motifs to be extracted and visualised.

Note: The MK algorithm can only extract one motif, so the topk parameter is not required. In addition, the distance calculation method 'euclidean' or 'dtw' can be selected in the MK algorithm.
      
## MK(ts, windowsize = 150, metrix = 'euclidean')

```python
from mdalgos import algos

mkal(ts, path='/home/')
```
## STOMP, SCRIMP++, STUMP, GPU_STUMP(ts, windowsize=200, topk=3)

```python
stomp(ts, path='/home/')
scrimp_2plus(ts, path='/home/')
stump(ts, path='/home/')
spu_stump(ts, path='/home/')
```

## Grammarintroduction(ts, w=6, n=100, k=10, topk=3) 
__w:__ Subsequence character length

__n:__ subsequences length

__k:__ stride of k

```python
grammarintroduction(ts, path='/home/')
```

# License
[MIT LICENSE](https://github.com/7SiebenPunch/MotifDiscovery-Algos/blob/main/LICENSE)

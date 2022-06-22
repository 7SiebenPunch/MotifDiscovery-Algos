# MotifDiscovery-Algos
Integration of some of the best existing Motif discovery algorithms in time series.

# Background
The discovery of Motif in time series has large implications for the study of time series data.However, there are numerous existing algorithms with different principles and no one has integrated them.

By integrating the different algorithms, it will be easier to use different algorithms for motif discovery in time series and to compare the algorithms in order to choose the best one for motif discovery.

The algorithms chosen here are MK, STOMP, STUMP, SCRIMP++, Grammarintroduction.
Based on raw sequence violence search, Matrirx profile and SAX.

# Advantages

This algorithm package has the following Advantages:

1. The call of different algorithm is simple and straightforward
2. Direct visualisation of results for easy comparison

# limitations:
1. The input is restricted to univariate time series, so if the original data is a multivariate time series, the data needs to be extracted separately for motif discovery.
2. The algorithms used are not necessarily optimal, as some of them do not yet have good conversion code in python.

# Requirements
pandas

numpy

matplotlib

RePair (https://github.com/axelroques/RePair)

SAX (https://github.com/axelroques/SAX)

Grammar Induction (https://github.com/axelroques/GrammarInduction)

mkalgo ([MK Algorithm](https://github.com/saifuddin778/mkalgo#mkalgo-mk-algorithm)) 

STUMPY (https://github.com/TDAmeritrade/stumpy）

Matrix Profile (https://github.com/matrix-profile-foundation/matrixprofile)

# Example
## Data read
```python
ts = pd.read_csv('Euro_pfund.csv')
```
## Import

ts: One column is time (must be named as 't')， One column for data. 

![Image text](https://github.com/7SiebenPunch/img-folder/blob/main/Testdata.png)

windowsize: The time window, in this case, can be thought of as the length of the motif.

topk: The first k motifs to be extracted and visualised.

Note: The MK algorithm can only extract one motif, so the topk parameter is not required. In addition, the distance calculation method 'euclidean' or 'dtw' can be selected in the MK algorithm.
      
## MK(ts, windowsize = 150, metrix = 'euclidean')

```python
from mdalgos import algos

mkal(ts)
```
## STOMP, SCRIMP++, STUMP, GPU_STUMP(ts, windowsize=200, topk=3)

```python
stomp(ts)
scrimp_2plus(ts)
stump(ts)
spu_stump(ts)
```

## Grammarintroduction(ts, df, w=6, n=100, k=10, topk=3) 
w:Subsequence character length

n:subsequences length

k:stride of k

```python
grammarintroduction(ts)
```

# License
MIT License

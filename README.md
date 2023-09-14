# TS-clustering-using-graph

In this project, a Graph based Time series data mining framework (Graft) for time series (TS) dataset, is shown.  The TS dataset consists of multiple TS samples from varying domain like- engineering and utilities. 

To obtain a single weighted directed graph structure for entire TS dataset,  the Piecewise aggregate approximation (PAA) techniques is used that gives multiple non-overlapping fixed length subsequences. 
The features extracted from the subsequences are converted into symbols using Symbolic aggregate approximation technique. The weighted directed graph structure help to represent the TS into a unique structure 
of nodes and edges that can capture the temporal nature of co-occurring patterns. 

The whole TS clustering approach is proposed based in the features obtained from the graph paths. The graph path based clustering can cluster the TS data based on frequently appearing longest temporal patterns.
A simple yet effective graph based approach for identification of temporal dependent rare events is also shown.

More information about the work and the algorithms are discussed in the given link  https://www.sciencedirect.com/science/article/pii/S0952197622000215

The dataset details are given in:

1. London datastore-- https://data.london.gov.uk/dataset/smartmeter-energy-use-data-in-london-households
2. Ausgrid solar homes-- https://www.ausgrid.com.au/Industry/Our-Research/Data-to-share/Solar-home-electricity-data
3. Stock Market-- https://www.kaggle.com/datasets/rohanrao/nifty50-stock-market-data
4. Web Traffic history-- https://www.kaggle.com/c/web-traffic-time-series-forecasting

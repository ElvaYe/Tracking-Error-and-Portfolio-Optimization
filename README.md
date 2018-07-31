# Tracking-Error-and-Portfolio-Optimization
The input data file has to be the same format as the given example file. The first column of the return matrix is the index of rows, starting at 0. The second column is the corresponding date. Then, each column is for each asset's daily return.

When using the code, change the "startRow" and "endRow" to the row with the date you want to calculate. (Notice the oldest date starts from the bottom.) Also, change the integer "n" to the number of assets in the portfolio.

The output will be stored into a new excel file, including weights of all assets under tracking error minimization as well as MPT. Furthermore, it also gives the value of tracking error and portfolio volatility for a fixed vector of expected (excess) return, so that efficient frontior can be easily generated in excel.

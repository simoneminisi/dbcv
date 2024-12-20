# DBCV
R implementation of Density-Based Clustering Validation based on the [DBCV Python module](https://github.com/christopherjenness/DBCV).

## Source

Moulavi, Davoud, et al. "Density-based clustering validation." Proceedings of the 2014 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2014.

[PDF](http://epubs.siam.org/doi/pdf/10.1137/1.9781611973440.96)

## What is DBCV

How do you validate clustering assignmnets from unsupervised learning algorithms?  A common method is the [Silhoette Method](https://en.wikipedia.org/wiki/Silhouette_(clustering)), which provides an objective score between -1 and 1 on the quality of clustering.  The silhouette value measures how well an object is classified in its own cluster instead of neighboring clusters.  The silhouette (and most other popular methods) work very well on globular clusters, but can fail on non-glubular clusters such as:

![non-globular](http://hdbscan.readthedocs.io/en/latest/_images/advanced_hdbscan_5_1.png)

Here, we implement DBCV which can validate clustering assignments on non-globular, arbitrarily shaped clusters (such as the example above).  In essence, DBCV computes two values:

* The density **within** a cluster
* The density **between** clusters

High density within a cluster, and low density between clusters indicates good clustering assignments.

## Example

To assess the quality of clustering, using Density-Based Clustering Validation, we call `dbcv`.

```r
# Install devtools if not already installed
install.packages("devtools")
devtools::install_github("simoneminisi/dbcv")

library(dbcv)

# Create sample data
data <- matrix(rnorm(100), ncol = 2)
labels <- sample(1:3, 50, replace = TRUE)

# Calculate DBCV index
result <- dbcv(data, labels)
print(result)
```  

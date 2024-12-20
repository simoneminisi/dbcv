#' Density Based Clustering Validation
#' 
#' Implementation of DBCV (Density-Based Clustering Validation)
#' as described in Moulavi et al. (2014)
#' 
#' @param X Matrix with dimensions [n_samples, n_features]
#' @param labels Vector of cluster assignments
#' @param dist_function Function to calculate distance between points (default: euclidean)
#' @return Cluster validity score in range [-1, 1]
#' @export
DBCV <- function(X, labels, dist_function = function(x, y) sqrt(sum((x - y)^2))) {
  graph <- mutual_reach_dist_graph(X, labels, dist_function)
  mst <- mutual_reach_dist_MST(graph)
  cluster_validity <- clustering_validity_index(mst, labels)
  return(cluster_validity)
}

cdist <- function(XA, XB) {
    if(is.null(dim(XA))) XA <- matrix(XA, nrow=1)
    if(is.null(dim(XB))) XB <- matrix(XB, nrow=1)
    
    na <- nrow(XA)
    nb <- nrow(XB)
    
    result <- matrix(0, na, nb)
    
    for(i in 1:na) {
        for(j in 1:nb) {
            result[i,j] <- sqrt(sum((XA[i,] - XB[j,])^2))
        }
    }
    return(result)
}

#' Calculate core distance of a point
#' 
#' @param point Vector of length n_features
#' @param neighbors Matrix of dimensions [n_neighbors, n_features]
#' @param dist_function Distance function
#' @return Core distance (float)
core_dist <- function(point, neighbors, dist_function) {
  n_features <- length(point)
  n_neighbors <- nrow(neighbors)
  
  # Calculate distances excluding self
  distance_vector <- cdist(matrix(point, nrow=1), neighbors)
  distance_vector <- distance_vector[distance_vector != 0]
  
  # Calculate core distance using density formula
  numerator <- sum((1/distance_vector)^n_features)
  core_dist <- (numerator / (n_neighbors - 1))^(-1/n_features)
  return(core_dist)
}

#' Calculate mutual reachability distance between two points
#' 
#' @param point_i Vector representing point i
#' @param point_j Vector representing point j
#' @param neighbors_i Matrix of neighbors for point i
#' @param neighbors_j Matrix of neighbors for point j
#' @param dist_function Distance function
#' @return Mutual reachability distance
mutual_reachability_dist <- function(point_i, point_j, neighbors_i, neighbors_j, dist_function) {
  core_dist_i <- core_dist(point_i, neighbors_i, dist_function)
  core_dist_j <- core_dist(point_j, neighbors_j, dist_function)
  dist <- dist_function(point_i, point_j)
  mutual_reachability <- max(c(core_dist_i, core_dist_j, dist))
  return(mutual_reachability)
}

#' Create mutual reach distance graph
#' 
#' @param X Data matrix
#' @param labels Cluster assignments
#' @param dist_function Distance function
#' @return Matrix of mutual reachability distances
mutual_reach_dist_graph <- function(X, labels, dist_function) {
  n_samples <- nrow(X)
  graph <- matrix(0, n_samples, n_samples)
  
  for (row in 1:n_samples) {
    for (col in 1:n_samples) {
      point_i <- X[row,]
      point_j <- X[col,]
      class_i <- labels[row]
      class_j <- labels[col]
      members_i <- get_label_members(X, labels, class_i)
      members_j <- get_label_members(X, labels, class_j)
      
      dist <- mutual_reachability_dist(point_i, point_j,
                                     members_i, members_j,
                                     dist_function)
      graph[row, col] <- dist
    }
  }
  return(graph)
}

#' Calculate minimum spanning tree from distance matrix
#' 
#' @param dist_tree Distance matrix
#' @return Symmetric MST matrix
mutual_reach_dist_MST <- function(dist_matrix) {
  n <- nrow(dist_matrix)
    
    # Convert to edge list format (similar to scipy's csr_matrix)
    edges <- which(lower.tri(dist_matrix) & dist_matrix > 0, arr.ind = TRUE)
    weights <- dist_matrix[edges]
    
    # Sort edges by weight (Kruskal's algorithm)
    sorted_idx <- order(weights)
    edges <- edges[sorted_idx,]
    weights <- weights[sorted_idx]
    
    # Union-find data structure
    parent <- 1:n
    rank <- rep(0, n)
    
    # Find set with path compression
    find <- function(x) {
        if (parent[x] != x) {
            parent[x] <<- find(parent[x])
        }
        parent[x]
    }
    
    # Union by rank
    union <- function(x, y) {
        px <- find(x)
        py <- find(y)
        
        if (px == py) return(FALSE)
        
        if (rank[px] < rank[py]) {
            parent[px] <<- py
        } else if (rank[px] > rank[py]) {
            parent[py] <<- px
        } else {
            parent[py] <<- px
            rank[px] <<- rank[px] + 1
        }
        TRUE
    }
    
    # Build MST
    mst <- matrix(0, n, n)
    edges_added <- 0
    
    for (i in seq_along(weights)) {
        u <- edges[i, 1]
        v <- edges[i, 2]
        
        if (union(u, v)) {
            mst[u, v] <- weights[i]
            edges_added <- edges_added + 1
            if (edges_added == n - 1) break
        }
    }
    
    # Make symmetric
    mst <- mst + t(mst)
    return(mst)
}

#' Calculate cluster density sparseness
#' 
#' @param MST Minimum spanning tree matrix
#' @param labels Cluster assignments
#' @param cluster Cluster ID
#' @return Cluster density sparseness
cluster_density_sparseness <- function(MST, labels, cluster) {
  indices <- which(labels == cluster)
  cluster_MST <- MST[indices, indices]
  cluster_density_sparseness <- max(cluster_MST)
  return(cluster_density_sparseness)
}

#' Calculate density separation between clusters
#' 
#' @param MST Minimum spanning tree matrix
#' @param labels Cluster assignments
#' @param cluster_i First cluster ID
#' @param cluster_j Second cluster ID
#' @return Density separation value
cluster_density_separation <- function(MST, labels, cluster_i, cluster_j) {
  indices_i <- which(labels == cluster_i)
    indices_j <- which(labels == cluster_j)
    
    n <- nrow(MST)
    distances <- matrix(Inf, nrow=n, ncol=n)
    
    # Dijkstra's algorithm implementation to match scipy
    dijkstra <- function(start) {
        dist <- rep(Inf, n)
        dist[start] <- 0
        visited <- rep(FALSE, n)
        
        for(i in 1:n) {
            # Find min distance vertex not yet processed
            u <- which.min(dist * (!visited))
            visited[u] <- TRUE
            
            # Update dist value of adjacent vertices
            for(v in 1:n) {
                if(!visited[v] && MST[u,v] > 0 && 
                   dist[u] != Inf && 
                   dist[u] + MST[u,v] < dist[v]) {
                    dist[v] <- dist[u] + MST[u,v]
                }
            }
        }
        return(dist)
    }
    
    # Calculate shortest paths from all points in cluster_i
    all_paths <- t(sapply(indices_i, dijkstra))
    
    # Get relevant paths (to cluster_j points)
    relevant_paths <- all_paths[, indices_j]
    
    # Return minimum path length
    return(min(relevant_paths))
}

#' Calculate validity index for a single cluster
#' 
#' @param MST Minimum spanning tree matrix
#' @param labels Cluster assignments
#' @param cluster Cluster ID
#' @return Cluster validity index
cluster_validity_index <- function(MST, labels, cluster) {
  min_density_separation <- Inf
  unique_labels <- unique(labels)
  
  for (cluster_j in unique_labels) {
    if (cluster_j != cluster) {
      cluster_density_sep <- cluster_density_separation(MST, labels, cluster, cluster_j)
      if (cluster_density_sep < min_density_separation) {
        min_density_separation <- cluster_density_sep
      }
    }
  }
  
  cluster_density_sparse <- cluster_density_sparseness(MST, labels, cluster)
  numerator <- min_density_separation - cluster_density_sparse
  denominator <- max(c(min_density_separation, cluster_density_sparse))
  cluster_validity <- numerator / denominator
  return(cluster_validity)
}

#' Calculate overall clustering validity index
#' 
#' @param MST Minimum spanning tree matrix
#' @param labels Cluster assignments
#' @return Overall validity index
clustering_validity_index <- function(MST, labels) {
  n_samples <- length(labels)
  validity_index <- 0
  unique_labels <- unique(labels)
  
  for (label in unique_labels) {
    fraction <- sum(labels == label) / n_samples
    cluster_validity <- cluster_validity_index(MST, labels, label)
    validity_index <- validity_index + fraction * cluster_validity
  }
  
  return(validity_index)
}

#' Get members of a specific cluster
#' 
#' @param X Data matrix
#' @param labels Cluster assignments
#' @param cluster Cluster ID
#' @return Matrix of cluster members
get_label_members <- function(X, labels, cluster) {
  indices <- which(labels == cluster)
  members <- X[indices, , drop = FALSE]
  return(members)
}
# HELPER
# TRANSFORM M VALUES TO BETA VALUES

# GSE89253 and GSE50759 need transformation from M-values to Beta - done

checkBeta <- function(dataset_path, write=TRUE) {
  file_path <- paste0(beta_path, dataset_path)
  file_ext <- file_ext(file_path)
  if (write==TRUE) {
    if (file_ext=="csv") {
      temp <- fread(file_path) 
      dataset_name <- gsub("_beta.csv", "", dataset_path)
    } else {
      dataset_name <- gsub("_beta.xlsx", "", dataset_path)
      sheet_names <- excel_sheets(file_path)
      temp <- lapply(sheet_names, function(sheet) {
        read_excel(file_path, sheet=sheet)
      }) %>% bind_rows()
    }
    temp <- as.data.frame(temp)
    rownames(temp) <- temp[, 1]
    temp <- temp[, -1]
    btom <- function(m) {
      2^m/(1+2^m)
    }
    temp_t <- apply(temp, c(1, 2), btom)
    temp_t <- as.data.frame(temp_t)
    fwrite(temp_t, paste0(dataset_name, "_beta.csv"), row.names=TRUE)
    return(NULL)
  } else {
    if (file_ext == "csv") {
      temp <- fread(file_path, nrows = 100) 
      dataset_name <- gsub("_beta.csv", "", dataset_path)
    } else {
      dataset_name <- gsub("_beta.xlsx", "", dataset_path)
      sheet_names <- excel_sheets(file_path)
      temp <- lapply(sheet_names, function(sheet) {
        read_excel(file_path, sheet = sheet, range = cell_limits(c(1, 1), c(100, 10)))
      }) %>% bind_rows()
    }
    return(temp)
  }
}

# Check list of all beta files
#list.files(path=beta_path)

# Run
#check <- checkBeta("GSE137495_beta.xlsx", write=FALSE)

# HELPER
# mLiftOver
# Use mLiftOver to harmonize data from GSE114134 and GSE124076 - done,
# original files overwritten

temp <- fread(paste0(beta_path, "GSE124076_beta.csv"))
temp <- as.data.frame(temp)
rownames(temp) <- temp[, 1]
temp <- temp[, -1]
temp2 <- as.matrix(temp)
rownames(temp2) <- rownames(temp)

# Lift over betas
betas <- mLiftOver(temp2, "HM450", impute=F)

# Lift over site names
#sites <- rownames(betas)
#new_sites <- convertProbeID(sites, "HM450")

df_long <- melt(temp2)
colnames(df_long) <- c("Sample", "BetaValue")
ggplot(df_long, aes(x = BetaValue)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Original", x = "Beta Value", y = "Count") +
  theme_minimal()

df_long <- melt(betas)
colnames(df_long) <- c("Sample", "BetaValue")
ggplot(df_long, aes(x = BetaValue)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Adjusted", x = "Beta Value", y = "Count") +
  theme_minimal()


# Figure out how to liftover site names (if necessary)
# Finish liftovers

# HELPER
# Print metadata
# Print Dataset mean and SD ages - for description of datasets table
metadata_file_paths <- list.files(path="metadata18")

for (i in 1:length(metadata_file_paths)) {
  file_ext <- file_ext(metadata_file_paths[i])
  if (file_ext=="csv") {
    data <- read.csv(paste0("metadata18/", metadata_file_paths[i]))
    temp <- data[, c("age")]
    temp <- na.omit(temp)
    temp <- as.numeric(temp)
    temp <- na.omit(temp)
    
    print(paste0("Dataset: ", metadata_file_paths[i], " Max Age: ", max(temp), " Min Age: ", min(temp)))
  } else {
    data <- read_excel(paste0("metadata18/", metadata_file_paths[i]))
    temp <- data[, c("age")][[1]]
    temp <- na.omit(temp)
    temp <- as.numeric(temp)
    temp <- na.omit(temp)
    print(paste0("Dataset: ", metadata_file_paths[i], " Max Age: ", max(temp), " Min Age: ", min(temp)))
  }
}

# 
# for (i in 1:length(metadata_file_paths)) {
#   file_ext <- file_ext(metadata_file_paths[i])
#   if (file_ext=="csv") {
#     data <- read.csv(paste0("metadata18/", metadata_file_paths[i]))
#     temp <- data[, c("gender")]
#     print(paste0("Dataset: ", metadata_file_paths[i]))
#     print(table(temp, useNA="ifany"))
#   } else {
#     data <- read_excel(paste0("metadata18/", metadata_file_paths[i]))
#     temp <- data[, c("gender")][[1]]
#     print(paste0("Dataset: ", metadata_file_paths[i]))
#     print(table(temp, useNA="ifany"))
#   }
# }



# EVERYTHING BELOW IS OLD CODE

\*\* check which age windows to investigate after generating the
histogram

```{r}
filterCorrelationsAndChange <- function(num_splits, training=TRUE) {
  subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/combined_split3/"
  #subset_folder <- testing_subset_folder
  #if (training) {
  #  subset_folder <- training_subset_folder
  #}
  #initialize output list 
  results_by_window <- list()
  for (start_age in seq(5, 70, by=5)) {
    end_age <- start_age + 20
    if (end_age <= 100) {
      results_by_window[[paste(start_age, "to", end_age)]] <- vector("list")
    }
  }
  
  for (j in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_" , j, ".csv")) 
    temp <- as.data.frame(temp)
    #if (!training) {
    #  rownames(temp) <- temp$cpg...207 #466 for training, 207 for testing
    #}
    #if (training) {
    #  rownames(temp) <- temp$cpg...466 #466 for training, 207 for testing
    #}
    rownames(temp) <- temp$cpg...466
    temp <- temp[, grep("GSM", names(temp), value=TRUE)]
    ages <- as.numeric(temp[1, ])
    betas <- as.matrix(temp[-1, ])
    
    
    for (start_age in seq(5, 70, by=5)) {
      end_age <- start_age + 20 
      print(paste0("split: ", j, " startage: ", start_age, " endage: ", end_age))
      age_window <- which(ages>=start_age & ages<=end_age)
      
      #check if window is within bounds (sample present)
      if (end_age > 100 || length(age_window) == 0) {
        next
      }
      
      betas_sub <- betas[, age_window, drop = FALSE]
      ages_sub <- ages[age_window]
      
      if (ncol(betas_sub) < 2 || length(ages_sub) < 2) {
        next  # Skip if there aren't enough data points to compute correlation
      }
      
      cor_results <- apply(betas_sub, 1, function(x) {
        if (sum(!is.na(x)) < 2) {
          return(NA)
        } else {
          return(cor(x, ages_sub, use="complete.obs", method="pearson"))
        }
      })
      
      beta_changes <- apply(betas_sub, 1, function(x) {
        max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      })
      
      sig_indices <- which(abs(cor_results) > 0.4 & beta_changes>=0.2)
      if (length(sig_indices) > 0) {
        print("Filling results")
        selected_cpg_sites <- rownames(betas)[sig_indices]
        key <- paste(start_age, "to", end_age)
        results_by_window[[key]] <- unique(c(results_by_window[[key]], selected_cpg_sites))
      }
    }
  }
  return(results_by_window)
}

results_by_window <- filterCorrelationsAndChange(num_groups, training=TRUE)

#print length of each list
for (i in 1:length(results_by_window)) {
  print(paste0("Start age: ", (i-1)*5+5, " End age: ", (i-1)*5+25, " Num CpG sites r>0.4 & Abs Change>=0.2: ", length(results_by_window[[i]])))
}

#print total unique sites
all_cpg_sites <- unlist(results_by_window)
unique_cpg_sites <- unique(all_cpg_sites)
print(length(unique_cpg_sites))

#down to 167387 sitesa
sorted_unique_cpg_sites <- sort(unique_cpg_sites)
write.csv(sorted_unique_cpg_sites, "7_15_cpgfilter2.csv")

#[1] "Start age: 5 End age: 25 Num CpG sites r>0.4 & Abs Change>=0.2: 77370"
#[1] "Start age: 10 End age: 30 Num CpG sites r>0.4 & Abs Change>=0.2: 8669"
#[1] "Start age: 15 End age: 35 Num CpG sites r>0.4 & Abs Change>=0.2: 75266"
#[1] "Start age: 20 End age: 40 Num CpG sites r>0.4 & Abs Change>=0.2: 87682"
#[1] "Start age: 25 End age: 45 Num CpG sites r>0.4 & Abs Change>=0.2: 6005"
#[1] "Start age: 30 End age: 50 Num CpG sites r>0.4 & Abs Change>=0.2: 2"
#[1] "Start age: 35 End age: 55 Num CpG sites r>0.4 & Abs Change>=0.2: 0"
#[1] "Start age: 40 End age: 60 Num CpG sites r>0.4 & Abs Change>=0.2: 0"
#[1] "Start age: 45 End age: 65 Num CpG sites r>0.4 & Abs Change>=0.2: 0"
#[1] "Start age: 50 End age: 70 Num CpG sites r>0.4 & Abs Change>=0.2: 0"
#[1] "Start age: 55 End age: 75 Num CpG sites r>0.4 & Abs Change>=0.2: 637"
#[1] "Start age: 60 End age: 80 Num CpG sites r>0.4 & Abs Change>=0.2: 1217"
#[1] "Start age: 65 End age: 85 Num CpG sites r>0.4 & Abs Change>=0.2: 950"
#[1] "Start age: 70 End age: 90 Num CpG sites r>0.4 & Abs Change>=0.2: 688"
#[1] 127642

```

Filter for sites that have variance \> 0.01 #46,233 sites left after
filtering

#using all datasets

\*\* might not need this one

```{r}
varianceFilter <- function(num_splits, folder_with_pre_filter, threshold=0.01, training=TRUE) {
  sites <- read.csv(folder_with_pre_filter)
  sites <- sites$x
  #subset_folder <- testing_subset_folder
  #if (training) {
  #  subset_folder <- training_subset_folder
  #}
  subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/combined_split3/"
  
  
  all_sites <- list()
  
  for (j in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_" , j, ".csv")) 
    temp <- as.data.frame(temp)
    #des_rows <- which(temp$cpg...207 %in% sites)
    #if (training) {
    #  des_rows <- which(temp$cpg...466 %in% sites)
    #}
    des_rows <- which(temp$cpg...466 %in% sites)
    rownames(temp) <- temp$cpg...466 #466 for training, 207 for testing
    temp <- temp[des_rows, ]
    temp <- temp[, grep("GSM", names(temp), value=TRUE)]
    
    variances <- apply(temp, 1, function(x) var(x, na.rm=TRUE))
    selected_variances <- variances[variances>threshold]
    selected_cpg_sites <- names(selected_variances)
    all_sites <- c(all_sites, selected_cpg_sites)
    
  }
  return(all_sites)
}

sites_filtered2 <- varianceFilter(num_groups, "7_15_cpgfilter2.csv", training=TRUE)
print(length(sites_filtered2))
write.csv(sites_filtered2, "7_15_cpgfilter3.csv")
```

Generate an i x j matrix, where i are CpG sites and j is ages 0-100
Params: training = TRUE is training, training = FALSE is testing average
= TRUE is generate an avg_beta matrix, average = FALSE is generate a
variance matrix (actually standard deviations)

```{r}
generateMatrix <- function(file_with_pre_filter, num_splits=num_groups, training=TRUE, average) {
  sites <- read.csv(file_with_pre_filter)
  sites <- sites[, -1]
  sites <- as.character(sites)
  #length: 288962
  NUMCPGSITE2 <- length(sites)
  #initialize output matrixs
  B <- matrix(NA, nrow=NUMCPGSITE2, ncol=101)
  age_range <- 0:100
  rownames(B) <- sites
  column_names <- character(101)
  for (i in age_range) {
    column_names[i+1] <- paste0("age_", i)
  }
  colnames(B) <- column_names
  subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/combined_split3/"
  
  #subset_folder <- testing_subset_folder
  #if (training) {
  #  subset_folder <- training_subset_folder
  #}
  
  for (i in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_", i, ".csv"))
    # Remove extra first column (automatic indices)
    temp <- temp[, -1]
    #filter for sites from correlation function
    #466 for training, 207 for testing
    #if (training) {
    #  des_rows <- which(temp$cpg...466 %in% sites)
    #}
    #if (!training) {
    #  des_rows <- which(temp$cpg...207 %in% sites)
    #}
    des_rows <- which(temp$cpg...466 %in% sites)
    temp <- temp[c(1, des_rows), ] #keep ages in the first row
    # Extract CpG site list
    identifiers <- temp[-1, ncol(temp)]
    for (j in 0:100) {
      # DEBUG // print(paste0("Split: ", i, " Age: ", j))
      des_cols <- which(temp[1, ]==j)
      # DEBUG // print(paste0("Columns at Age ", j, " : ", length(des_cols)))
      # Check if there are samples at age j
      if (length(des_cols)>0) {
        # Filter for age j
        age_data <- temp[, des_cols, drop = FALSE]
        # Drop age row
        age_data <- age_data[-1, ]
        
        #indice of identifiers corresponds to indice of age_data
        #age_data has samples as column names, no other identifiers
        # DEBUG // print(paste0("Dim of age_data: ", dim(age_data)))
        #this step takes the longest
        setDT(age_data)
        # Vectorized averaging
        if (average) {
          aggregated <- age_data[, .(final = rowMeans(.SD, na.rm = TRUE)), by = .I]
        }
        if (!average) {
          aggregated <- age_data[, .(final = apply(.SD, 1, sd, na.rm = TRUE)), by = .I]
        }        
        # DEBUG // print(paste0("Averaging done"))
        aggregated <- as.data.frame(aggregated)
        #if (training) {
        #  rownames(aggregated) <- as.character(identifiers$cpg...2441) #2441 for training, 1512 for testing
        #}
        #else {
        #  rownames(aggregated) <- as.character(identifiers$cpg...1512) #2441 for training, 1512 for testing
        #}
        rownames(aggregated) <- as.character(identifiers$cpg...3952)
        # Reassignment to matrix B
        B[rownames(aggregated), paste0("age_", j)] <- aggregated$final
        # DEBUG //print(paste0("Split: ", i, ", Age: ", j, " updated"))
      }
    }
  }
  return(B)
}

V <- generateMatrix("7_15_cpgfilter3.csv", num_splits=num_groups, training=TRUE, average=FALSE)
write.csv(V, "7_15_variance_matrix_v1.csv", row.names=TRUE)
B <- generateMatrix("7_15_cpgfilter3.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(B, "7_15_training_matrix_v1.csv", row.names=TRUE)
```

Compare average and variance matrices (for training data) and filter out
sites that have SD \> cutoff (calculated by transforming the beta
                              matrix) at more than 1 age point 26,247 sites left after filtering

```{r}
transformAVG <- function(x) {
  if (!is.na(x) && x >= 0 && x <= 1) {
    return(2 * sqrt(x * (1 - x)))
  } else {
    return(NA)  # Handle non-numeric or out-of-bound values appropriately
  }
}

filterBySDComparison <- function(path_to_variance_matrix, path_to_avg_matrix, transformAVG=transformAVG) {
  final_site_list <- list()
  #assume already formatted correctly here
  var <- readr::read_csv(path_to_variance_matrix)
  #var <- var[, colSums(is.na(var)) != nrow(var)]
  cpg_site_names <- var[[1]]
  var <- var[, -1, with = FALSE]
  var <- as.data.frame(var)
  rownames(var) <- cpg_site_names
  avg <- readr::read_csv(path_to_avg_matrix)
  #avg <- avg[, colSums(is.na(avg)) != nrow(avg)]
  avg <- avg[, -1, with = FALSE]
  avg <- as.data.frame(avg)
  rownames(avg) <- cpg_site_names
  #t <- apply(avg, c(1, 2), transformAVG)
  #rownames(t) <- rownames(avg)
  #colnames(t) <- colnames(avg)
  for (i in 1:nrow(avg)) {
    #DEBUGGING
    print(paste("ROW", i))
    count <- 0
    for (j in 1:ncol(avg)) {
      avg_value <- avg[i, j]
      adj <- 0.5 * sqrt(avg_value * (1 - avg_value)) # TWEAK VALUE FUNCTION
      variance <- var[i, j]
      if (!is.na(avg_value) && !is.na(adj) && !is.na(variance) && (variance > adj)) {
        count <- count + 1
      }
    }
    if (count<2) {
      final_site_list <- append(final_site_list, rownames(avg)[i])
      print(paste("Added Site", rownames(avg)[i]))
    } else {
      print("Too High Variance")
    }
  }
  return(final_site_list)
}

second_filter <- filterBySDComparison("6_30_variance_matrix_v1.csv", "6_30_training_matrix_v1.csv") 
print(length(second_filter))
write.csv(second_filter, "6_30_cpgfilter4.csv")
```

Regenerate training and testing matrices using generateMatrix() function
to move onto clustering

```{r}
Z <- generateMatrix("6_30_cpgfilter4.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(Z, "6_30_training_matrix_v2.csv", row.names=TRUE)

P <- generateMatrix("6_30_cpgfilter4.csv", num_splits=num_groups, training=FALSE, average=TRUE)
write.csv(P, "6_30_testing_matrix_v2.csv", row.names=TRUE)
```

Format and clean up matrices 9678 CpG sites used in further clustering

```{r}
B <- readr::read_csv("6_20_training_matrix_v2.csv")
G <- readr::read_csv("6_20_testing_matrix_v2.csv")
# Remove ages with no samples
# reduced to 89 cols
B <- B[, colSums(is.na(B)) != nrow(B)]
G <- G[, colSums(is.na(G)) != nrow(G)]

# Remove rows with NA values
#reduced to 227963 rows, could add imputation
#97958 in testing
B <- na.omit(B)
G <- na.omit(G)

# Reassign CpG site names to rownames and convert back to "data frame"
cpg_site_names <- B[[1]]
B <- B[, -1, with = FALSE]
B <- as.data.frame(B)
rownames(B) <- cpg_site_names
#df <- B

cpg_site_names <- G[[1]]
G <- G[, -1, with = FALSE]
G <- as.data.frame(G)
rownames(G) <- cpg_site_names

scaled_B <- t(scale(t(B)))
df <- as.data.frame(scaled_B)

scaled_G <- t(scale(t(G)))
df2 <- as.data.frame(scaled_G)

orderr <- rownames(df2)
#df2 is test
df <- df[orderr,]

df <- na.omit(df)
orderr <- rownames(df)

df2 <- df2[orderr, ]
```

Old clustering code below Using hierarchical clustering and cutting
dendrogram arbitrarily fitting B-splines with 5 df to visually represent
each cluster in ggplot

```{r}
d <- dist(df, method = "euclidean")  # Compute the distance matrix
hc <- hclust(d, method = "ward.D2")  # Perform hierarchical clustering

# Plot the dendrogram
plot(hc, main = "Dendrogram", xlab = "Sample index", ylab = "Height")

num_clusters = 8 #change number
clusters <- cutree(hc, k = num_clusters)
df$Cluster <- as.factor(clusters)

df <- df[order(df$Cluster), ]
df_order <- rownames(df)
df2 <- df2[df_order, ]

library(splines)
library(dplyr)
library(tidyverse)

for (i in 1:num_clusters) {
  train_df <- df[df$Cluster == i, ]
  ident <- rownames(train_df)
  test_df <- df2[ident, ]
  
  long_df_test <- test_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  long_df_train <- train_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample, Cluster), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  #fit <- loess(Beta ~ Age, data = long_df, span = 0.3)
  fit <- lm(Beta ~ bs(Age, df = 5), data = long_df_train)
  fit2 <- lm(Beta ~ bs(Age, df = 5), data = long_df_test)
  long_df_test$Fitted <- predict(fit2)
  long_df_train$Fitted <- predict(fit)
  
  
  
  ggplot(long_df_test, aes(x = Age, y = Beta)) +
    geom_point(size=1) +
    geom_line(aes(y = Fitted), color = "blue") +
    labs(title = paste("Test Cluster", i), x = "Age", y = "Scaled Beta Values") +
    theme_minimal()
  ggsave(
    paste0("test_cluster_", i, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300,
  )
  
  ggplot(long_df_train, aes(x = Age, y = Beta)) +
    geom_point(size=1) +
    geom_line(aes(y = Fitted), color = "blue") +
    labs(title = paste("Train Cluster", i), x = "Age", y = "Scaled Beta Values") +
    theme_minimal()
  ggsave(
    paste0("train_cluster_", i, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300,
  )
}

```

New clustering code here (in progress)

```{r}
df_matrix <- as.matrix(df)

# Set seed for reproducibility
set.seed(123)

# Initialize variables to store results
results <- list()
gaps <- numeric()

# Test for a range of K values
for (k in 2:10) {
  km.perm <- KMeansSparseCluster.permute(df_matrix, K = k, wbounds = seq(3, 7, len = 15), nperms = 5)
  results[[k]] <- km.perm
  gaps[k] <- max(km.perm$gaps)
}

optimal_k <- which.max(gaps)
cat("Optimal number of clusters:", optimal_k, "\n")

best_w <- results[[optimal_k]]$bestw
km.out <- KMeansSparseCluster(df_matrix, K = optimal_k, wbounds = best_w)
print(km.out)
plot(km.out)
```

Try time series clustering

```{r}
time_series_list <- tslist(df2)
n_clusters <- 7
clustering_result <- tsclust(time_series_list, type = "partitional", k = n_clusters, distance = "dtw_basic", seed = 123, trace = TRUE)
plot(clustering_result)
print(clustering_result@cluster)

cluster_centers <- clustering_result@centroids
print(cluster_centers)
summary(clustering_result)



```

Assign clusters back to df and order

```{r}
df$cluster <- clustering_result@cluster
df <- df[order(df$cluster),]
df_order <- rownames(df)
df2 <- df2[df_order, ]
```

Plot cluster results

```{r}
pheatmap(as.matrix(df[, -ncol(df)]),  # Exclude cluster column
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45,
         main="Training",
         breaks = seq(-5, 5, length.out = 101)
)
```

Order testing exactly like training

```{r}
df_order <- rownames(df)
df2 <- df2[df_order,]
```

Plot testing

```{r}
pheatmap(as.matrix(df2),  # Exclude cluster column
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45,
         main="Testing",
         breaks = seq(-5, 5, length.out = 101)
)
```

Generate Plots

```{r}
for (i in 1:num_clusters) {
  temp_df <- df[df$cluster == i, ]
  ident <- rownames(temp_df)
  temp_df2 <- df2[ident, ]
  
  long_df_test <- temp_df2 %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  long_df_train <- temp_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample, cluster), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  #fit <- loess(Beta ~ Age, data = long_df, span = 0.3)
  fit <- ls(Beta ~ bs(Age, df = 5), data = long_df_train)
  fit2 <- lm(Beta ~ bs(Age, df = 5), data = long_df_test)
  long_df_test$Fitted <- predict(fit2)
  long_df_train$Fitted <- predict(fit)
  
  q <- ggplot(long_df_train, aes(x = Age, y = Beta)) +
    geom_point(size=1) +
    geom_line(aes(y = fitted), color = "blue") +
    labs(title = paste("train cluster", i), x = "age", y = "scaled beta values") +
    theme_minimal()
  print(q)     
  
  
  p <- ggplot(long_df_test, aes(x = age, y = beta)) +
    geom_point(size = 1) +
    geom_line(aes(y = fitted), color = "blue") +
    labs(title = paste("test cluster", i), x = "age", y = "scaled beta values") +
    theme_minimal()
  print(p)
  
}
```


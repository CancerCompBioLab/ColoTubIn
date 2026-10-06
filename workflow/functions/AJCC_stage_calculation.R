ajcc_stage<- function(T, N, M) {
  # Convert to character to allow comparison
  T <- as.character(T)
  N <- as.character(N)
  M <- as.character(M)
  
  # Start as NA (Unknown)
  stage <- rep(NA_character_, length(T))
  
  # ---- Stage 4 ----
  stage[M == "1"] <- "4"
  
  # ---- Stage 3 ----
  stage[M == "0" & N %in% c("1", "2")] <- "3"
  
  # ---- Stage 2 ----
  stage[M == "0" & N == "0" & T %in% c("3", "4")] <- "2"
  
  # ---- Stage 1 ----
  stage[M == "0" & N == "0" & T %in% c("1", "2")] <- "1"
  
  # Return as ordered factor (1 < 2 < 3 < 4)
  factor(stage, levels = c("1","2","3","4"), ordered = TRUE)
}
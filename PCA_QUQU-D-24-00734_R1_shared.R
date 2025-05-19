## Shared script for: QUQU-D-24-00734 
## -- accompanying data and code
## -- replicates numerical contents and code listings in the manuscript
## -- For review process

## 19 May 2025 




## clear
rm(list = ls())
gc()


## PLEASE: make sure you set a working directory, first
## ----2) in RStudio just click Session--> Set Working Directory --> ...
## --- 1) ...or specify manually e.g. setwd("C:/Users/....")
##                                                 ####
## PREAMBLE                                        ####
##                                                 ####
## Functions declaration                           ####
## ---- 1 Own functions: PCA from scratch          ####


## ----  1.1 YAPCAB: own take on PCA                ####
YAPCAB <- function(mytab, myscale = "none", fix_eigen_sign = FALSE){
  
  ## Validate arguments thread: https://blog.r-hub.io/2022/03/10/input-checking/; https://stackoverflow.com/questions/19363827/validating-input-to-a-function-in-r
  in_arglist <- myscale %in% c("none", "log", "norm", "unit")
  stopifnot("admissible values for myscale: none, log, norm, unit" = in_arglist)
  
  ## centre or transform or normalise the data:               ####
  if(myscale == "norm"){
    ## -- center on the mean AND normalise to unit vector. WARNING: this isn't the default case always
    Y <- apply(mytab,2,function(x) (x - mean(x))/sqrt(sum((x - mean(x))^2)/(length(x)-1)))                  # std deviation (most used)
    # Y <- apply(mytab,2, scale)                                                                            # done in one go using "scale": see https://uc-r.github.io/pca
    message("data has been centred & normalised")
    ## CHECKS FOR NORMALISED DATA
    #round(apply(A_norm,2,mean))                                     # mean of each feature should be 0
    #round(apply(A_norm,2,var))                                      # variance of each feature should be 1
  } 
  if(myscale == "log"){
    log_data_tab <- apply(mytab,2,function(x){ifelse(x==0,0,log(x))})
    Y <- apply(log_data_tab,2,function(x) (x - mean(x)))
    message("data has been log-transformed and then centred")
    #round(apply(Y,2,mean),7)                                        # mean of each feature should be 0
  }
  if(myscale == "none"){
    Y <- apply(mytab,2,function(x) (x - mean(x)))
    #Y <- t(t(mytab) - apply(mytab,2,mean))                          # equivalent alternative
    message("data has only been centred")
    #round(apply(Y,2,mean),7)                                        # mean of each feature should be 0
  }
  if(myscale == "unit"){
    Y <- unitise_cols(mytab, ordinal = FALSE)
    Y <- apply(Y,2,function(x) (x - mean(x)))
    message("data has been unitised and then centred")
  }
  test_centred_data <- ifelse(as.numeric(sum(round(apply(Y,2,mean),7))) == 0, T, F)

  
  ## STAGE 1: PCA by eigenvalue decomposition                 ####
  ## -- caveat: pre-built scripts add a "symmetric" argument to the eigen decomposition
  ## -- however, even just mentioning "symmetric" regardless of FALSE or TRUE generates sign issues
  S <- 1/(nrow(Y)-1) * t(Y) %*% Y                              # covariance matrix of centered data
  eigenS <- eigen(S)                                           
  #eigenS_sym <- eigen(S, symmetric = T)                       # equivalent
  
  # CAVEAT: signs issue when using from scratch S and eigen
  if(fix_eigen_sign){
    ## debug only: show that when using cov() the first eigenvector is always positive
    # S1 <- cov(Y)
    # S2 <- cov(Y, use = "pairwise")                               # for comparison you must use THIS (as in package psych)
    # S_equal_test <- all.equal(S, S1)
    # eigen0 <- eigen(S)
    # eigen1 <- eigen(S1)
    # eigen2 <- eigen(S2)
    # isTRUE(all.equal(eigen1$vectors, eigen0$vectors))           ## will fail
    # isTRUE(all.equal(eigen2$vectors, eigen0$vectors))           ## will fail
    # eigen0$vectors[,1] <- eigen0$vectors[,1] * sign(eigen0$vectors[,1])               ## force the first eigenvector to be positive
    # isTRUE(all.equal(eigen1$vectors, eigen0$vectors))             ## will work

    ## enforce fixes to help compare with other implementations
    eigenS$vectors[,1] <- eigenS$vectors[,1] * sign(eigenS$vectors[,1])
    eigenS$values[eigenS$values < .Machine$double.eps] <- .Machine$double.eps
  }
  
  ## assign output
  lambda <- eigenS$values                                    ## eigenvalues (should be identical to variance of scores)    
  V_eig <- eigenS$vectors                                    ## "loadings" per component. Will get a differnt version later via SVD
  Z <- Y %*% V_eig                                           ## Scores
  
  ## check:  scores' variance = eigenvalues of cov matrix
  test_eigval_var <- all.equal(apply(Z,2,var), lambda)       ##  Eq.5 in paper
  
  ## name rows and cols of loadings and scores matrix
  colnames(Z) <- paste0("PC_",seq(1:ncol(Z))) 
  colnames(V_eig) <- paste0("PC_",seq(1:ncol(V_eig))) 
  rownames(V_eig) <- colnames(Y)
  names(lambda) <- colnames(V_eig)
  
  ## variance "explained" by each PC (reported along axis in PCA plot)
  var_retain_scores1 <- var(Z[,1])/sum(apply(Y,2,var))
  var_retain_scores2 <- var(Z[,2])/sum(apply(Y,2,var))
  
  ## Check equivalence w some pre-built options               ####
  scale_arg <- ifelse(myscale == "norm", T, F)
  if(myscale == "log") X_arg <- log_data_tab else X_arg <- mytab
  Y_prebuilt <- scale(X_arg, center = T, scale = scale_arg)[,]                                                                     
  S_prebuilt <- cov(Y)
  prebuilt_test_S <- all.equal(S_prebuilt, S)  
  prebuilt_test_Y <- all.equal(Y_prebuilt, Y)
  
  ## STAGE 2: SVD PCA                                         ####
  svd_Y <- svd(Y)                                                                               # singular value decomposition of a rectangular matrix of centred features
  U <- svd_Y$u                                                                                  # left eigenvectors
  V <- svd_Y$v                                                                                  # right eigenvectors equivalent to loading V_eig but different signs
  ell <- svd_Y$d                                                                                # singular values (NOT eigenvalues)
  D <- diag(as.vector(ell))                                                                     
  Z_svd <- U %*% D                                                                              # Paper eq. 10: scores matrix, SVD style
  
  ## label rows and cols
  colnames(D) <- rownames(D) <- names(lambda)
  colnames(U) <- colnames(V) <- names(lambda)
  rownames(U) <- rownames(Y)
  rownames(V) <- colnames(Y)
  colnames(Z_svd) <- colnames(Z)
  rownames(Z_svd) <- rownames(Z)
  
  ## check: singular values vs eigenvalues vs scores St.Dev.
  all.equal(abs(Z_svd), abs(Z))
  test_sv_ev <- all.equal(ell, sqrt((nrow(Y)-1)*eigen(S)$values))                               # Eq. 7 
  #ctest_sv_ev_v2 <- all.equal(ell, sqrt(nrow(Y)-1)*apply(unname(Z),2,sd))                      # equivalent: Eq. 7 assuming Eq. 5 holds
  
  
  ## STAGE 3: PC Biplot coord                                 ####
  A <- U                                                                                        # Eq.11 (i.e., Eq. 8 for alpha = 0): observations coordinates
  B_trsp <- D %*% t(V)                                                                          # Eq.12 (i.e., Eq. 8 for alpha = 0): features coordinates
  B <- t(B_trsp)
  
  ## check equations
  check_B_VD <- isEqualUnname(B, V%*%D)                                                          # checking that B = VD because later we'll need  v_jk*l_k
  check_A_from_Z  <- all.equal(A, sweep(Z_svd,2,ell, "/"))                                       # checking that A = ZD^{-1} = U from scaling scores (Eq.10).
  check_A_from_Zeig <-  all.equal(A, sweep(Z,2,svd_Y$d, "/"))                                    # won't work with scores from eigen decomposition unless in abs value
  
  ## SATGE 4: Geom properties tests                  ####
  ## ---- 4.1a feat vector length                    ####
  eq14a <- apply(B_trsp,2,function(x){norm(x, "2")/sqrt(nrow(Y)-1)})      # Eq.14 , second row
  eq14b <- sapply(1:ncol(Y),function(k){
    B_trsp[k,]^2                                                          # third line in eq.14  
  })
  eq14b <- sqrt(apply(eq14b,1,sum)/(nrow(Y)-1))
  eq14c <- sapply(1:ncol(Y),function(k){
    (V_eig[,k]*sqrt(lambda[k]))^2                                          # last line in eq.14
  })
  eq14c <- sqrt(apply(eq14c,1,sum))
  
  ## check equations
  check_eq13 <- all.equal(B %*% B_trsp, (nrow(Y)-1)*S)                    # Eq. 13
  check_eq14a <- all.equal(sqrt(diag(S)), eq14a)                          # Eq. 14
  check_eq14b <- all.equal(eq14a, eq14b)
  check_eq14c <- all.equal(eq14b, eq14c)
  
  if(check_eq13 & check_eq14a & check_eq14b & check_eq14c){
    check_eq14 <- T  
  } else {
    check_eq14 <- F
  }
  
  
  ## ---- 4.1b replicate Tab 7 in paper              ####
  length_b <- apply(B_trsp,2,function(x){norm(x, "2")})                   # length (2-norm) of the feature's coordinates vector
  length_b_r2 <- apply(B_trsp[1:2,],2,function(x){norm(x, "2")})          # length of rank-2 approximation of feature's coordinates vector
  check_length_b <-all.equal(length_b/sqrt(nrow(Y)-1), sqrt(diag(S)))     
  length_B <- cbind(length_b_r2, length_b, sqrt(diag(S)))
  rownames(length_B) <- colnames(Y)
  colnames(length_B) <- c("rank2", "allPCs", "sigma_y")
  
  
  
  
  
  ## ---- 4.3 Correlation vs cosine                  ####
  
  ## Prepare correlations grid
  idx_grid <- expand.grid(1:nrow(B_trsp), 1:nrow(B_trsp))                                       # rows or columns, will do (square matrix)
  idx_grid <- idx_grid[,c(2:1)]
  idx_grid <- idx_grid[which(idx_grid[,2] > idx_grid[,1]),]                                     # just keep the upper triangular part...
  
  ## cosine of angle between features vectors   
  ## part 1 - Rank 2 approximation (what is diplayed in a biplot)
  cos_theta_ab <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]
    ## only first two PCs (since B': rows = PCs, cols = features)
    load_a <- B_trsp[1:2,i]
    load_b <- B_trsp[1:2,j]
    (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
  })
  
  ## part 2:  replicating Eq. 15 (cosine of angle between features vectors)  requires full rank
  cos_theta_ab_allPCs <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]        #
    ## ALL pcs to compare to corr
    load_a <- B_trsp[,i]
    load_b <- B_trsp[,j]
    (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 
  })
  theta_ab <- acos(cos_theta_ab)                                                                #  angle from 0 to pi radians 
  theta_ab_deg <- theta_ab*180/pi                                                               # "degrees" instead of radians
  
  ## Pearsons product-moment correlation                                                        # Eq. 16 in paper but comapct version
  correl_feat <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]
    cor_ab <- cor(Y[,i], Y[,j])                                                                 # pre-built
    test_cor_ab <- cor.test(Y[,i], Y[,j])                                                       # get significance level of the correlation matrix:
    c(cor_ab, test_cor_ab$p.value)
    ## how to read the p-value https://www.jmp.com/en_gb/statistics-knowledge-portal/what-is-correlation/correlation-coefficient.html#404f1893-ae56-43ed-b84c-f6c99f313eca
    ## --- "The p-value is the probability of observing a non-zero correlation coefficient in our sample data when in fact the null hypothesis is true. 
    ## --- "the null hypothesis is typically that the observed relationship between the variables is the result of pure chance i.e., the correlation coefficient is really zero
    ## --- "A low p-value would lead you to reject the null hypothesis.
  })
  correl_feat <- as.data.frame(t(correl_feat))
  colnames(correl_feat) <- c("corr_coeff", "p-value")                                           #  this tests the H_0: true correlation IS equal to 0 [here is no linear relationship between the two variables]; alternative hypothesis: true correlation is not equal to 0
  
  ## Eq 16 from scracth
  correl_feat_fromScratch <- apply(idx_grid,1,function(x){                                      # this is actually Eq. 16 in paper
    i <- x[1]
    j <- x[2]       
    ## ALL pcs to compare to corr
    y_a <- Y[,i]
    y_b <- Y[,j]
    (y_a %*% y_b) / (norm(y_a, "2")*norm(y_b, "2"))                                             # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
  })
  test_corr_scratch <- isTRUE(all.equal(unname(correl_feat_fromScratch), correl_feat$corr_coeff))
  
  ## Check: cos angle vs correlation
  Eq18_corrcos_r2 <- all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab),6))          # won't work for rank 2, it's just an approx
  Eq18_corrcos <- all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab_allPCs),6))      # this should works
  
  ## table 4 export
  grid_labels <- t(apply(idx_grid, 1, function(x){
    colnames(B_trsp)[x]
  }))
  coscorr_export_tab <- cbind(grid_labels, 
                              round(theta_ab_deg,2), 
                              round(cos_theta_ab,2), 
                              round(cos_theta_ab_allPCs,2), 
                              round(correl_feat$corr_coeff,2), 
                              round(correl_feat$`p-value`,4))
  coscorr_export_tab <- as.data.frame(coscorr_export_tab)
  colnames(coscorr_export_tab ) <- c("a", "b",
                                     "theta_ab_deg",
                                     "cos_theta_ab_rank2",
                                     "cos_theta_ab",
                                     "corr_coeff",
                                     "p-value") 
  
  ## only 5% significant coeff
  coscorr_export_tab_signif <- coscorr_export_tab[which(coscorr_export_tab$`p-value` <= 0.05),]
  
  
  
  
  
  
  ## ---- 4.4 Correlation wiht PC                    ####
  ## ---- based on Johnson and Wichern Ch8 Result 8.3, 
  ## ---- 4.4.1 right hand side of Eq.19 
  eq19_rhs <- (1/sqrt(diag(S)))*t(t(V_eig)*sqrt(lambda))                                # Eq.19, right hand side, compact operation for all j and k
  
  ## check remark 6: Eq. 14 contains the necessary elements
  alt_eq14_sqrtArg2 <- sapply(1:ncol(Z),function(k){
    (V_eig[,k]*sqrt(lambda[k]))                                                         # last row in eq.14 'as is'
  })
  eq19_rhs2 <- sweep(alt_eq14_sqrtArg2,1,sqrt(diag(S)),"/")                             # Eq.19 in a way that confirms remark 6
  remark_6 <- all.equal(unname(eq19_rhs), unname(eq19_rhs2))
  
  ## ---- 4.4.2 left hand side of Eq.19
  ##  corr coeff
  idx_grid3 <- expand.grid(1:ncol(Y),1:2)                                               # a) Prepare correlations grid
  idx_grid3 <- idx_grid3[,c(2:1)]
  correl_feat_PC <- apply(idx_grid3,1,function(x){                                      # this is actually Eq. 16 in paper
    k <- x[1]
    j <- x[2]       
    z_k <- Z[,k]
    y_j <- Y[,j]
    (z_k %*% y_j) / (norm(z_k, "2")*norm(y_j, "2"))                                     # Pearsons product-moment correlation from scratch
  })
  check_eq19 <- all.equal(correl_feat_PC, c(eq19_rhs[,1:2]))                            # verify equivalence
  
  ## redo with pre-built corr (returns p-value)
  correl_feat_PC_prebuilt <- apply(idx_grid3,1,function(x){                             # this is actually Eq. 16 in paper
    k <- x[1]
    j <- x[2] 
    z_k <- Z[,k]
    y_j <- Y[,j]
    cor_ab <- cor(z_k , y_j)
    test_cor_ab <- cor.test(z_k , y_j)
    c(cor_ab, test_cor_ab$p.value)
  })
  feat_PC_cor_tab <- as.data.frame(cbind(idx_grid3, t(correl_feat_PC_prebuilt)))
  colnames(feat_PC_cor_tab ) <- c("feat", "PC", "corr", "p-value")
  all.equal(correl_feat_PC,  feat_PC_cor_tab$corr ) 
  
  
  
  ## ---- 4.5 obs vec: Mahalanobis distance          ####
  ## Jolliffe P.77-78; Gabriel 1971; and also mentioned by Venables and Ripley
  
  ## small version for arbitrary observations pair: what I use in the paper due to space
  # x <- as.matrix(Y[1,] - Y[2,])
  # x2 <- as.matrix(A[1,] - A[2,])
  # mymahala <- as.numeric(t(x) %*% solve(S, x))
  # adj_dist2 <-as.numeric((nrow(Y)-1)*t(x2)%*%x2) 
  # all.equal(adj_dist2, mymahala)
  
  ## all observation pairs
  idx_grid2 <- expand.grid(1:nrow(Y), 1:nrow(Y))
  idx_grid2 <- idx_grid2[,c(2:1)]
  idx_grid2 <- idx_grid2[which(idx_grid2[,2] > idx_grid2[,1]),]                                     # just keep the upper triangular part...
  
  mymahala_all <- apply(idx_grid2, 1, function(k){
    i <- k[1]
    j <- k[2]
    x <- as.matrix(Y[i,] - Y[j,])
    as.numeric(t(x) %*% solve(S, x))
  })
  
  ## now the squared euclidean distance between biplot coordinates
  sqdist_biplot_coord <- apply(idx_grid2, 1, function(k){
    i <- k[1]
    j <- k[2]
    x <- as.matrix(A[i,] - A[j,])
    t(x)%*%x
  })
  
  ## check the second is "proportional" to the first
  check_eq20 <- all.equal((nrow(Y)-1)*sqdist_biplot_coord, mymahala_all)
  
  ## dist table to export
  observ_pairs <- matrix(rownames(Y)[as.matrix(idx_grid2)], nco=2, byrow = F)
  mahala_tab <- as.data.frame(cbind( observ_pairs, round(mymahala_all,2), round(sqdist_biplot_coord,2)))
  colnames(mahala_tab) <- c("observ_h", "observ_i", "mahalanobis_dist2_data", "euclidean_dist2_biplot")
  
  
  
  ## prepare output                                  ####
  list(data_modified = Y, 
       covar_mat = S,
       PCA_eigen = list(eigenvalues = lambda, 
                        loadings = V_eig, 
                        scores = Z),
       var_explained = list(variance_expl_PC1 = var_retain_scores1,
                            variance_expl_PC2 = var_retain_scores2),
       PCA_svd = list(U = U,
                      V = V,
                      D = D,
                      Z_svd = Z_svd),
       biplot = list(A = A,
                     B_trsp = B_trsp,
                     b_len_sd = length_B,
                     coscorr = coscorr_export_tab,
                     feat_PC_cor_tab = feat_PC_cor_tab,
                     dist_mahala = mahala_tab),
       checks = list(main_checks = list(Eq5_eigval_is_var = test_eigval_var,
                                        Eq7_sv_ev = test_sv_ev,
                                        Eq11_A_Z = check_A_from_Z,
                                        Eq14_1_b_len = check_eq14,
                                        Eq18_corrcos = Eq18_corrcos,
                                        Eq19_corrPC = check_eq19,
                                        Eq20_mahala = check_eq20 ), 
                     peace_of_mind_checks = list(
                       data_is_centred = test_centred_data,
                       same_as_prebuilt_S = prebuilt_test_S,
                       same_as_prebuilt_Y = prebuilt_test_Y)
       ),
       flag_transform = myscale)
  
}

## ----  1.2 Reconstruct data form PCs              ####
CoordApprox <- function(YAPCAB_out, n_PC = 2){
  #### Reconstruction (prediction) from PC scores ####
  ## Think of image compression: we seek to reconstruct the p-dimensional [p = n of features] coordinates of an observation [i.e. a row of the original data matrix]
  ## from a lower-dimensional projection of the original data onto "fewer than p" principal components [e.g. only PC1 and PC2, which are usually meant to be a "good enough" approx]
  ## The idea of "projection matrix"
  # -- thread: https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
  # -- According to  Govwer and Hand (1996) Biplots. p.12 the "scores" matrix Z [centred data matrix * eigenvectors of the covariance matrix] is the "best display of the n points"
  # ----- and the projection of a point x in a 2-dimensional space is given by z = x * [matrix consisting of the first 2 principal eigenvectors]
  # ----- and the prediction of a point z in the original, say 3-dimensional space is given by x = z * TRANSPOSE[matrix consisting of the first 2 principal eigenvectors]
  
  ## how many components to use in the reconstruction?
  if( (n_PC < 2) | (n_PC >= ncol(PCA_out$data_modified)) ){
    ## if the number specified is stupid, revert back to 2
    n_PC = 2
  }
  use_n_PCs <- n_PC
  
  ## for illustration: all observations
  ## -- reconstruct p-coordinates [p is the number of features] for all observations using 2 (i.e, just PC1 and PC2) or more of their projections
  ## -- this will reconstruct ALL higher-dimentional points from their "lower dimensional" approximations (first two or three, four etc. eigenvectors)
  coord_reconstruct_2PC <- PCA_out$scores[,1:n_PC] %*% t(PCA_out$loadings[,1:n_PC])
  
  ## -- NOW for a benchmark use ALL PCs i.e. recreates the original higher dimensional coordinates - "perfect reconstruction"
  coord_perfect_reconstr <- PCA_out$scores%*% t(PCA_out$loadings)
  test_perfect_reconstr <- round(coord_perfect_reconstr - PCA_out$data_modified,6)
  
  ## how much resolution have we lost by using fewer PCs than the number of features?
  ## -- please notice I've made this up: to be refined in case there's something better in the literature
  length_a <- norm(c(as.matrix(coord_perfect_reconstr)), "2")
  length_b <- norm(c(coord_reconstruct_2PC), "2")
  reconstruciotn_loss <- abs(length_a - length_b) / length_b
  
  ## OUTPUT
  list(reconstructed_coordinates = coord_reconstruct_2PC, 
       PCs_used = n_PC, 
       reconstruciotn_loss = reconstruciotn_loss)
}

## ----  1.3 Draw PC biplot, base R                 ####
showmetheBplot <- function(YAPCAB_out, hires = FALSE){
  
  if(hires){
    ## The next bit is to export a higher figure resolution, goes before the actual plot 
    ## thread: https://stackoverflow.com/questions/22813723/how-can-i-increase-the-resolution-of-my-plot-in-r
    ## https://youtu.be/1SyLtXskq2g?si=LxLg6_B8ZZlS7JhJ
    png(paste0(Out_path,"/biplot_extd_higherres.png",collapse = " "), res=300, 
        width=800, 
        height=800)
  }
  
  ## actual biplot                                      
  biplot_feat_coord <- t(YAPCAB_out$biplot$B_trsp[1:2,])
  biplot_observ_coord <- YAPCAB_out$biplot$A[,1:2]
  Y <- YAPCAB_out$data_modified
  lambda <- YAPCAB_out$PCA_eigen$eigenvalues
  
  ## composite figure - 3-part layout
  layout(mat = matrix(c(2, 1, 0, 3), 
                      nrow = 2, 
                      ncol = 2),
         heights = c(0.5, 2),    # Heights of the two rows
         widths = c(2, 0.5))     # Widths of the two columns
  
  ## Plot 1: Scatterplot
  ## -- fix margins so to show the axis' labels
  par(mar = c(4, 4, 0, 0))
  ## -- Scale down features' coordinates as they may be way larger than the observations'
  adjust_viz <- abs(min(biplot_observ_coord)/min(biplot_feat_coord))
  biplot_feat_coord_viz <- biplot_feat_coord*adjust_viz
  ## -- Shorten features' labels for display
  ## -- thread https://stackoverflow.com/questions/11776287/remove-pattern-from-string-with-gsub
  attr_short_label <- gsub(".*_","", colnames(Y))
  
  
  ## -- Set axis' range considering both features and observations
  lim_max1ob <- max(biplot_observ_coord[,1])*1.3
  lim_min1ob <- min(biplot_observ_coord[,1])*1.3
  lim_max2ob <- max(biplot_observ_coord[,2])*1.3
  lim_min2ob <- min(biplot_observ_coord[,2])*1.3
  lim_max1feat <- max(biplot_feat_coord_viz[,1])*1.3
  lim_min1feat <- min(biplot_feat_coord_viz[,1])*1.3
  lim_max2feat <- max(biplot_feat_coord_viz[,2])*1.3
  lim_min2feat <- min(biplot_feat_coord_viz[,2])*1.3
  
  ## keep ratio 1:1 as recommended by Gower
  axis_range_1 <- c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat))
  axis_range_2 <- c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat))
  
  axis_ranges <- rbind(axis_range_1, axis_range_2)
  axis_range_idx <- apply(abs(axis_ranges),2,which.max)
  common_axis_range <- axis_ranges[cbind(axis_range_idx,1:ncol(axis_ranges))]
  
  ## -- plot the actual biplot
  par(mar = c(4, 4.2, 0, 0))     
  ## -- observations' plot
  plot(biplot_observ_coord[,1:2],
       pch = 18,                # point shape
       cex = 2,               # point size
       #xlim = c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat)),
       #ylim = c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat)),
       xlim = common_axis_range,
       ylim = common_axis_range,
       xlab = paste0("PC 1 (var expl: ", round(lambda[1]/sum(lambda),2)*100,"%)"),
       ylab = paste0("PC 2 (var expl: ", round(lambda[2]/sum(lambda),2)*100,"%)")
  )
  abline(h = 0, v = 0, col = "gray60")
  text(biplot_observ_coord[,1], biplot_observ_coord[,2], cex=1.5, pos = 3, labels = rownames(Y))
  ## Overlay "loading plot"
  arrows(rep(0,2),rep(0,2), biplot_feat_coord_viz[,1], biplot_feat_coord_viz[,2], length = .1, angle =10, col = "purple")
  text(biplot_feat_coord_viz[,1], biplot_feat_coord_viz[,2], pos = 4, cex=1.5, labels = attr_short_label)
  
  ## instead of a dual axis I will plot "around" the main plot
  ## -- features' axes range
  # lim_max1 <- max(biplot_feat_coord[,1])
  # lim_min1 <- min(biplot_feat_coord[,1])
  # lim_max2 <- max(biplot_feat_coord[,2])
  # lim_min2 <- min(biplot_feat_coord[,2])
  # x_range_feat <- c(lim_min1, lim_max1)
  # y_range_feat <- c(lim_min2, lim_max2)
  
  x_range_feat <- c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat))/adjust_viz
  y_range_feat <- c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat))/adjust_viz
  
  axis_ranges2 <- rbind(x_range_feat, y_range_feat)
  axis_range_idx2 <- apply(abs(axis_ranges2),2,which.max)
  common_axis_range2 <- axis_ranges2[cbind(axis_range_idx2,1:ncol(axis_ranges2))]
  
  
  ## Plot 2: Top - this will represent the horizontal scale for the feature's coordinates 
  par(mar = c(0, 4.2, 2.2, 0))
  x <- data.frame(biplot_feat_coord[,1],1)                      # notice multiplication by -1 to help overlaying the two chard
  plot(x, 
       # xlim = c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat))/adjust_viz,
       xlim = common_axis_range2, 
       axes = F,
       ylab = "",
       xaxt="n",                                   #remove tick marks
       type = "o",
       col = "purple",
       #pch = 15                                    # shape
  )
  text(x, pos=2, labels=attr_short_label, cex= 1.5)
  axis(side=3, cex.axis = 1, col = "grey", las=1)
  
  ## Plot 3: Right - feature's coordinates
  par(mar = c(4, 0, 0, 2.2))
  x2 <- data.frame(1,biplot_feat_coord[,2])                      # notice multiplication by -1 to help overlaying the two chard
  plot(x2, 
       #ylim = c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat))/adjust_viz,
       ylim =common_axis_range2,
       axes = F,
       xlab = "",
       xaxt="n",                                   #remove tick marks
       type = "o",
       col = "purple",
       #pch = 15                                    # shape
  )
  text(x2, pos=3, labels=attr_short_label, cex= 1.5) 
  axis(side=4, cex.axis = 1, col = "grey", las=1)
  
  if(hires){
    ## close device (when saving higher res)
    dev.off()
  }
}

## ----  1.4 Unitise data (or ordinal to numeral)   ####
unitise_cols <- function(mytab, ordinal = FALSE){
  ## Original version converted ordinal to numeric (assumes input is ranking of ordinal variable)
  ## -- ref: Kaufman & Rousseeuw 2005 Finding Groups in Data  p29-31
  ## Adapted here for Biplot: based on la Grange et al (biplotGUI JSS article Appendix) - called "unitising" transformation
  ## -- all same except min
  z_f <- apply(mytab,2,function(r){
    M <- max(r)
    if(ordinal){
      m1 <- 1
      m2 <- 0
    } else {
      m1 <- m2 <- min(r)
    }
    (r - m1) / (M - m2)  
  })
  return(z_f)
}



## ---- 2 Own functions for testing properties     ####

## test 1 length of feat vector and feat std dev
## Check eq. 14: length of features vectors in biplot vs feats St.Dev
lenTest <- function(myB_trsp, Y, S = NULL){
  sd_y <- apply(Y,2,sd)                                     # to avoid passing the covariance mat S as an input. Otherwise same as sqrt(diag(S)) as per equation
  if(!is.null(S)){
    testS <- isTRUE(all.equal(unname(sd_y),
                       unname(sqrt(diag(S))))
    )
    sd_y <- sqrt(diag(S))                                   # now use the diagonal elements of S, not the sd of the features
    if(!testS){
      warning("Test failed: the diagonal elements of S are not equal to the sd of matrix Y's columns")
    } 
  } else  {
    testS <- NA
    warning("No covariance matrix provided. Using St.Dev. of features instead")
  }
  length_B <- apply(myB_trsp,2,function(x){norm(x, "2")})                           # Eq.14 , second row
  eq14a <- length_B/sqrt(nrow(Y)-1)
  test_len <- isTRUE(all.equal(unname(sd_y), unname(eq14a)))                      # Eq. 14
  
  if(!is.null(S)){
    check_eq13 <- isTRUE(all.equal(t(myB_trsp) %*% myB_trsp, (nrow(Y)-1)*S))          # Eq. 13
  } else {
    check_eq13 <- F
  }
  
  ## table output
  if(!is.null(S)){
    sqdiagS <- sqrt(diag(S))
  } else {
    sqdiagS <-rep(0,length(eq14a))
  }
  tab_out <- cbind(length_B, 
                   eq14a, 
                   apply(Y,2,sd), 
                   sqdiagS
                   )
  colnames(tab_out) <- c("length B", "RHS eq14", "LHS eq14", "LHS eq4 S")
  
  ## output
  list(length_B = length_B,
       eq14 = eq14a,
       tab_out = tab_out,                 # puts RHS and LHS of Eq.14 next to each other
       testS = testS,                     # this could fail if the covariance matrix is not computed correctly
       test_len = test_len,               # this could still work if we consider the sd of the columns of the data matrix
       check_eq13 = check_eq13)
}


## test 2 Correlation vs cosine of features                
corcosTest <- function(myB_trsp, myY){
  ## Prepare correlations grid
  idx_grid <- expand.grid(1:nrow(myB_trsp), 1:nrow(myB_trsp))
  idx_grid <- idx_grid[,c(2:1)]
  idx_grid <- idx_grid[which(idx_grid[,2] > idx_grid[,1]),]                                     # just keep the upper triangular part...
  
  ## cosine of angle between features vectors                                                   # Eq. 15 in paper
  ## -- Rank 2 approximation (diplayed in a biplot)
  cos_theta_ab <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]
    ## only first two PCs
    load_a <- myB_trsp[1:2,i]
    load_b <- myB_trsp[1:2,j]
    (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 
  })
  
  ## -- no dimension reduction: the equivalence between cos theta and corr requires this
  cos_theta_ab_allPCs <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]        #
    ## ALL pcs to compare to corr
    load_a <- myB_trsp[,i]
    load_b <- myB_trsp[,j]
    (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 
  })
  
  theta_ab <- acos(cos_theta_ab)                                                                #  angle from 0 to pi radians 
  theta_ab_deg <- theta_ab*180/pi                                                               # "degrees" instead of radians
  
  ## Pearsons product-moment correlation                                                        # Eq. 16 in paper but comapct version
  correl_feat <- apply(idx_grid,1,function(x){
    i <- x[1]
    j <- x[2]
    cor_ab <- cor(myY[,i], myY[,j])
    ## we might be interested in the significance level of the correlation matrix:
    ## --- https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/#correlation-test
    ## --- for nice viz see thread https://rpubs.com/MajstorMaestro/240657
    test_cor_ab <- cor.test(myY[,i], myY[,j])
    c(cor_ab, test_cor_ab$p.value)
    ## how to read the p-value https://www.jmp.com/en_gb/statistics-knowledge-portal/what-is-correlation/correlation-coefficient.html#404f1893-ae56-43ed-b84c-f6c99f313eca
    ## --- "The p-value is the probability of observing a non-zero correlation coefficient in our sample data when in fact the null hypothesis is true. 
    ## --- "the null hypothesis is typically that the observed relationship between the variables is the result of pure chance i.e., the correlation coefficient is really zero
    ## --- "A low p-value would lead you to reject the null hypothesis.
  })
  correl_feat <- as.data.frame(t(correl_feat))
  colnames(correl_feat) <- c("corr_coeff", "p-value")                                           #  this tests the H_0: true correlation IS equal to 0 [here is no linear relationship between the two variables]; alternative hypothesis: true correlation is not equal to 0
  
  ## Alternative, equivalent 
  test_correl_feat_fromScratch <- apply(idx_grid,1,function(x){                                      # this is actually Eq. 16 in paper
    i <- x[1]
    j <- x[2]       
    ## ALL pcs to compare to corr
    y_a <- myY[,i]
    y_b <- myY[,j]
    (y_a %*% y_b) / (norm(y_a, "2")*norm(y_b, "2"))                                             # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
  })
  test_corr_scratch <- isTRUE(all.equal(unname(test_correl_feat_fromScratch), correl_feat$corr_coeff))
  
  ## Check n.4: cos angle vs correlation
  coscorr_test_r2 <- isTRUE(all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab),6)))        # won't work for rank 2, it's just an approx
  coscorr_test <- isTRUE(all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab_allPCs),6)))    # works
  
  
  ## table to export
  #coscorr_export_tab <- cbind(idx_grid, theta_ab_deg, cos_theta_ab, cos_theta_ab_allPCs, correl_feat)
  
  grid_labels <- t(apply(idx_grid, 1, function(x){
    colnames(myB_trsp)[x]
  }))
  coscorr_export_tab <- cbind(grid_labels, 
                              round(theta_ab_deg,2), 
                              round(cos_theta_ab,2), 
                              round(cos_theta_ab_allPCs,2), 
                              round(correl_feat$corr_coeff,2), 
                              round(correl_feat$`p-value`,4))
  coscorr_export_tab <- as.data.frame(coscorr_export_tab)
  colnames(coscorr_export_tab ) <- c("a", "b",
                                     "theta_ab_deg",
                                     "cos_theta_ab_rank2",
                                     "cos_theta_ab",
                                     "corr_coeff",
                                     "p-value") 

  ## output
  list(coscorr_export_tab = coscorr_export_tab,
       coscorr_test =  coscorr_test,
       coscorr_test_r2 = coscorr_test_r2,
       test_corr_scratch)
  
}

## test 3 Correlation between feature and PC
corPCTest <- function(V, Z, Y, lambda, S = NULL){
  ## even if S is not provided, we can work things out from Y
  sd_y <- apply(Y,2,sd)                                     # to avoid passing the covariance mat S as an input. Otherwise same as sqrt(diag(S)) as per equation
  if(!is.null(S)){
    testS <- isTRUE(all.equal(unname(sd_y),
                              unname(sqrt(diag(S))))
    )
    sd_y <- sqrt(diag(S))                                   # now use the diagonal elements of S, not the sd of the features
    if(!testS){
      warning("Test failed: the diagonal elements of S are not equal to the sd of matrix Y's columns")
    } 
  } else  {
    testS <- NA
    warning("No covariance matrix provided. Using St.Dev. of features instead")
  }
  ## based on Johnson and Wichern Ch8 Result 8.3, 
  eq19_rhs <- (1/sd_y)*t(t(V)*sqrt(lambda))                                # Eq.19, right hand side, compact operation for all j and k
  eq19_rhs_alt <- (1/apply(Y,2,sd))*t(t(V)*sqrt(lambda))
  
  idx_grid3 <- expand.grid(1:ncol(Y),1:ncol(Z))                            # Prepare correlations grid
  idx_grid3 <- idx_grid3[,c(2:1)]
  correl_feat_PC <- apply(idx_grid3,1,function(x){                         # this is actually Eq. 16 in paper
    k <- x[1]
    j <- x[2]       
    z_k <- Z[,k]
    y_j <- Y[,j]
    (z_k %*% y_j) / (norm(z_k, "2")*norm(y_j, "2"))                              # Pearsons product-moment correlation from scratch
  })
  check_eq19 <- isTRUE(all.equal(correl_feat_PC, c(eq19_rhs[,1:ncol(Z)])))       # verify equivalence
  
  feat_PC_cor_tab <- as.data.frame(cbind(idx_grid3, correl_feat_PC, c(eq19_rhs), c(eq19_rhs_alt)))
  colnames(feat_PC_cor_tab ) <- c("PC", "feat", "corr", "eq19RHS_S", "eq19RHS_sd_y")                          
  all.equal(correl_feat_PC,  feat_PC_cor_tab$corr ) 
  
  ## outout
  list(corr_feat_PC =  feat_PC_cor_tab,
       testS = testS,                            # this will fail if the covariance matrix is incorrectly specified
       check_eq19 = check_eq19)                  # this might still work if we use directly the St.Dev of matrix Y's columns
}


## test 4 Euclidean distance between observ in biplot = Mahalanobis distance between observations in data
MahalaTest <- function(A, Y, S){
  ## MUST PROVIDE S (sample covariance matrix of Y)
  
  ## Jolliffe P.77-78; Gabriel 1971; and also mentioned by Venables and Ripley
  ## -- one main issue is whether the relaionship should be "exact2" or "proportional
  ## -- Jolliffe explains the trick: Gabriel and bpca multiply and divide by sqrt(n-1)
  ## -- the same trick has impact on the relationship between feature vector's length and feature St.dev
  if(!is.matrix(A)){A <- as.matrix(A)}
  if(!is.matrix(Y)){Y <- as.matrix(A)}
  ## all rows pair from the centred data matrix
  idx_grid2 <- expand.grid(1:nrow(Y), 1:nrow(Y))
  idx_grid2 <- idx_grid2[,c(2:1)]
  idx_grid2 <- idx_grid2[which(idx_grid2[,2] > idx_grid2[,1]),]                                     # just keep the upper triangular part...
  
  mymahala_all <- apply(idx_grid2, 1, function(k){
    ## recall that there is a built-in function:   mahalanobis(Y[i,],center = Y[j,], cov = S) 
    i <- as.numeric(k[1])
    j <- as.numeric(k[2])
    x <- as.matrix(Y[i,] - Y[j,])                                    
    t(x) %*% solve(S, x)
  })
  
  ## now the squared euclidean distance between biplot coordinates
  sqdist_biplot_coord <- apply(idx_grid2, 1, function(k){
    i <- as.numeric(k[1])
    j <- as.numeric(k[2])
    x <- as.matrix(A[i,] - A[j,])
    t(x)%*%x
  })
  ## check the second is "proportional" to the first
  mTest <- all.equal((nrow(Y)-1)*sqdist_biplot_coord, mymahala_all)
  
  ## table of square distances
  dist_tab <- cbind(idx_grid2, mymahala_all, (nrow(Y)-1)*sqdist_biplot_coord, sqdist_biplot_coord)
  colnames(dist_tab) <- c("from", "to", "mahala dist sq (data)", "euclid dist sq (biplot) x const", "euclid dist sq (biplot)")
  
  ## output
  list(mahalaTest = mTest,
       dist_tab = dist_tab)
}


## ---- 3 Replica of selected pre-built functions  ####

replica_biplot <- function (x, myscale = 1, original = TRUE,  pc.biplot = FALSE, my_princomp = FALSE){
  ## my version of base R biplot with some fixes
  ## stats::biplot does not return an output for the users hence the need for this replica
  ## The following is based on direct scripts inspection
  ## --- View(stats:::biplot.default)
  ## --- View(stats:::biplot.prcomp)
  ## --- View(stats:::biplot.princomp)
  
  ## allow scale to be only 0 or 1 (the original function does not take the complement)
  if(is.na(match(myscale, c(0,1)))){
    stop("the argument scale can only be 0 or 1 sorry")
  }
  
  ## I have added this bit so that it handles multiple objects
  if(inherits(x, "princomp") | my_princomp){
    ## PCA done using eigenvalue/eigenvector decomposition
    Z <- x$scores 
    V <- x$loadings
  }
  if(inherits(x, "prcomp")){
    ## PCA done using SVD
    Z <- x$x                                       # Called scores in  stats:::biplot.prcomp                                 
    V <- x$rotation                                
  }
  n <- NROW(Z)
  
  ## basic arguments
  choices = 1:ncol(Z)                              # default 1L:2L. Some of my tests require all dimension so I am imposing this
  scale = myscale                                  # default = 1. Corresponds to alpha in my paper's notation but scale = 1 means alpha = 0 (PC biplot) reminds me of Venables and Ripley 2002 Ch.1
  pc.biplot = FALSE                                # a bit misleading. It aims to implements Eq. 52 in Gabriel 1971 (from ?stats::biplot.prcomp). (see Jolliffe Ch.5 who explains what this workaround is for in terms of geometrical properties)
  
  if(original){
    ## will genrate discrepancies
    lam <- x$sdev[choices]                         # if prcomp, then sdev = sv/sqrt(n-1) = scores' st.dev., consistently with Eq.7 in my paper 
    lam <- lam * sqrt(n)                           # gets sv. however, it should be: sv = sdev * qrt(n-1). This is a discrepanacy
    if (scale != 0)                                # this is not quite correct unless we are considering just 1 vs 0 (hence why I restricted the input)
      lam <- lam^scale                             # D^1 = D. Leads to A = ZD^{-1} = U and B^T = DV^T i.e., a PC biplot (yet it is not enforced by pc.biplot as the documentation says)
    else lam <- 1                                  # D^0 = Identity. in this case A = scores and B = loadings
    if (pc.biplot) lam <- lam/sqrt(n)              # From the documentation, this implements Eq. 52 in Gabriel 1971 - but does not enforce scale = 1 (as the documentation says)
    A <- t(t(Z[,choices])/lam)                     # if pc.plot = FALSE and scale = 1 this is Eq.11 in paper; A = ZD^{-1} = U (yet D is wrong because of the multiplication by sqrt(n) not sqrt(n-1)); otherwise we get A = Z for scale = 0
    B <- t(t(V[,choices])*lam)                     # if pc.plot = FALSE and scale = 1 this is Eq.12 in paper; B^T = DV^T (yet D is wrong because of the multiplication by sqrt(n) not sqrt(n-1)); otherwise we get B = V
  } else {
    ## Apply "fixes" to address discrepancy
    lam <- x$sdev[choices]*sqrt(n-1)               # fixed sqrt(n-1) see eq.7 in my paper
    if (pc.biplot){
      scale <- 1                                   # added: needs to be enforced
      lam <- lam/sqrt(n-1)                         # Eq. 52 in Gabriel 1971 but consistent with sample covariance matrix of Y. 
    } else {
      if (scale == 1)                                # assume it's either 1 or zero: 1 happens to be correct for a PC biplot 
        lam <- lam^scale                             # D^1 = D. In this case A = A = ZD^{-1} = U and B^T = DV^T 
      else lam <- 1                                  # D^0 = Identity. in this case A = scores and B = loadings (alpha = 1)
    }
    A <- t(t(Z[,choices])/lam)                     # observation coordinates as per Eq.11 in my paper; A = ZD^{-1} = U which is the case for alpha = 0 (i.e. pc biplot)
    B <- t(t(V[,choices])*lam)                     # features coordinates as per Eq.12 in my paper; B^T = DV^T
   
    # ## equivalent ways of computing Eq.11
    # all.equal(t(t(Z)/lam), 
    #           Z%*%solve(diag(lam)))
    # all.equal(t(t(Z)/lam), sweep(Z,2,x$sdev,"/"))
  }
  ## out
  list(obser_coord = A, feat_coord = B)
}

replica_princomp <- function(X, original = T, test_symmetric = T, fix_sign = T){
  ## my version of base R princomp with some fixes
  Y <- scale(X, scale = F)[,]
  if(original){
    ## make the same mistakes if needed but allow to play with some arguments for testing
    ## cov matrix step 1 (raw data first)
    pca10_internal_covmat <- cov.wt(X)                            # not a mistake yet although this is an unusual step: covariance matrix  of the raw data FIRST...
    test_cov_ok1 <- all.equal(pca10_internal_covmat$cov, cov(X))                    # NOTICE: this means at this stage it divides by n - 1, not n
    ## cov matrix step 2 (center the data)
    pca10_internal_n.obs <- pca10_internal_covmat$n.obs
    pca10_internal_cv <- pca10_internal_covmat$cov * (1 - 1/pca10_internal_n.obs)   # This is where the division by n occurs
    ## covariance matrix oddity: same same yet different...
    S_wrong <- (1/nrow(Y))*t(Y)%*%Y                                  # from function help: "The calculation is done using eigen on the correlation or covariance matrix, as determined by cor [...]" 
    test_cov_ok2 <- isTRUE(all.equal(S_wrong, pca10_internal_cv))    # and yet.. THE DIFFERENCE WILL SHOW UP IN EIGEN
    pca10_internal_edc <- eigen(pca10_internal_cv, symmetric = test_symmetric)                  ## the eigendecomposition is NOT the same IF SYMMETRIC = TRUE
    pca10_internal_edc_2 <- eigen(S_wrong, symmetric = test_symmetric)                          # on paper the covariance matrices are equivalent bu this produce different results if "symmetric = T"...
    test_evec_ok <- all.equal(pca10_internal_edc$vectors, pca10_internal_edc_2$vectors)         # If symmetric = TRUE, this will fail even starting from the same covariance matrix
    test_ev_ok <- all.equal(pca10_internal_edc$values, pca10_internal_edc_2$values)             # eigenvalues unaffected
    ## then fiddles with signs: default = TRUE will generate further discrepancies but can be disabled
    fix <- if (fix_sign){
      function(A) {
        mysign <- function(x) ifelse(x < 0, -1, 1)
        A[] <- apply(A, 2L, function(x) x * mysign(x[1L]))
        A
      }
    } else {
      identity
    }
    pca10_internal_evect <- fix(pca10_internal_edc$vectors)                                      # original code reuses "ev" for both eigenvalues and eigenvectors
    test_evec_sign_ok <- all.equal(pca10_internal_evect, pca10_internal_edc_2$vectors)
  } else {
    ## bypass everything that is wrong (by default)
    test_symmetric = F
    fix_sign = F
    S_right <- (1/(nrow(X)-1))*t(Y)%*%Y  
    #all.equal(S_right, cov(X))
    ## -- NOTE: if eigen decomposition is called, just mentioning "symmetric" regardless of FALSE or TRUE generates sign issues
    pca10_internal_edc <- eigen(S_right)
    pca10_internal_cv <- S_right
    
    ## void the tests
    test_cov_ok1 <- NA
    test_cov_ok2 <- NA
    test_evec_ok <- NA
    test_ev_ok <- NA
  }
  
  ## in any case
  scores = Y %*% pca10_internal_edc$vectors
  
  list(loadings = pca10_internal_edc$vectors,
       eigenvalues = pca10_internal_edc$values,
       sdev = sqrt(pca10_internal_edc$values),
       scores = scores,
       cov_mat = pca10_internal_cv,
       tests = list(test_cov_ok1 = test_cov_ok1,
                    test_cov_ok2 =test_cov_ok2,
                    test_evec_ok =test_evec_ok,
                    test_ev_ok = test_ev_ok)
  )
}

replica_dude <- function(X, original = T, test_symmetric = T){
  ## my version of ade4 dudi.pca
  arg_scale = FALSE                   ## default: TRUE
  df <- X
  df <- as.data.frame(df)
  nf = ncol(df)                       ## default: 2
  nc <- ncol(df)
  if(original){
    test_symmetric <- TRUE            ## default
    ## Part 1 (dudi.pca) centres and normalises -- adapted from View(ade4::dudi.pca)
    row.w = rep(1, nrow(df))/nrow(df)
    col.w = rep(1, ncol(df)) 
    f1 <- function(v) sum(v * row.w)/sum(row.w)
    f2 <- function(v) sqrt(sum(v * v * row.w)/sum(row.w))
    center <- apply(df, 2, f1)
    df <- sweep(df, 2, center)
    if (arg_scale) {
      norm <- apply(df, 2, f2)
      norm[norm < 1e-08] <- 1
      df <- sweep(df, 2, norm, "/")
    } else {
      norm <- rep(1, nc)
    }
    
    ## checkpoint 1: same centred data matrix we would have obtained convenitonally
    all.equal(as.matrix(df), scale(X, center = TRUE, scale = arg_scale)[,])
    
    ## Part 2: variance covariance matrix -- Adapted from View(ade4::as.dudi)
    tol = 1e-07
    df <- as.matrix(df)
    df.ori <- df                                     # backing up data matrix Y
    df <- df * sqrt(row.w)
    df <- sweep(df, 2, sqrt(col.w), "*")
    internal_cov <- df <- crossprod(df, df)                          # covariance matrix to send to eigen
  
    ## checkpoint 2: the covariance matrix that will be decomposed is 1/n Y'Y
    S_Wrong_fromScratch <- 1/nrow(Y)*(t(Y) %*% Y)
    test_cov_orig <- all.equal(as.matrix(df), S_Wrong_fromScratch)
    
    ## Part 3: eigendecomposition
    eig1 <- eigen(internal_cov, symmetric = test_symmetric)
    eig <- eig1$values                                   # eigenvalues
    dval <- sqrt(eig)[1:nf]                              # sqrt(eigenvalues) should be the score's std. dev (assuming Eq.5 holds)
    c1 <- eig1$vectors                                   # I cut short the original version (involves weights but doesn't affect the result)
    li <- data.frame(df.ori %*% c1)                      # scores: Z = Y %*% V
    V <- c1                                           ## loadings matrix
    Z <- li                                           ## scores matrix
                                                      
    ## checkpoint 3: same covariance matrix, yet different eignvectors
    test_eig <- eigen(S_Wrong_fromScratch, symmetric = test_symmetric)
    test_lambda <- test_eig$values
    test_V <- test_eig$vectors
    #all.equal(test_lambda, eig1$values)                    # should be the same
    test_same_cov_diff_load <- all.equal(test_V, c1)        # not the same, even if check point 2 says the covariance matrix is the same!! eigen is overly sensible...                         
    
    ## checkpoint 4: compare with my princomp replica
    test_eig_myprin <- replica_princomp(X, original = TRUE, test_symmetric = test_symmetric, fix_sign = TRUE)
    test_myprin <- all.equal(test_eig_myprin$loadings, eig1$vectors)
    
    
    ## PC Biplot coordinates 
    ## -- unlike base-R biplot, issue with the singular values:  dudi.pca doesn't attempt to get from sqrt(ev) to sv
    ## -- instead, base-R does it, although it wrongly uses sqrt(n)
    biplot_feat_coord <- co <- sweep(c1, 2, dval, "*")                       # this should be B' = DV' but D is made from sqrt(lambda) so it's incorrect
    biplot_obser_coord <- l1 <- sweep(li, 2, dval, "/")                      # this hsould be A = ZD^{-1} but D is incorrect
    
  } else {
    Y = scale(X, center = TRUE, scale = FALSE)
    df <- as.matrix(Y)
    #test_symmetric <- FALSE
    S_right <- (1/(nrow(X)-1))*t(Y)%*%Y  
    #all.equal(S_right, cov(Y))
    internal_cov <- S_right
    
    ## The following "repetition" is needed toreplicate my from scratch results
    ## -- if eigen decomposition is called, just mentioning "symmetric" regardless of FALSE or TRUE generates sign issues
    #eig1 <- eigen(internal_cov, symmetric = test_symmetric) 
    eig1 <- eigen(internal_cov) 
    V <- eig1$vectors
    Z <- Y %*% V
    dval <- sqrt(eig1$values)*sqrt(nrow(X)-1)                                # Eq. 7 in my paper
    
    # Principal component biplot coordinates are computed here as per Gabriel 1971
    ## -- issue: 
    biplot_feat_coord <- sweep(V, 2, dval, "*")
    biplot_obser_coord <-  sweep(Z, 2, dval, "/")
    
    test_cov_orig <- NA
    test_myprin <- NA
    test_same_cov_diff_load <- NA
  }
  
  
  ## output
  list(Z = Z,
       V = V,
       ev = eig1$values,
       sv = dval,
       internal_cov = internal_cov,
       biplot_obser_coord = biplot_obser_coord,
       biplot_feat_coord = biplot_feat_coord,
       test_myprin = test_myprin,
       test_cov_orig = test_cov_orig,
       test_same_cov_diff_load = test_same_cov_diff_load
       )
}

replica_acp <- function(X, original = T, test_symmertic = FALSE){
  center_arg <- T                                               # default
  reduce_arg <- F                                               # default is TRUE
  if(original){
    ## A) center and normalise first. Uses cols/row weights;
    X <- as.matrix(X)
    wI = rep(1, nrow(X))                                      # rows weights set to 1
    wV = rep(1, ncol(X))                                      # cols weights set to 1
    if (center_arg) 
      Y <- t(t(X) - as.vector((wI %*% X)/sum(wI)))
    if (reduce_arg) 
      Y <- apply(Y, 2, function(u) {
        u/stats::sd(u)
      })
    ## checkpoint 1: equivalent to usual centering and scaling
    Y_checks <- isTRUE(all.equal(Y, scale(X, center = center_arg, scale = reduce_arg)[,]))
    
    ## B) covariance matrix
    S <- (t(Y) * wI) %*% (Y * wV)                              # original code, using weights
    S_equiv <- (t(Y)%*%Y)                                      # covariance matrix as I understand it
    ## checkpoint 2
    S_checks <- isTRUE(all.equal(S, S_equiv))
    
    ## C) eigendecomposition
    eig_pca6 <- eigen(S, symmetric = test_symmertic)           # "symmetric" = F set by default and cannot be changed
    eig_pca6_alt <- eigen(S_equiv, symmetric = test_symmertic)
    lambda <- eig_pca6$values
    V <- eig_pca6$vectors
    ## checkpoint 3
    V_checks <- isTRUE(all.equal(eig_pca6$vectors, eig_pca6_alt$vectors))    # check if different cov matrix equivalent and if symmetric argument (even if F) will impact results 
    
    ## D) mess up with signs (original only)
    V <- V %*% diag(sign(lambda))                              # signs of V "fixed" using the signs of the eigenvalues
  } else {
    Y <- scale(X, center = T, scale = F)[,]
    S <- (1/(nrow(X)-1))*t(Y)%*%Y  
    eig_pca6 <- eigen(S)
    lambda <- eig_pca6$values
    V <- eig_pca6$vectors
    ## checks not needed
    Y_checks <- NULL
    S_checks <- NULL
    V_checks <- NULL
  }
  Z = Y %*% V
  sdev_Z <- apply(Z, 2, sd)
  list(Z = Z,
       V = V,
       lambda = lambda,
       cov_mat = S,
       Y_checks = Y_checks,
       S_checks = S_checks,
       V_checks = V_checks
        )
}

replica_svdTriplet <- function(Y, original = T, sign_leave = F){
  ## taken from View(FactoMineR::svd.triplet)
  ## Main outcome: D and U will be respectively mutliplied and divided by sqrt(n)
  if(original){
    row.w <- rep(1/nrow(Y), nrow(Y))
    row.w <- row.w/sum(row.w)
    col.w <- rep(1, ncol(Y))
    Y_weigh <- sweep(Y, 2, sqrt(col.w), "*")
    Y_weigh <- sweep(Y, 1, sqrt(row.w), "*")
    #Y_weigh <- t(t(Y) * sqrt(col.w)) * sqrt(row.w)   # original formulation
    
    ## checkpoint 1: the matrix being decomposed is just Y / sqrt(n)
    test_Y_weigh <- all.equal((1/sqrt(nrow(Y)))*Y, Y_weigh)
  } else {
    Y_weigh <- Y                                      # do nothing
  }
  #svd.usuelle <- svd(Y_weigh, nu = ncp, nv = ncp)    # original formulation
  svd.usuelle <- svd(Y_weigh)                         # turns out this "generalised" svd mean svd of Y / sqrt(n). package Vegan does something simliar too
  U <- svd.usuelle$u
  V <- svd.usuelle$v
  vs <- svd.usuelle$d
  #if (ncp > 1) {
  if(original){
    if(!sign_leave){
      ## messes around with sings. I give the possibiilty to switch this off
      mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
      mult[mult == 0] <- 1
    } else {
      mult <-  rep(1, nrow(V))
    }
    U <- t(t(U) * mult)
    V <- t(t(V) * mult)
    U <- U/sqrt(row.w)
    V <- V/sqrt(col.w)
    #vs <- svd.usuelle$d[1:min(ncol(Y_weigh), nrow(Y_weigh) - 1)]
  }
  #}

  ## out
  list(U = U, 
       V = V, 
       sv = vs, 
       test_Y_weigh = test_Y_weigh)
}

replica_psych <- function(Y, original = T, signs_thing = T){
  ## adapted from View(psych::principal)
  mess_w_signs <- function(loadings){
    ## For replicability I parked part of the original code in this auxiliary function, which modifies the signs
    sign.tot <- vector(mode = "numeric", length = length(eigens$values))
    sign.tot <- sign(colSums(loadings))
    sign.tot[sign.tot == 0] <- 1
    loadings <- loadings %*% diag(sign.tot)
    ## the followign bit  will mess things about beyond the signs so I'm switching it off
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated, decreasing = TRUE)
    #loadings <- loadings[, ev.order]                      
    signed <- sign(colSums(loadings))
    signed[signed == 0] <- 1
    loadings <- loadings %*% diag(signed)
    loadings
  }
  cor = "cov"                                  # default = "cor"
  rotate = "none"                              # more of a factor analysis thing
  ## A) compute covariance matrix
  S <- cov(Y, use = "pairwise")                # MUST USE THIS. looks like it won't make any difference but it DOES when using eigen
  sd_feat <- sqrt(diag(S))                     # this is the St.Dev. of the FEATURES, I use it on LHS of Eq.14
  eigens <- eigen(S)
  values <- eigens$values
  if(original){
    eigens$values[eigens$values < .Machine$double.eps] <- .Machine$double.eps
    ## B) Get the "weighted" version of the loadings
    #loadings <- loadings %*% sqrt(diag(eigens$values))
    loadings <- eigens$vectors %*% sqrt(diag(eigens$values, nrow = length(eigens$values)))  
    if(signs_thing){
      ## C) mess about with signs. not optional in the original script
      loadings <- mess_w_signs(loadings)
      evec_4replica <- mess_w_signs(eigens$vectors)             # added by me for reproducibility
    }
    
    ## checkpoint
    isEqualUnname(loadings, test_pca7$loadings[,])
    
    ## D) get the scores in the original script
    Structure <- loadings
    weights <- try(solve(S, Structure), silent = TRUE)          # this should be comparable with eigenvectors / sqrt(eigenvalues)
    scores <- scale(Y, scale = FALSE) %*%  weights              # the first term is just Y
    
    ## My alternative to get the same output but requires that eigenvectors' signs are messed about via the function mess_w_signs()
    scores2 <- sweep(Y %*% evec_4replica, 2, sqrt(values), "/")       
    loadings2 <- sweep(evec_4replica, 2, sqrt(values), "*")
    
    ## checkpoints
    test_scores <- isEqualUnname(scores2, scores)
    test_loadings <- isEqualUnname(loadings2, loadings)
    test_weigths <- isEqualUnname(sweep(evec_4replica, 2, sqrt(values), "/") , weights)  # under the assumptions we've messed about the signs of the loadings
    
    ## biplot coordinates (not in the original script)
    biplot_feat_coord = loadings
    biplot_observ_coord = scores
  } else {
    ## Reset to classic
    loadings <- eigens$vectors
    scores <- Y %*% loadings
    biplot_feat_coord <- sweep(loadings, 2, sqrt(values), "*")            # keeping the weighted version
    biplot_observ_coord  <- sweep(Y %*% loadings, 2, sqrt(values), "/")   # keeping the weighted version
    
    ## set tests to null
    test_scores <- NULL
    test_loadings <- NULL
    test_weigths <- NULL
    
  }
  ## output
  list (scores = scores,
        loadings = loadings,
        test_scores =  test_scores,
        values = values,
        sd_feat = sd_feat,                                                  # not an output in the original script
        biplot_feat_coord = biplot_feat_coord,                              # not an output in the original script
        biplot_observ_coord = biplot_observ_coord,                          # not an output in the original script
        cov_mat = S,                                                        # not an output in the original script
        test_scores = test_scores,
        test_loadings =test_loadings,
        test_weigths = test_weigths
        )
}

replica_svd_vegan <- function(X){
  ## To understand pca done via "rda" you have to "drill down" through several scripts and inspect the following scripts in sequence
  ## --- a) View(vegan:::rda.default
  ## --- b) View(vegan:::ordConstrained)
  ## --- c) View(vegan:::initPCA)
  ## --- d) View(vegan:::ordConstrain)
  Y_veg <- scale(X, scale = FALSE)
  Y_veg <- Y_veg/sqrt(nrow(Y_veg) - 1)      # happens in vegan:::initPCA notice weighting of Y is similar to (but not the same as) svd.triplet in FactoMiner
  svd_veg <- svd(Y_veg)                     # this is just me cutting to the chase. Details in vegan:::ordConstrain
  D_veg <- diag(svd_veg$d)                  # this is not recorded by vegan. 
  eig <- svd_veg$d^2                        # this happens to be actually equal to the eigenvalues of S 
  
  ## output
  list(u = svd_veg$u,
       v = svd_veg$v,
       eig = eig,
       d = svd_veg$d,
       D = D_veg)
}

replica_vegan_scores_FAIL <- function(test_pca_x1){
  ## COULD NOT FULLY REPLICATE the mystery constant  (my take on some excerpts of the relevant code below but not quite there yet)
  ## based on View(vegan:::scores.rda)
  eigval <- test_pca_x1$CA$eig                                   # x$sdev^2
  sumev <- test_pca_x1$CA$tot.chi                                # sum(eigenvalues)
  constant <- sqrt(sqrt((nrow(test_pca_x1$CA$u) - 1) * sumev))   # IMPORTANT FOR LATER: MYSTERY CONSTANT
  const <- c(constant, constant)
  wa <- test_pca_x1$CA$u
  v <- test_pca_x1$CA$v
  slam <- sqrt(eigval/sumev)
  
  v2 <- sweep(v, 2, sqrt(slam), "*")
  v1 <- v2 * sqrt(sumev/(nrow(test_pca_x1$CA$u) - 1))
  v <- const[1] * v1
  species <- v
 
  wa1 <- sweep(wa, 2, sqrt(slam), "*")
  wa <- const[2] * wa1
  
  ## output
  list(sites_replica = wa,
       species_replica = species,
       costant = constant,
       sites_pre_const = wa1,
       species_pre_cosnt = v1)
} 

replica_svd_PCAmix <- function(X, original = T, sign_leave = F){
  ## for comparison
  we_veg <- 1/sqrt(nrow(X) - 1)                    # the weight used in vegan, to show the difference
  Y_vegan <- scale(X, scale = F)[,]*we_veg
  
  ## reproduces svd triplet for PCAmix
  we <- sqrt(nrow(X)/(nrow(X) - 1))                # different from both Vegan and FactoMiner pkgs
   if(original){
    myscale = apply(X, 2, sd)                      # unlike other packages, PCA mix forces normalisation
  } else {
    myscale = FALSE                                # The following fix is not available in the original package
  }
  ## for output
  Y <- scale(X, scale = myscale)[,]
  Y_w <- Y*we
  
  ## Approach 1: call "generalised" SVD as the original function does
  svd_facto <- replica_svdTriplet(Y_w, sign_leave = sign_leave)       # call my replica of svd triplet i.e. the matrix being decomposed is now Y_w / sqrt(n)= Y / sqrt(n-1) which is like vegan     
  D <- diag(svd_facto$sv)
  
  ## Approach 2: equivalent but in the style of package "vegan" 
  ## performs a conventional svd but on a "revised" weighted matrix, knowing what "svd.triplet" actually does
  Y_w2 <- Y_w / sqrt(nrow(Y_w))  # revised weighted matrix
  svd_veg <- svd(Y_w2 )
  
  ## recall that svd.triplet messes about
  U_triplet_to_vegan <- svd_facto$U/sqrt(nrow(Y))
  
  ## checks equivalences (for debug)
   if(!original){ all.equal(Y_w2,  Y_vegan)      # if not standarised should be the same weighted matrix used in vegan
     }
   isEqualUnname(U_triplet_to_vegan, svd_veg$u)
   isEqualUnname(svd_facto$V, svd_veg$v)
   isEqualUnname(svd_facto$sv, svd_veg$d)
  
  ## output
  list(Y_weigh = Y_w,            # returned as Z
       u = svd_facto$U,          #  U will be divided by sqrt(n)
       v = svd_facto$V,
       d = svd_facto$sv,         #  D will be divided by sqrt(n-1)
       D = D,
       svd_veg = svd_veg)
}

replica_biplotPCAtools <- function (pcaobj, showLoadings = TRUE){
  ## selected lines from View(PCAtools::biplot)
  x = "PC1"
  y = "PC2"
  lengthLoadingsArrowsFactor = 1.5                         # it's a constant
  plotobj <- NULL
  
  ## Observations coordinates
  plotobj$x <- pcaobj$rotated[, x]
  plotobj$y <- pcaobj$rotated[, y]
  plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)
  
  ## Features coordinates (in the original code they're conditional on users selecting "showLoadings")
  if(showLoadings){
    xidx <- order(abs(pcaobj$loadings[, x]), decreasing = TRUE)
    yidx <- order(abs(pcaobj$loadings[, y]), decreasing = TRUE)
    r <- min((max(pcaobj$rotated[, x]) - min(pcaobj$rotated[, x])/(max(pcaobj$loadings[, x]) - min(pcaobj$loadings[, x]))), 
             (max(pcaobj$rotated[, y]) - min(pcaobj$rotated[, y])/(max(pcaobj$loadings[, y]) - min(pcaobj$loadings[,y]))))
    xend = pcaobj$loadings[, x] * r * lengthLoadingsArrowsFactor
    yend = pcaobj$loadings[, y] * r * lengthLoadingsArrowsFactor
  }
  feat_coord <- cbind(xend, yend)
  rownames(feat_coord) <- rownames(pcaobj$loadings)
  
  ## output
  list(observ_coord = plotobj, 
       feat_coord = pcaobj$loadings[, c(x, y)],
       feat_coord_scales = feat_coord)
}

replica_biplotEZ <- function(X, correlation.biplot = T, original = T){
  ## selected code from package biplotEZ, based on 
  ## --- https://rdrr.io/cran/biplotEZ/src/R/PCA.R
  ## --- https://rdrr.io/cran/biplotEZ/src/R/biplot.R
  ## --- View(biplotEZ:::PCA.biplot)
  ## underlying math: https://cran.r-project.org/web/packages/biplotEZ/vignettes/biplotEZ.html
  center = TRUE                             # default
  scaled = FALSE                            # default
  Y_EZ <- scale(X, scale = scaled)          # interanlly calls "biplot" to centre and scale: this is my shortcut
  e.vects = 1:ncol(Y_EZ)                    # n of PCs to retaub
  
  ## performs svd
  svd.out <- svd(Y_EZ)
  U.mat <- svd.out$u
  Sigma.mat <- diag(svd.out$d)                             # this is D
  V.mat <- svd.out$v
  Lmat <- svd.out$v                                        # duplications, but unlike "V", Lmat is reported as output keeping all dimension
  Vr <- svd.out$v[, e.vects, drop = FALSE]                 # in principle, coordinates of the VARIABLES see https://cran.r-project.org/web/packages/biplotEZ/vignettes/biplotEZ.html
  eigenvalues <- svd.out$d^2                               # incorrect based on Eq.7 in my paper
 
  ## PC biplot vs scores and loadings biplot
  if (correlation.biplot){                                 # equivalent to is what I call PC biplot
    if(original){                                          # reproduces with what I think are MISTAKES
      ## Note: the original formula is busy and I break it down, but the logic is as follows
      ## -- A = U =ZD^{-1} assuming Z = YV which leads to A = YVD^{-1} 
      ## -- replacing D = diag(sv) = sqrt(n-1)*diag(sqrt(eigenvalues of S)) where S = Y'Y 
      
      ## step 1) 
      S_EZ_replica <- t(Y_EZ) %*% Y_EZ                                                                     # a kind of covariance matrix
      lambda_EZ_replica <- svd(S_EZ_replica)$d                                                             # seems to (incorrectly) treat these as "eigenvalues" of S 
      D_EZ_replica <-  diag(sqrt(lambda_EZ_replica))/sqrt(nrow(Y_EZ)-1)                                    # incorrect but preserves the original logic: when isolated, this step disagrees with Eq.7 
      
      ## step 2
      Z_EZ_replica <- Y_EZ %*% V.mat                                                                        # scores are not returned as such in the original code. The original approach uses Vr to limit to 2 dimensions
      
      ## step 3: correct replication of incorrect PC biplot coordinates
      EZ_biplot_obscoord_replica <- Z_EZ_replica %*% solve(D_EZ_replica)                                    # returned as Z by the original function. Actually it's A = ZD^{-1} except D is not diag(sv)
      EZ_biplot_featcoord_replica <- V.mat %*% solve(D_EZ_replica)                                          # returned as Lmat by the original function. Violates Eq.12 i.e., B = VD not B = VD^{-1} 
    } else {
      ## now with "fixes" / from scratch understanding
      lambda_EZ_replica <- svd(Y_EZ)$d^2/(nrow(Y_EZ) - 1)                                                   # added by me: correct relationship between ev and sv based on my Eq.7
      D_EZ_replica <- Sigma.mat                                                                             # what I expect to be sv
      
      ## checkpoint (optional)
      #isEqualUnname(my_D, D_EZ_replica)
      #isEqualUnname(lambda_EZ_replica, my_lambda_svd) 
      
      Z_EZ_replica <- Y_EZ %*% V.mat    
      EZ_biplot_obscoord_replica <- Z_EZ_replica %*% solve(D_EZ_replica)                                           # Eq.11
      EZ_biplot_featcoord_replica <-  D_EZ_replica %*% t(V.mat)                                                    # Eq.12
      
      S_EZ_replica <- NULL
    } 
    ## Coordinates passed to plotting function: View(biplotEZ:::plot.biplot)
    ## I omit the axis calibration as Gower 2011 orignally does not do it for PC biplots as understood here
    ax.one.unit <- NULL
  } else {
    ## scores and loadings biplot
    Z_EZ_replica <- Y_EZ %*% V.mat
    EZ_biplot_obscoord_replica <- Z_EZ_replica                                                              # returned as Z by the original function
    EZ_biplot_featcoord_replica <- V.mat                                                                    # returned as Vr or Lmat by the original function
    
    S_EZ_replica <- NULL
    lambda_EZ_replica <- NULL
    
    ## Biplot axes calibration as described by Gower for a scores and loadings biplt
    ## See Gower book: "calibrated axis biplot" uses some kind of reconstruction / projection of the points onto the variables "axes"
    Xhat <- Z_EZ_replica %*% t(Lmat[,e.vects])
    ax.one.unit <- 1/(diag(Vr %*% t(Vr))) * Vr                                                              # "one unit in the positive direction of the biplot axis." accroding to https://cran.r-project.org/web/packages/biplotEZ/vignettes/biplotEZ.html
  }
  
  ## output list
  list(biplot_obscoord = EZ_biplot_obscoord_replica,
       biplot_featcoord = EZ_biplot_featcoord_replica,
       S_EZ_replica  = S_EZ_replica,
       D_EZ_replica = D_EZ_replica,
       lambda_EZ_replica = lambda_EZ_replica,
       Z_EZ_replica = Z_EZ_replica,
       biplot_axis_unit = ax.one.unit
       )
}

replica_biplotGUI <- function(X, predict = TRUE, mysymmetric = TRUE, signs_thing = T){
  ## based on View(BiplotGUI::Biplots)
  ## extensive descripiton at https://www.jstatsoft.org/article/view/v030i12
  ## The code is difficult to navigate and poorly annotated hence I rely most on Appendix A.4 of the paper
  ## Below I focus on what they cal "covariance / correlation biplot" specifically covariance biplot see JSS paper A.4 appendix
  ## based on functions (each with a menu entry in the GUI)
  ## -- "Joint.CovarianceCorrelation.determine"
  
  ## Caveats for comparison
  ## -- similar to EZ biplots (possibly some of the same mistakes)
  ## -- The package is trying to obtain A = ZD^{-1} see eq.11 assuming Z = YV so that A = YVD^{-1} 
  ## -- in principle, since D = diag(sv) = sqrt(n-1)*diag(sqrt(eigenvalues of S)). The original formula crams all of this together with some imprecisions
  
  
  ## A) centre data matrix
  n <- nrow(X)
  p <- ncol(X)
  Biplot.Xtransformed <- scale(X, scale = F)[,]
  isEqualUnname(Biplot.Xtransformed, Y)
  
  ## B) eigen-decomposition of cov matrix
  S_guy <- t(Biplot.Xtransformed) %*% Biplot.Xtransformed
  temp1 <- eigen(S_guy, symmetric = mysymmetric)                                              ## eigendecomposition. symmetric = T default
  eigenval <- temp1$values
  if(signs_thing){                                                                            ## not optional in original script
    temp1$vectors <- (apply(temp1$vectors, 2, function(x) x * sign(x[which.max(abs(x))])))    ## messess with signs
  }
  V_guy <- temp1$vectors                                # V       loadings i.e. features coordinates
  Z_guy <- Biplot.Xtransformed %*% V_guy                # Z = YV  Scores
  
  ## C) "scores and loadings biplot", which you get under "Joint.PCA.determine"      
  guy_biplot_PCA_feat_coord <- V_guy                    # should match GUI export output: "Basis vectors", joint PCA
  guy_biplot_PCA_observ_coord <- Z_guy                  # should match GUI export output: "matrix of point coordinates", Joint PCA
  
  ## D) "PC biplot", which you get for export under Joint.CovarianceCorrelation.determine
  ## -- WARNING: here I replicate what I think is an incorrect result
  ## -- fixes provided later for information
  sv_guy_wrong <- sqrt(eigenval)/sqrt(nrow(Y_EZ)-1)                              # same mistake as EZ biplot
  sv_guy_wronger <-  sqrt(eigenval)/sqrt(nrow(Y_EZ))
  D_guy_wrong <- diag(sv_guy_wrong)
  D_guy_wronger <- diag(sv_guy_wronger)
  A_guy_wrong <- Z_guy %*% solve(D_guy_wrong)                                    # in original code: Biplot.Y3D_
  B_guy_wrong_interpolate <- V_guy %*% solve(D_guy_wrong)                        # in original code: Biplot.Binterpolate_ regardless of "predict"
  if (predict){
    ## alternative values for Biplot.B3D_
    B_guy_wrong <- V_guy %*% D_guy_wronger 
  } else {
    B_guy_wrong <- B_guy_wrong_interpolate                                  
  }
  guy_biplot_Cov_feat_coord <- A_guy_wrong      #  should match GUI export output: "matrix of point coordinates", covariance/correlation 
  guy_biplot_Cov_observ_coord <-  B_guy_wrong   #  should match GUI export output: "Basis vectors", covariance/correlation
  
  ## -- For information: these are suggested fixes
  sv_guy_fixed <- sqrt(eigenval)*sqrt(nrow(Y_EZ)-1)                              # based on Eq.7
  D_guy_fixed <-  diag(sv_guy_fixed)                    
  A_guy_fixed <- Z_guy %*% solve(D_guy_fixed)                                    # based on Eq.11 A = ZD^{-1} Eq.11, the coord of observ in a PC biplot    
  B_guy_fixed <- V_guy %*% D_guy                                                 # based on Eq.12 
  
  ## E) Axis calibration? PLEASE CHECK GOWER 2011
  Xhat <- Biplot.Xtransformed %*% Biplot.B_ %*% t(Biplot.B_)                        # projection
  Xhat1 <- Biplot.Xtransformed %*% Biplot.B_[, 1] %*% t(Biplot.B_[, 1])
  PointsTab.predictivities <- rowSums(Xhat^2)/rowSums(Biplot.Xtransformed^2)
  PointsTab.predictivities1dim <- rowSums(Xhat1^2)/rowSums(Biplot.Xtransformed^2)
  AxesTab.predictivities <- colSums(Xhat^2)/colSums(Biplot.Xtransformed^2)
  AxesTab.predictivities1dim <- colSums(Xhat1^2)/colSums(Biplot.Xtransformed^2)
  
  ## ouput
  list(biplot_obscoord_pc = A_guy_wrong,         
       biplot_featcoord_pc = B_guy_wrong,       
       S_guy  = S_guy,
       V_guy = V_guy,
       Z_guy = Z_guy,
       eigenval = eigenval,
       sv_guy_wrong = sv_guy_wrong,
       sv_guy_wronger = sv_guy_wronger
  )
}

replica_bpca <- function(X, original = T){
  ## package bpca
  ## View(bpca:::bpca.default)  & ?bpca:::bpca.default & https://rdrr.io/cran/bpca/src/R/bpca.default.R & 
  ## View(bpca:::bpca.prcomph & https://rdrr.io/cran/bpca/src/R/bpca.prcomp.R
  x <- as.matrix(X)
  scale = FALSE                               # default is TRUE
  x.cent <- sweep(x, 2, apply(x, 2, mean))    # default, centred data matrix 
  isEqualUnname(x.cent, Y)
  if (scale){
    x.scal <- sweep(x.cent, 2, apply(x.cent, 2, sd), "/")
  } else {
    x.scal <- x.cent
  } 
  
  ## part 1: perfrom svd
  svdx.scal <- svd(x.scal)
  
  ## unfortunately they get the eigenvalues wrong and so the variance explained
  eigenvalues = svdx.scal$d 
  
  ## part 2: PCA biplot a la Gabriel (it respects eq. 8 in our paper but swaps the multiplication by the constant sqrt(nrow(x) - 1))
  ## NOTICE: this has consequences: the correlation coefficients won't be correct.... i.e. equation in my paper 18 won't hold
  ## corresponds to argument: method = "gh": this is the PCA biplot as intended by Gabriel
  s2.scal <- diag(svdx.scal$d)                                     # D in our notation
  
  if(original){
    g.scal <- sqrt(nrow(x) - 1) * svdx.scal$u                        # Eq.52 in Gabriel 1971? our A = U but times sqrt(nrow(x) - 1)
    h.scal <- 1/sqrt(nrow(x) - 1) * s2.scal %*% t(svdx.scal$v)       # Eq.52 in Gabriel 1971? our B^T = DV^T but divided by sqrt(nrow(x) - 1)
  } else {
    g.scal <- svdx.scal$u                                            # A = U 
    h.scal <- s2.scal %*% t(svdx.scal$v)                             # B^T = DV^T
  }
  
  #isEqualUnname(h.scal, B_trsp/sqrt(nrow(Y)-1))                        # check...
  hl.scal <- t(h.scal)                                             # from B' to B, however our B is now divided by sqrt(nrow(x) - 1)
  pc.names <- paste("PC", 1:ncol(hl.scal), sep = "")
  
  ## Part 3: correlation coefficients
  ## -- from function var.rbf()
  d = 1:2
    ## AS PER SCRIPT
    x <- hl.scal[, d[1]:d[length(d)]]                             # matrix B select as many COLUMNS as there are axis in the biplot
    ## HOWEVER, TO MAKE THE CORRELATION WORK WE NEED TO SET IT AS MY MODEL
    #x <- hl.scal*sqrt(nrow(Y)-1)                                  # for comparison with our matrix B
    lv <- function(x) sqrt(t(x) %*% x)                            # 2-norm of vector
    l <- apply(x, 1, lv)                                          # 2-norm of each ROW: compare eq. 15 in my paper
    n <- nrow(x)
    var.rb <- diag(n)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        cost <- (t(x[i, ]) %*% x[j, ])/(l[i] * l[j])              # cosine 
        var.rb[j, i] <- cost
        var.rb[i, j] <- var.rb[j, i]
      }
    }
    dimnames(var.rb) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
    var.rb.res <- var.rb
  
  ## this is the correlations coeffs using ALL PCS ...BUT ONLY WORKS IF x <- hl.scal*sqrt(nrow(Y)-1) i.e. my version  !!!!
    vector_varrbres <- c(var.rb.res[lower.tri(var.rb.res)])
    #all.equal(unname(vector_varrbres), unname(cos_theta_ab_allPCs))
    #all.equal(unname(vector_varrbres), unname(correl_feat_fromScratch))
    
  ## part 3: visualise
    list(coord = list(objects = g.scal, variables = hl.scal),
         correl = vector_varrbres)
}

replica_ggbiplot <- function(pcobj, pcbiplot_alpha = 1, pcbiplot_option = T){
  ## step 1) based on View(ggbiplot::get_SVD)
  ## Based on our review issue arise already here
  if (inherits(pcobj, "prcomp")) {
    n <- nrow(pcobj$x)
    D <- pcobj$sdev                                           ## incorrect notation suggests these are  singular values from SVD; but here sdev is just sqrt(eigenvalues [of S])  or equivalently sqrt(variance of scores)
    U <- sweep(pcobj$x, 2, 1/(D * sqrt(n)), FUN = "*")        ## Assumes Z = UD -> U = ZD^{-1}; also aware that D is sqrt(eigenvalues), not sv hence follows  Eq.7 in our paper and assumes: D*sqrt(n) =  singular values -- although it should be sqrt(n-1)
    V <- pcobj$rotation                                       ## retrieve loadings
  }
  if (inherits(pcobj, "princomp")) {
    n <- pcobj$n.obs
    D <- pcobj$sdev                                           ## same notation problem: this are not the singular values
    U <- sweep(pcobj$scores, 2, 1/(D * sqrt(n)), FUN = "*")   ## again Z = UD -> U = ZD^{-1} but since D = sqrt(eigenvalues), assumes D*sqrt(n) =  singular values
    V <- pcobj$loadings
  }
  if (inherits(pcobj, "pca") & inherits(pcobj, "dudi")) {
    n <- nrow(pcobj$tab)
    D <- sqrt(pcobj$eig)                                      ## this is incorrect. dval, not eig
    U <- pcobj$li                                             ## seems incorrect, this is the scores matrix 
    V <- pcobj$co                                             ## this seems incorrect: c1 is V
  }
  
  ## step 2) based on  View(ggbiplot::ggbiplot); please also check https://rdrr.io/cran/ggbiplot/src/R/ggbiplot.r
  ## -- notation suggests the assumption is svd = UDV' but as shown above D is not s.v.
  ## -- default values
  choices = 1:ncol(U)                                          ## default 1:2
  pc.biplot =  pcbiplot_option                                 ## principal component biplot?
  scale = pcbiplot_alpha                                       ## default: 1; similar to alpha in my paper
  obs.scale = 1 - scale
  var.scale = scale
  var.factor = 1
  
  ## just relabels
  d <- D
  u <- U                                                      
  v <- V
  
  ## if scale = 1 then we get a PC biplot else we get a scores and loadings biplot
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))  ## if obs.scale =1 this goes back to the scoring matrix since UD = Z (with all the caveats); if obs.scale = 0 then we get a PC biplot since  A = U
  v <- sweep(v, 2, d^var.scale, "*")                                              ## if obs.scale =1 then var.scale = 0 and we get B = V; otherwise we get B = VD (with all the caveats)
  df.v <- var.factor * v
  df.v <- as.data.frame(v[, choices])
  df.v <- var.factor * df.v                                                       
  #names(df.u) <- c("xvar", "yvar")
  #names(df.v) <- names(df.u)
  if (pc.biplot) {                                                                ## yes is default. From documentation, this is trying Gabriel Eq.52 but not sure it's correct
    nobs.factor <- ifelse(inherits(pcobj, "prcomp"), sqrt(n - 1), sqrt(n))        ## I suppose because prcomp is the only one that gets the sample variance matrix assumption right...
    df.u <- df.u * nobs.factor                                                    ## shouldn't this affect both u and v?
  }
  
  ## additional manipulations on the features coordinates (commentary on https://rdrr.io/cran/ggbiplot/src/R/ggbiplot.r but not on source code)
  circle.prob = 0.68
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v_2 <- r * df.v/sqrt(max(v.scale))
  
  ## output
  list(svd = list(D = d,
                  U = u,
                  V = v),
       observ_coord = df.u, 
       feat_coord = df.v,
       feat_coord_processed = df.v_2)
}

replica_multBiplotR <- function(Y, alpha){
  ## PCA.Biplot in  https://cran.r-project.org/web/packages/MultBiplotR/MultBiplotR.pdf
  ## based on 
  ## ---View(MultBiplotR::PCA.Analysis) and View(PCA.Biplot)
  ## ---https://rdrr.io/github/villardon/MultBiplotR/src/R/PCA.Biplot.R
  # Calculating the Biplot
  
  dimension = ncol(Y)                                       # default: 3
  #Scaling = 5                                              # I ignore and start from Y which corresponds to 4 but defualt is 5
  #alpha = 1                                                 # default is 1 and I think this should be the PCA biplot. Compare with bpca package, same options?
  
  n = nrow(X)
  p = ncol(X)
  
  ## SVD as usual
  ## -- warning: oddly enough the svd is done on the data matrix. For comparison I am using Y
  SD = svd(Y, nu = dimension, nv = dimension)
  EV = SD$d^2                                               # wrong, these are not eigenvalues
  
  a = SD$u %*% diag(SD$d[1:dimension])                      #  scores matrix Z = UD see eq. 10 in my paper, not A as in eq.11 in my paper 
  b = SD$v %*% diag(SD$d[1:dimension])                      # this should be matrix B according to eq. 12 of my paper
  Inertia = round((EV/sum(EV)) * 100, digits = 3)           # same language as French school e.g. Factominer
  CorrXCP = cor(Y, a, method = "pearson")                   # shortcut for the correlation coefficients
  
  # I skip a bunch of stuff that remind me of FactoMiner - cannot be bothered
  
  ## Biplot coordinates: correct in principle
  a = a %*% diag((1/SD$d[1:dimension])^(1 - alpha))       # if alpha = 0 eq 11: A = ZD^{-1} if alpha = 1 A = Z
  b = b %*% diag((1/SD$d[1:dimension])^alpha)             # if alpha = 0 eq 12: B = VD; if alpha = 1 B = V
 
  ## then they do this and messes up: is this lambda scaling as in Gower 2011?
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/n
  scb = scb/p
  scf = sqrt(sqrt(scb/sca))
  wrong_a = a * scf
  wrong_b = b / scf
  
  list(EigenValues = EV,
       sv = SD$d,                    # normally not returned
       EV = SD$v,
       ## coordinates
       RowCoordinates = wrong_a,
       ColCoordinates = wrong_b,
       correct_coord_observ = a,
       correct_coord_feat = b,
       Structure = CorrXCP,
       Inertia = Inertia)
}



replica_bipl5 <- function(X_extd){
  ## This function is meant to replicate relevant functions from the package  bipl5
  ## It is my own streamlined reading of what is relevant for my analysis and is not meant to just reproduce the original scripts
  ## The remainder is based on these sources for the original script: 
  ## -- View(bipl5::PCAbiplot)
  ## -- https://rdrr.io/cran/bipl5/src/R/classbipl5.R
  ## -- https://rdrr.io/cran/bipl5/src/R/PCAbiplot_Helper.R
  
  ## Part 0 feature's means and sd                        ####
  n <- nrow(X_extd)
  p <- ncol(X_extd)
  mu <- colMeans(X_extd)                                  # mean for each feature
  all.equal(apply(X_extd,2,mean), colMeans(X_extd))       # equivalent using apply
  stddev <- apply(X_extd,2,sd)                            # st.dev for each feature
  scale_spec = FALSE                                      # default true
  if (scale_spec) {
    stddev <- rep(1, p)                                   # as in original script: sets the st.dev to be 1 but does not change the mean to zero...weird
  }
  
  ## Part 1 calls SVD on centred or standardised data     #####
  Y <- scale(X_extd, scale = scale_spec)[,]               # built-in equivalent of centering or standardising
  PCA <- svd(Y)                                           # eq. 6 in my paper
  V.mat <- PCA$v
  U.mat <- PCA$u
  stddev.mat <- diag(PCA$d)                               # WRONG: not a standard deviation (...of ??)
  eigval <- PCA$d^2                                       # WRONG: squaring the singular values is not equal to the eigenvalues (...of var-cov mat of Y)
  lambda.mat <- diag(eigval)
  
  ## impose a rank 2-approx
  basis = 1:2                                             # by default two components: therefore it's a rank-2 approx see eq. 9 in my paper
  D <- diag(PCA$d)[basis, basis]                              
  U <- PCA$u[, basis]
  V <- PCA$v[, basis]                                     # calls this "Matrix of vector loadings from SVD"
  Z <- U %*% D                                            # eq. 10 in my paper. PCA scores from SVD...but from a rank-2 approximation
  
  ## warning when using a rank-2 approx
  back_to_Y <- U.mat %*% stddev.mat %*% t(V.mat)
  approx_Y <- Z %*% t(V)
  
  ## Part 2 creates the biplot from loadings              ####
  ## -- the following is what I undersetand from the function "make_biplot" 
  ## -- note: in the original argument "pc12" means "the first two principal components" but really the input is a big list of stuff called "x" (just like the raw data) and created by the "construct_biplot" function, 
  ## -- original call (overloads the "x" object): 
  ## ---- biplot_details<-make_biplot(x, colorpalette=NULL, symbol = "circle")
  ## ------- which in turns calls add_Vector_biplot: addpc12 <- add_Vector_biplot(x = pc12)
  ## ------- which returns
  ## --------- a plotly lyout p_ly    addpc12[[1]]
  ## --------- xhat (?): addpc12[[3]]
  ## --------- counters: addpc12[[4]]
  ## --------- angles:   addpc12[[5]]
  ## ---- x$bipl<-biplot_details[[1]]
 
  ## The following seems a reconstruction from 2 PCs
  ## --- based on function "add_vector_biplot" goes back to Z, p, n, mu, stddev, m,
  ## --- first output returned is also the first computation since SVD... 
  Xhat <- Z%*%t(V)                                       ## Z = UD so time V^T goes back to Y the matrix being decomposed ?
  all.equal(Xhat, approx_Y)                              ## (actually just an approximation)
  Xhat <- sweep(Xhat, MARGIN = 2,STATS=stddev,FUN="*")   ## then multiply by the precalculated "std.dev" 
  Xhat <- sweep(Xhat, MARGIN=2,STATS=mu,FUN="+")         ## then add the mean 
  
  ## The Observation coordinates are just the score Z 
  ## -- the code in  add_Vector_biplot gives it away https://rdrr.io/cran/bipl5/src/R/PCAbiplot_Helper.R
  ## -- the below it's just my take. In the original code the output is returned in such a way to added a plotly trace
  obs_coordinates <- Z
  
  ## The features coordiantes are just the loadings. To be visualised as "arrows"
  ## not easy at allo to spot but the function "add_annotations" gives it away (works with "PC12" but that's just the rank-2 obtained earlier)
  ax = V[,1]
  ay = V[,2]
  x = rep(0,p)                   # I guess this is the origin of each
  y = rep(0,p)
  
  
  ## Part 3 does the layout - for information only        ####
  ## -- this is just for information, I couldn't quite figure what this accomplishes
  ## -- ???? angles from slopes???
  m <- V[, 2]/V[, 1]                                                                              # calls this "vector of slopes" 
  for(i in 1:p){
    angles[[i]]<-list(x=-10*sin(atan(m[i])),y=10*cos(atan(m[i])))
  }
  ## -- ???? Generate a vector of quadrants (reproduced from https://rdrr.io/cran/bipl5/src/R/classbipl5.R)
  getquad<-function(V,m){
    quads<-numeric(length(m))
    for(i in 1:length(m)){
      if(m[i]>0 && V[i,1]>0)
        quads[i]<-1
      if (m[i]>0 && V[i,1]<0)
        quads[i]<-3
      if(m[i]<0 && V[i,1]<0)
        quads[i]<-2
      if(m[i]<0 && V[i,1]>0)
        quads[i]<-4
    }
    return(quads)
  }
  
  ## Part 4 re-does the PCA picking first and third PCs - again not sure ####
  ## -- just a note that "PCAbiplot"  is called within the function "make_biplot", but this time with basis = c(1,3) and calls the result "pc13" and "pc23".
  

  
  
  
  
  
  
}









## ---- 4 aux functions                            ####
isEqualUnname <- function(a, b){
  a1 <- zapsmall(unname(a))
  b1 <- zapsmall(unname(b)) 
  isTRUE(all.equal(a1, b1))
}

isEqualAbs <- function(a, b){
  a1 <- abs(zapsmall(unname(a)))
  b1 <- abs(zapsmall(unname(b)))
  isTRUE(all.equal(a1, b1))
}

## Load data                                       ####
## ---- 0) set relative paths                      ####
mypath <- "./"
Out_path <- paste0(mypath, "Analysis_OUT")                                                             ## save output generated by the script in a separate folder
## ---- 1) Toy example - load data                 ####

## this example is based on Starmer J. https://youtu.be/FgakZw6K1QQ?si=JFFlcEKSg5lNiJ7G
mydata <- read.csv("data/PCA_toy_example.csv", header=TRUE, row.names=1)
X <- (mydata[,1:2])                                                                           # this is the "raw" data in Tab. 2 in the paper
#X <- mydata                                                                                  # extend to 4 features  

## ---- 2) Real(-ish) example - load data          ####
mydata_real <- read.csv("data/X_unitized.csv", header=TRUE, row.names=1)
X_real <- mydata_real




##                                                 ####
## PART I Building blocks + simple example         ####
##                                                 ####
## Section 3.2                                     ####

## ---- 1) Main computations                       -----------------------
## replicates results in sections 2 incl. code listings, and Tab.2
# rownames(X) <- NULL
# colnames(X) <- NULL
# rownames(Y) <- NULL
# colnames(Y) <- NULL

Y <- apply(X,2,function(x) (x - mean(x)))                                                     # matrix of centered features as shown in Tab. 2 in the paper
S <- 1/(nrow(Y)-1) * t(Y) %*% Y                                                               ## covariance matrix of centred data
lambda <- eigen(S)$values                                                                     ## eigenvalues (should be identical to variance of scores)    
V_eig <- eigen(S)$vectors                                                                     ## "loadings" per component. Will get a differnt version later via SVD
Z <- Y %*% V_eig                                                                              ## Scores

## check that scores' variance = eigenvalues of cov matrix, centred data
all.equal(apply(Z,2,var), lambda)                                                             ##  Eq.5 in paper

## variance "explained" by each PC (reported along axis in PCA plot)
var_retain_scores1 <- var(Z[,1])/sum(apply(Y,2,var))
var_retain_scores2 <- var(Z[,2])/sum(apply(Y,2,var))

## eigenvalues and eigenvectors as displayed in code listing
round(lambda,2)
round(V_eig, 2)

## equivalent pre-built options (comments in code listings) 
Y_prebuilt <- scale(X, center = T, scale = F)[,]                                                                     
S_prebuilt <- cov(Y)
all.equal(S_prebuilt, S)  
all.equal(Y_prebuilt, Y)

## ---- 2) Tab 2 - reconstruct proj. coordinates   ####
## reconstruction of orthogonal projections on 2 PCs
orthProj1 <- outer(Z[,1], V_eig[,1]) 
orthProj2 <- outer(Z[,2], V_eig[,2])
colnames(orthProj2) <- colnames(Y)
colnames(orthProj1) <- colnames(Y)

## variance retained by reconstructing data using only the first PC
var_retain_proj <- sum(apply(orthProj1,2,var))/sum(apply(Y,2,var))
all.equal(var_retain_proj, var_retain_scores1)                                   # indeed you can see from Tab.2 that all.equal(var(Z[,1]), sum(apply(orthProj1,2,var)))

## Perfect reconstruction (for information)
Y_recontrcted_perfect <- Z %*% t(V_eig)
all.equal(unname(Y_recontrcted_perfect), unname(Y))

## Generic less than perfect reconstruction with 2 or more PCs
#chose_PCs <- 1
chose_PCs <- 2
#chose_PCs <- 3
Y_recontrcted <- Z[,1:chose_PCs] %*% t(V_eig[,1:chose_PCs])

## compute distance between each point and its projection (complete Tab.2)
## --- PCA "scores" are the distance between the origin and the projectoon
Y_proj <- cbind(Y[,1:2], orthProj1[,1:2])                                        ## in a larger example, there may be more than 2 columns
dist_from_proj1 <- apply(Y_proj, 1, function(x){
  sqrt((x[1] - x[3])^2 + (x[2] - x[4])^2) 
})
Y_proj_2 <- cbind(Y[,1:2], orthProj2[,1:2])
dist_from_proj2 <- apply(Y_proj_2, 1, function(x){
  sqrt((x[1] - x[3])^2 + (x[2] - x[4])^2) 
})




## ---- 2) Tab 2 - assemble Table                  ####
Tab_2 <- cbind(
  round(X, 2),         
  round(Y, 2),                                                                                  ## Tab 2 paper cols 3-4
  round(orthProj1[,1:2], 2),                                                                    ## Tab 2 paper cols 5 & 6
  round(cbind(Z[,1], dist_from_proj1),2),                                                       ## Tab 2 paper cols 7 & 8
  round(orthProj2[,1:2], 4),                                                                    ## Tab 2 paper cols 9 & 10
  round(cbind(Z[,2], dist_from_proj2),4)                                                        ## Tab 2 paper cols 11 & 12
)
colnames(Tab_2) <- c(paste0("feat",1:ncol(Y)," raw"),
                     paste0("feat",1:ncol(Y)," centr"),
                     paste0("feat",1:2," PC1 proj"),
                     "PC1 proj dist orig",
                     "PC1 proj dist observ",
                     paste0("feat",1:2," PC2 proj"),
                     "PC2 proj dist orig",
                     "PC2 proj dist observ"
                    )

## compute mean and sample variance for each column of Tab.2
Tab2_col_means <-  apply(Tab_2,2, mean)                                                          # mean of raw data
Tab2_col_s.var  <- apply(Tab_2,2,var)                                                            # sample variance of raw data, see section 2.1 in the paper
Tab2_col_ms <- rbind(Tab2_col_means, Tab2_col_s.var)

## Display Table 2
round(Tab_2,2)
round(Tab2_col_ms,2)


## ---- 3) Fig 1 - Reproduce                       ####

## THIS WON'T LOOK GOOD IF THE ORIGINAL DATA HAS >2 FEATURES

## optimal slopes of principal components (lines on Cartesian plane)
slope_opt <- unique(round(as.numeric(apply(orthProj1,1,function(x){x[2]/x[1]})),7))           ## Fig 1 LHS in paper
slope_opt2 <- unique(round(as.numeric(apply(orthProj2,1,function(x){x[2]/x[1]})),7))          ## Fig 1 LHS in paper




## FIGURE 1
graphics::layout(mat = matrix(c(2, 1, 0, 3), 
                    nrow = 2, 
                    ncol = 2),
       heights = c(0.5, 2),    # Heights of the two rows
       widths = c(2, 0.5))     # Widths of the two columns
par(mar = c(4, 4.5, 0, 0))
plot(Y[,1], Y[,2],
     cex = 2,
     xlim = c(-6,6),
     ylim = c(-6,6),
     xlab = colnames(Y)[1],
     ylab = colnames(Y)[2],
     cex.lab=1.5, cex.axis=1.5
)
abline(h = 0, v = 0, col = "gray60")
curve(expr = x*slope_opt, add = TRUE, col = "red")                                                                                              # thread: 
curve(expr = x*slope_opt2, add = TRUE, col = "blue")
arrows(Y[,1],Y[,2], orthProj1[,1], orthProj1[,2], length = .05, angle =8, col = "purple")
arrows(Y[,1],Y[,2], orthProj2[,1], orthProj2[,2], lty="dashed", length = .05, angle =8, col = "purple")
text(Y[,1], Y[,2], labels=rownames(Y), cex= 2, pos=3)


## FIG 2 A
## Multiplication by -1 helps keep the visual analogy with original plot but alters scores. Here I make it optional
#pca_plot_coorz <- -1*Z
pca_plot_coorz <- Z

plot.new()
layout(mat = matrix(c(2, 1, 0, 3),
                    nrow = 2,
                    ncol = 2),
       heights = c(0.5, 2),    # Heights of the two rows
       widths = c(2, 0.5))     # Widths of the two columns

# Plot 1: Scatterplot
par(mar = c(4.5, 4.5, 0, 0))
plot(pca_plot_coorz,
     pch = 18,             # point shape
     cex = 2,              # point size
     xlim = c(-6,6),
     ylim = c(-1.8,1.8),
     xlab = paste0("PC 1 (var expl: ", round(var(Z[,1])/sum(apply(Y,2,var))*100,2),"%)"),
     ylab = paste0("PC 2 (var expl: ", round(var(Z[,2])/sum(apply(Y,2,var))*100,2),"%)"),
     cex.lab=1.5, cex.axis=1.5
)
abline(h = 0, v = 0, col = "gray60")
text(pca_plot_coorz[,1], pca_plot_coorz[,2], labels=rownames(Y), cex= 2, pos=3)

## FIG2 B
# Plot 2: Top
par(mar = c(0, 4.5, 0, 0))
x <- data.frame(pca_plot_coorz[,1],1)                      # notice multiplication by -1 to help overlaying the two chard
plot(x,
     cex = 1.5,
     xlim=c(-6,6),
     axes = F,
     ylab = "",
     xlab = "projections on PC1",
     xaxt="n",                                   #remove tick marks
     type = "o",
     #bty = "n",                                 # remove bounding box
     col = "red",
     #pch = 19                                    # shape
)
#axis(side=1)
text(x, labels=rownames(pca_plot_coorz), cex= 1.5, pos=3)

# Plot 3: Right
par(mar = c(4, 2, 0, 0))
#par(fig = c(0.65,1,0,0.8, new=TRUE))
x2 <- data.frame(1,pca_plot_coorz[,2])                      # notice multiplication by -1 to help overlaying the two chard
#par(mar = c(8,3,8,3))
plot(x2,
     cex = 1.5,
     ylim=c(-1.8,1.8),
     axes = F,
     xlab = "",
     ylab = "projections on PC2",
     xaxt="n",                                   #remove tick marks
     type = "o",
     #bty = "n",                                 # remove bounding box
     col = "blue",
     #pch = 15                                    # shape
)
#axis(side=4)
text(x2, labels=rownames(pca_plot_coorz), cex= 1.5, pos=2)






## Section 3.3                         ####
## ---- 4) SVD                         ---------------------------------------------------------------
svd_Y <- svd(Y)                                                                               # singular value decomposition of a rectangular matrix of centred features
U <- svd_Y$u                                                                                  # left eigenvectors
V <- svd_Y$v                                                                                  # right eigenvectors equivalent to loading V_eig but different signs
ell <- svd_Y$d                                                                                # singular values (NOT eigenvalues)
D <- diag(as.vector(ell))                                                                     
Z_svd <- unname(U %*% D)                                                                      # Paper eq. 10: scores matrix, SVD style

## check: singular values vs eigenvalues vs scores St.Dev.
all.equal(ell, sqrt(nrow(Y)-1)*apply(Z,2,sd))                                                 # paper Eq. 7
all.equal(sqrt((nrow(Y)-1)*eigen(S)$values), sqrt(nrow(Y)-1)*apply(Z,2,sd))                   # combining  Eq.5 & 7 in paper

## Display (code shown in paper, Section 2.3)
round(U,2)                                                                                    
round(V,2)                                                                                    
round(D,2)
round(Z_svd, 2)



## ---- 5) Biplot SVD                  ---------------------------------------------------------------
A <- U                                                                                        # Eq. 8 & 11 in paper: observations coordinates
all.equal(A, sweep(Z_svd,2,ell, "/"))                                                         # equivalent, but redundant. You find it in biplot script because it takes PCA scores as input
B_trsp <- D %*% t(V)                                                                          # Eq. 12 in paper:  features coordinates
colnames(B_trsp) <- colnames(Y)
rownames(B_trsp) <- paste0("dim_",seq(1:nrow(B_trsp)),sep = "")
B <- t(B_trsp)

## check equialence in eq. 11 (holds in absolute values when using Z instead of Z_svd)
scaled_scores <- sweep(Z,2,svd_Y$d, "/")
all.equal(abs(unname(scaled_scores)), abs(A))

## Display feat = COLUMNS OF B' (= ROWS OF B), but just 2 PCs (= ROWS of B')
round(B_trsp[1:2,], 2) 



## ---- 6) biplot simple example       ####
## The next bit is to export a higher figure resolution, goes before the actual plot 
## thread: https://stackoverflow.com/questions/22813723/how-can-i-increase-the-resolution-of-my-plot-in-r
## https://youtu.be/1SyLtXskq2g?si=LxLg6_B8ZZlS7JhJ
#png(paste0(Out_path,"/biplot1_higherres.png",collapse = " "), res=300, 
#    width=1300, 
#    height=700)
 

## actual biplot                                      
biplot_feat_coord <- t(B_trsp)                                      # if we use B: cols 1 and 2 = x, y coords, else rows
#biplot_observ_coord <- A
biplot_observ_coord <- scaled_scores                                # comparable with Fig 1

## composite figure - 3-part layout
layout(mat = matrix(c(2, 1, 0, 3), 
                    nrow = 2, 
                    ncol = 2),
       heights = c(0.5, 2),    # Heights of the two rows
       widths = c(2, 0.5))     # Widths of the two columns

## Plot 1: Scatterplot
## -- fix margins so to show the axis' labels
par(mar = c(4, 4, 0, 0))
## -- Scale down features' coordinates as they may be way larger than the observations'
adjust_viz <- abs(min(biplot_observ_coord)/min(biplot_feat_coord))
biplot_feat_coord_viz <- biplot_feat_coord*adjust_viz
## -- Shorten features' labels for display
## -- thread https://stackoverflow.com/questions/11776287/remove-pattern-from-string-with-gsub
attr_short_label <- gsub(".*_","", colnames(Y))


## -- Set axis' range considering both features and observations
lim_max1ob <- max(biplot_observ_coord[,1])*1.3
lim_min1ob <- min(biplot_observ_coord[,1])*1.3
lim_max2ob <- max(biplot_observ_coord[,2])*1.3
lim_min2ob <- min(biplot_observ_coord[,2])*1.3
lim_max1feat <- max(biplot_feat_coord_viz[,1])*1.3
lim_min1feat <- min(biplot_feat_coord_viz[,1])*1.3
lim_max2feat <- max(biplot_feat_coord_viz[,2])*1.3
lim_min2feat <- min(biplot_feat_coord_viz[,2])*1.3
## -- plot the actual biplot
par(mar = c(4, 4.2, 0, 0))     
## -- observations' plot
plot(biplot_observ_coord[,1:2],
     pch = 18,                # point shape
     cex = 2,               # point size
     xlim = c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat)),
     ylim = c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat)),
     xlab = paste0("PC 1 (var expl: ", round(lambda[1]/sum(lambda),2)*100,"%)"),
     ylab = paste0("PC 2 (var expl: ", round(lambda[2]/sum(lambda),2)*100,"%)")
)
abline(h = 0, v = 0, col = "gray60")
text(biplot_observ_coord[,1], biplot_observ_coord[,2], cex=1.5, pos = 3, labels = rownames(Y))
## Overlay "loading plot"
arrows(rep(0,2),rep(0,2), biplot_feat_coord_viz[,1], biplot_feat_coord_viz[,2], length = .1, angle =10, col = "purple")
text(biplot_feat_coord_viz[,1], biplot_feat_coord_viz[,2], pos = 4, cex=1.5, labels = attr_short_label)

## instead of a dual axis I will plot "around" the main plot
## -- features' axes range
lim_max1 <- max(biplot_feat_coord[,1])
lim_min1 <- min(biplot_feat_coord[,1])
lim_max2 <- max(biplot_feat_coord[,2])
lim_min2 <- min(biplot_feat_coord[,2])
x_range_feat <- c(lim_min1, lim_max1)
y_range_feat <- c(lim_min2, lim_max2)

## Plot 2: Top - this will represent the horizontal scale for the feature's coordinates 
par(mar = c(0, 4.2, 2.2, 0))
x <- data.frame(biplot_feat_coord[,1],1)                      # notice multiplication by -1 to help overlaying the two chard
plot(x, 
     xlim = c(min(lim_min1ob, lim_min1feat), max(lim_max1ob,lim_max1feat))/adjust_viz,
     axes = F,
     ylab = "",
     xaxt="n",                                   #remove tick marks
     type = "o",
     col = "purple",
     #pch = 15                                    # shape
)
text(x, pos=2, labels=attr_short_label, cex= 1.5)
axis(side=3, cex.axis = 1, col = "grey", las=1)

## Plot 3: Right - feature's coordinates
par(mar = c(4, 0, 0, 2.2))
x2 <- data.frame(1,biplot_feat_coord[,2])                      # notice multiplication by -1 to help overlaying the two chard
plot(x2, 
     ylim = c(min(lim_min2ob, lim_min2feat), max(lim_max2ob, lim_max2feat))/adjust_viz,
     axes = F,
     xlab = "",
     xaxt="n",                                   #remove tick marks
     type = "o",
     col = "purple",
     #pch = 15                                    # shape
)
text(x2, pos=3, labels=attr_short_label, cex= 1.5) 
axis(side=4, cex.axis = 1, col = "grey", las=1)


## close device (when saving higher res)
#dev.off()

## Section 3.4                         ####
## ---- 7) feat vec length - eq14      ####
## Check eq. 14: length of features vectors in biplot vs feats St.Dev
eq14 <- apply(B_trsp,2,function(x){norm(x, "2")/sqrt(nrow(Y)-1)})       # this is the second row in eq 14
all.equal( unname(sqrt(diag(S))), unname(eq14))  # Eq. 14 in paper


## now check the next lines in eq 14 as it fully develops
all.equal(unname(B), V%*%D)                     # checking that B = VD because that's  v_jk*l_k in the third line
alt_eq14_a <- sapply(1:ncol(Y),function(k){
  #B_trsp[,k]^2/(nrow(Y)-1)                          # this is the output of the third row in eq.14  
  (V_eig[,k]*sqrt(lambda[k]))^2                 # equivalent, this is the last row in eq.14
  })
alt_eq14_b <- sqrt(apply(alt_eq14_a,1,sum))
all.equal(alt_eq14_b, unname(sqrt(diag(S))))    # should work



## ---- 8) feat vec length - alt w eigen         #### 
## alternative features coordinates using eigenvalues not singular values
l_diag <- diag(as.vector(lambda))
B_alt <- V_eig %*% sqrt(l_diag)
all.equal( unname(sqrt(diag(S))), apply(B_alt,1,function(x){norm(x, "2")}))

## ---- 9) feat vec length - rank2 approx fail   ####
## rank-2 approximation (2-dimenional biplot) in the presence of more than 2 features
U_r2 <- U[,1:2]
D_r2 <- D[1:2,1:2]
V_r2 <- V[,1:2]
A_r2 <- U_r2                                                                                  # biplot observation coords in 2 dimension, 
B_trsp_r2 <- D_r2 %*% t(V_r2)                                                                 # biplot feature coords in 2 dimensions
Y_r2 <- A_r2 %*% B_trsp_r2                                                                    # Eq. 9 in paper

## verify again length of features vector vs feat. St.Dev (no identity)
length_feature_vec_r2 <- apply(B_trsp_r2,2,function(b){                                             
  ## for each column of B'
  norm(b, "2")
}) 
modif_length_feat_vec_r2 <- unname(length_feature_vec_r2*(1/sqrt(nrow(Y)-1)))

## display how Eq.14 fails when we use a Rank-2 approx
all.equal(apply(Y,2,sd), modif_length_feat_vec_r2)


## ---- 9B) length B only vs rank-2              ####
length_b <- apply(B_trsp,2,function(x){
  norm(x, "2")
})
length_b_r2 <- apply(B_trsp[1:2,],2,function(x){
  norm(x, "2")
})

## check
all.equal(unname(length_b)/sqrt(nrow(Y)-1), unname(sqrt(diag(S))))

## for display
length_B <- cbind(length_b_r2, length_b, sqrt(diag(S)))
rownames(length_B) <- colnames(Y)
colnames(length_B) <- c("rank2", "allPCs", "sigma_y")

## ----10) Correlation vs cosine                 -----------------------------------------------------
## A bit overcooked for a small example, but comes in handy when we have more than 2 features

## Prepare correlations grid
idx_grid <- expand.grid(1:nrow(B_trsp), 1:nrow(B_trsp))
idx_grid <- idx_grid[,c(2:1)]
idx_grid <- idx_grid[which(idx_grid[,2] > idx_grid[,1]),]                                     # just keep the upper triangular part...

## cosine of angle between features vectors                                                   # Eq. 15 in paper
## -- Rank 2 approximation (diplayed in a biplot)
cos_theta_ab <- apply(idx_grid,1,function(x){
  i <- x[1]
  j <- x[2]
  ## only first two PCs
  load_a <- B_trsp[1:2,i]
  load_b <- B_trsp[1:2,j]
  (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
})

## -- no dimension reduction: the equivalence between cos theta and corr requires this
cos_theta_ab_allPCs <- apply(idx_grid,1,function(x){
  i <- x[1]
  j <- x[2]        #
  ## ALL pcs to compare to corr
  load_a <- B_trsp[,i]
  load_b <- B_trsp[,j]
  (load_a %*% load_b) / (norm(load_a, "2")*norm(load_b, "2"))                                 
})



theta_ab <- acos(cos_theta_ab)                                                                #  angle from 0 to pi radians 
theta_ab_deg <- theta_ab*180/pi                                                               # "degrees" instead of radians

## Pearsons product-moment correlation                                                        # Eq. 16 in paper but comapct version
correl_feat <- apply(idx_grid,1,function(x){
  i <- x[1]
  j <- x[2]
  cor_ab <- cor(Y[,i], Y[,j])
  test_cor_ab <- cor.test(Y[,i], Y[,j])
  c(cor_ab, test_cor_ab$p.value)
})
correl_feat <- as.data.frame(t(correl_feat))
colnames(correl_feat) <- c("corr_coeff", "p-value")                                           #  this tests the H_0: true correlation IS equal to 0 [here is no linear relationship between the two variables]; alternative hypothesis: true correlation is not equal to 0

## Alternative, equivalent 
correl_feat_fromScratch <- apply(idx_grid,1,function(x){                                      # this is actually Eq. 16 in paper
  i <- x[1]
  j <- x[2]       
  ## ALL pcs to compare to corr
  y_a <- Y[,i]
  y_b <- Y[,j]
  (y_a %*% y_b) / (norm(y_a, "2")*norm(y_b, "2"))                                             # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
})
all.equal(unname(correl_feat_fromScratch), correl_feat$corr_coeff)

## Check n.4: cos angle vs correlation
all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab),6))          # won't work for rank 2, it's just an approx
all.equal(round(as.numeric(correl_feat$corr_coeff),6), round(as.numeric(cos_theta_ab_allPCs),6))   # works


## table to export
coscorr_export_tab <- cbind(idx_grid, theta_ab_deg, cos_theta_ab, cos_theta_ab_allPCs, correl_feat)





## ----11) Correlation with PC                   ####

## based on Johnson and Wichern Ch8 Result 8.3

## this is convenient to link up to eq14 but LOSES THE SIGN
#eq19_rhs <- sweep(sqrt(alt_eq14_a),1,sqrt(diag(S)),"/")                                  ## right hand side of Eq. 19 in paper, for every j (rows) and k (Cols)

## this version avoids squaring and then taking the sqrt to keep the SIGNS
alt_eq14_sqrtArg <- sapply(1:ncol(Y),function(k){
  (V_eig[,k]*sqrt(lambda[k]))                 # equivalent, this is the last row in eq.14
})
eq19_rhs <- sweep(alt_eq14_sqrtArg,1,sqrt(diag(S)),"/")

## tiny example
k <- 1
j <- 2
z_k <- Z[,k]
y_j <- Y[,j]
corr_feat_PC <- cor(z_k, y_j) 
all.equal(eq19_rhs[j,k], corr_feat_PC)

## Now extend to ALL corr between PC and features
## a) Prepare correlations grid
idx_grid3 <- expand.grid(1:ncol(Y), 1:ncol(Z))
idx_grid3 <- idx_grid3[,c(2:1)]
#idx_grid3 <- idx_grid3[which(idx_grid3[,2] > idx_grid3[,1]),]                        # just keep the upper triangular part...

## b) Pearsons product-moment correlation                                             # Eq. 16 in paper but comapct version
correl_feat_PC <- apply(idx_grid3,1,function(x){                                      # this is actually Eq. 16 in paper
  k <- x[1]
  j <- x[2]       
  ## ALL pcs to compare to corr
  z_k <- Z[,k]
  y_j <- Y[,j]
  (z_k %*% y_j) / (norm(z_k, "2")*norm(y_j, "2"))                                             # equivalent to the correlation coefficient between "centred" variables (see my Poiwerpoint notes Jan 2024)
})

all.equal(correl_feat_PC, c(eq19_rhs))


## ----12) obs vec: Mahalanobis distance         ####

## Jolliffe P.77-78; Gabriel 1971 also mentioned by Venables and Ripley

## example for an arbitrary rows pair from the centred data matrix
x <- as.matrix(Y[1,] - Y[5,])
mymahala <- as.numeric(t(x) %*% solve(S, x))
rmahala <- mahalanobis(Y[1,],center = Y[5,], cov = S)  ## built-in R function
all.equal(mymahala, rmahala)

## small for use in the paper due to space
x <- as.matrix(Y[1,] - Y[2,])
x2 <- as.matrix(A[1,] - A[2,])
mymahala <- as.numeric(t(x) %*% solve(S, x))
adj_dist2 <-as.numeric((nrow(Y)-1)*t(x2)%*%x2) 
all.equal(adj_dist2, mymahala)



## scale up...all observation pairs
idx_grid2 <- expand.grid(1:nrow(Y), 1:nrow(Y))
idx_grid2 <- idx_grid2[,c(2:1)]
idx_grid2 <- idx_grid2[which(idx_grid2[,2] > idx_grid2[,1]),]                                     # just keep the upper triangular part...

mymahala_all <- apply(idx_grid2, 1, function(k){
  i <- k[1]
  j <- k[2]
  x <- as.matrix(Y[i,] - Y[j,])
  as.numeric(t(x) %*% solve(S, x))
})

## now the squared euclidean distance between biplot coordinates
sqdist_biplot_coord <- apply(idx_grid2, 1, function(k){
  i <- k[1]
  j <- k[2]
  x <- as.matrix(A[i,] - A[j,])
  t(x)%*%x
})


## check the second is "proportional" to the first
all.equal((nrow(Y)-1)*sqdist_biplot_coord, mymahala_all)



##                                               ####
## PART II Compare + realistic example           ####
##                                               ####
## ---- Z   back-to-basic framework              ####
## ---- Z.0 Get framework results (for comparison)  ----------
ext_PCA <- YAPCAB(X_real, myscale = "none")                      # shared data has been desensities using "unit" upfront and requires no further processing           

## run properties' tests from my grid on my function (unnecessary but shows the process)
my_Y <- ext_PCA$data_modified
my_S <- ext_PCA$covar_mat
my_D <- ext_PCA$PCA_svd$D                                        # diagonalised singular values from SVD
my_sv <-  apply(my_D,2,sum)
my_U <- ext_PCA$PCA_svd$U
my_V <- ext_PCA$PCA_svd$V
my_V_eig <- ext_PCA$PCA_eigen$loadings
my_Z_svd <- ext_PCA$PCA_svd$Z_svd
my_Z <- ext_PCA$PCA_eigen$scores
my_lambda_svd <- my_sv^2/(nrow(my_Y)-1)                          # based on Eq. 7 in my paper
my_lambda <- ext_PCA$PCA_eigen$eigenvalues                       
my_biplot_feat_coord <- ext_PCA$biplot$B_trsp                    # B transposed
my_biplot_observ_coord <- ext_PCA$biplot$A                       # A

## alternative coordinates to compare for instance with princomp
my_biplot_observ_coord_eig <- t(t(my_Z)/sqrt(my_lambda*(nrow(my_Y)-1)))
my_biplot_feat_coord_eig <- my_D %*% t(my_V_eig)

## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only) -- V.A
ext_lentest <- lenTest(my_biplot_feat_coord, my_Y, my_S)
ext_lentest$test_len

## (E.2) test 2: correlation between features -- VI.A
ext_corcostest <- corcosTest(my_biplot_feat_coord , my_Y)
ext_corcostest$coscorr_test
ext_corcostest$coscorr_export_tab

## (E.3) test 3: correlation between scores/PCs and features -- VI.B 
ext_corPC <- corPCTest(V = my_V, Z = my_Z_svd, Y = my_Y, lambda = my_lambda_svd, S = my_S)
ext_corPC$check_eq19
ext_corPC$corr_feat_PC


## (E.4)  dist between observ on PC biplot = Mahalanobis distance between observ in data matrix -- VII
ext_mahala <- MahalaTest(A = my_biplot_observ_coord, Y = my_Y,  S = my_S)
ext_mahala$mahalaTest









## ---- Z.1 Draw biplot                          ####
showmetheBplot(YAPCAB_out = ext_PCA)



## ---- Z.3 starting data matrices for packages  ####
## Assuming here I have "unitised" to anonymise the data
X <- unitise_cols(X_real)
Y <- scale(X, scale = FALSE)[,]

##                                               ####
## ---- A benchmark vs base-R          ####
## ---- A.1 base R prcomp + biplot     ---------------------
## View(stats:::prcomp.default)
## View(stats:::biplot.prcomp)

## (A) call pre-built function 
test_pca1 <- prcomp(X, center = TRUE, scale = FALSE)
Z_pca1 <- test_pca1$x                                                                         # in the script this is Y %*% V
V_pca1 <- test_pca1$rotation                                                                  # loadings (eigenvectors)
stDev_pca1 <- test_pca1$sdev                                                                  # instead of singular values (sv) this returns  sv/sqrt(n - 1) see Eq.7 in paper
sv_pca1 <- test_pca1$sdev*sqrt(nrow(Y)-1)                                                     # available internally but nor returned, hence had to be reproduced 


## (B) compare with from scratch SVD PCA (sv not returned)
isEqualUnname(my_Z_svd, Z_pca1)                                                               # scores as per Eq. 10 in paper
isEqualUnname(V_pca1, my_V)                                                                   # loadings = right eigenvectors (SVD)
isEqualUnname(sv_pca1, my_sv)
isEqualUnname(stDev_pca1^2, my_lambda)

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
isEqualUnname(stDev_pca1^2, my_lambda_svd)                                                    # verifies Eq.7 i.e., sv/sqrt(n - 1) = sqrt(eigenvalues of S) 
isEqualUnname(stDev_pca1, apply(test_pca1$x,2,sd))                                            # verifies Eq.5 i.e., sv/sqrt(n - 1) = St.Dev of scores 


## (D.1) Biplot coordinates: the function biplot does not return them
test_biplot1 <- biplot(test_pca1)                                                             # does the plot without returning the output of the computations


## (D.2) Replication of biplot coordinates (call own function)
## -- First reproduce coordinates for a PC biplot (to show discrepancy)
biplot_test_pca1_orig <- replica_biplot(test_pca1, original = T) 
isEqualUnname(biplot_test_pca1_orig$obser_coord, my_biplot_observ_coord)                                               
isEqualUnname(biplot_test_pca1_orig$feat_coord, t(my_biplot_feat_coord))     

## -- Re-run replica, this time applying "fixes" to amend the source of discrepancies (see function annotations for details)
biplot_test_pca1_fix <- replica_biplot(test_pca1, original = F) 
isEqualUnname(biplot_test_pca1_fix$obser_coord, my_biplot_observ_coord)                                               
isEqualUnname(biplot_test_pca1_fix$feat_coord, t(my_biplot_feat_coord))   

## -- Finally, show that the coordinates for a scores and loading biplot are generated correctly (alpha = 1 in the paper's notation)
biplot_test_pca1_orig_sl <- replica_biplot(test_pca1, original = T, myscale = 0) 
isEqualUnname(biplot_test_pca1_orig_sl$obser_coord, my_Z_svd)                                             
isEqualUnname(biplot_test_pca1_orig_sl$feat_coord, my_V)   



## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only)
## -- caveat 1: the test requires all the dimensions (in all cases, the function does not produce a covariance matrix for Y)
## -- caveat 2: SVD functions do not compute a covariance matrix for Y
pca1_test_len_orig <- lenTest(t(biplot_test_pca1_orig$feat_coord), Y, S = NULL)             # original
pca1_test_len_orig$test_len
pca1_test_len_orig$tab_out

## notice that it would have worked if we had used sqrt(n) instead of sqrt(n-1)
all.equal(pca1_test_len_orig$tab_out[,1]/sqrt(nrow(Y)), pca1_test_len_orig$tab_out[,3])


## (E.2) test 2: correlation between features
pca1_test_corcostest_orig <- corcosTest(myB_trsp = t(biplot_test_pca1_orig$feat_coord), myY = Y)
pca1_test_corcostest_orig$coscorr_test
pca1_test_corcostest_orig$coscorr_export_tab


## (E.3) test 3: correlation between scores and features
pca1_test_featPC <- corPCTest(V = V_pca1, Z = Z_pca1, Y = Y, lambda = stDev_pca1^2)
pca1_test_featPC$check_eq19
pca1_test_featPC$corr_feat_PC

## (E.4) test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- caveat: in principle I couldn't run this test because neither prcomp nor biplot give me the variance-covariance matrix of Y
pca1_test_mahala_orig <- MahalaTest(A = biplot_test_pca1_orig$obser_coord, Y = Y, S = cov(Y))
pca1_test_mahala_orig$mahalaTest
pca1_test_mahala_orig$dist_tab

## notice: it would have worked if we had used n instad of n-1
all.equal(pca1_test_mahala_orig$dist_tab[,3], pca1_test_mahala_orig$dist_tab[,5]*nrow(Y))






## ---- A.2 base R princomp + biplot         ####

## -- THIS ONE USES EIGEN() not SVD(). From function help:
## ----- [...] note that the default calculation uses divisor N for the covariance matrix"
## -- to see the source code: stats:::princomp.default 
test_pca10 <-  princomp(X, fix_sign = FALSE)        # default is TRUE but might not match my eigendecomposition
sd_pca10 <- test_pca10$sdev                         # sqrt of eigenvalues of (1/n)Y'Y NOT S
lambda_pca10 <- sd_pca10^2                          # eigenvalues computed but not returned, only the square root
V_pca10 <- test_pca10$loadings[,]                   # loadings (eigenvectors)
Z_pca10 <- test_pca10$scores                        # PCA scores
all.equal(Z_pca10, Y %*% V_pca10)                   # as expected

## compare with "from scratch" eigendecomposition PCA
isEqualUnname(my_Z, Z_pca10)                        # signs may differ
isEqualUnname(my_V_eig, V_pca10)                    # signs may differ
isEqualUnname(my_lambda, lambda_pca10)              # eigenvalues differ due to different covariance matrix 

## scores and loadings ok in absolute value
isEqualAbs(my_Z, Z_pca10)
isEqualAbs(my_V_eig, V_pca10)

## Test Eq.5 and 7
all.equal(apply(Z_pca10,2,var),lambda_pca10 )             # Eq 5 doesn't work
isEqualUnname(sd_pca10, my_sv/sqrt(nrow(Y)-1))            # Eq.7 doesn't work

## Now I run my replica to extract the covariance matrix and to enable "fixes" compared to my from-scratch computations
## -- check replica is correct
test_pca10_replica1 <- replica_princomp(X, original = T, test_symmetric = T, fix_sign = F)             # meant to reproduce the original results
isEqualUnname(test_pca10_replica1$loadings, V_pca10)
isEqualUnname(test_pca10_replica1$scores, Z_pca10)
isEqualUnname(test_pca10_replica1$sdev,sd_pca10)
## -- notice however that eigen is very sensitive to how the (wrong) covariance matrix is arrived at, even if the output is identical
test_pca10_replica1$tests$test_cov_ok2              # same (wrong) covariance matrix, from scratch vs pre-built approach 
test_pca10_replica1$tests$test_ev_ok                # same eigenvalues
test_pca10_replica1$tests$test_evec_ok              # ...and yet eigen produces eigenvectors with different signs. Unless test_symmetric = F the output will differ even if the cov matrix is the same on paper. But then I won't be able to replicate
## -- now retrieve covariance matrix for default case using own replica
S_pca10 <- test_pca10_replica1$cov_mat


## I can address the source discrepancies by running princomp with my fixes
test_pca10_fix <- replica_princomp(X, original = F, fix_sign = F)   # use the "correct" covariance matrix; avoid manipulating the signs
all.equal(unname(my_Z), unname(test_pca10_fix$scores))              # due to issues explained above (n instead of n-1; symmetric in eigen and weird sample cov matrix)
all.equal(unname(my_V_eig), unname(test_pca10_fix$loadings))                             # due to issues explained above (n instead of n-1; symmetric in eigen and weird sample cov matrix)
S_pca10_fix <- test_pca10_fix$cov_mat

## Test Eq.5 and 7 with "fixes
all.equal(apply(unname(Z_pca10),2,var),test_pca10_fix$eigenvalues)       # Now Eq 5 should work
isEqualUnname(sqrt(test_pca10_fix$eigenvalues), my_sv/sqrt(nrow(Y)-1))   # Now Eq.7 should work


## PC Biplot (alpha = 0) ##
## -- base-R Biplot takes princomp output as arguemt
## -- just like before, I run my replica of the original script to get the ourput
biplot_test_pca10 <- replica_biplot(x = test_pca10, original = T,  my_princomp = TRUE)            # PCA biplot but WITHOUT Gabriel's 1971 correction in Eq.52
isEqualUnname(biplot_test_pca10$obser_coord, my_biplot_observ_coord_eig)                          # discrepancy in observations coordinates
isEqualUnname(biplot_test_pca10$feat_coord, t(my_biplot_feat_coord_eig))                          # discrepancy in features coordinates
isEqualAbs(biplot_test_pca10$obser_coord, my_biplot_observ_coord_eig)
isEqualAbs(biplot_test_pca10$feat_coord, t(my_biplot_feat_coord_eig))

## Again, I can address the source discrepancies by using my replica instead, with "fixes" on 
biplot_test_pca10_fix <- replica_biplot(x = test_pca10_fix, original = F,  my_princomp = TRUE)                             
isEqualUnname(biplot_test_pca10_fix$obser_coord, my_biplot_observ_coord_eig)                     # now all checks out                                
isEqualUnname(biplot_test_pca10_fix$feat_coord, t(my_biplot_feat_coord_eig))  

## Scores and loading biplot (alpha = 1): signs issue
biplot_test_pca10_orig_sl <- replica_biplot(test_pca10, original = T, myscale = 0,  my_princomp = TRUE) 
isEqualUnname(biplot_test_pca10_orig_sl$obser_coord,  my_Z)                                      ## -- won't work due to signs                                     
isEqualUnname(biplot_test_pca10_orig_sl$feat_coord, my_V_eig)                                    ## -- won't work due to signs
isEqualAbs(biplot_test_pca10_orig_sl$obser_coord,  my_Z)
isEqualAbs(biplot_test_pca10_orig_sl$feat_coord, my_V_eig)

## Geometrical Properties
## test 1: features' vector length
## -- caveat: this test requires all the dimensions
pca10_test_len_origS <- lenTest(t(biplot_test_pca10$feat_coord), Y, S = S_pca10)     # here I reveal that this test can fail because of how the covariance matrix is computed in princomp
pca10_test_len_origS$testS
pca10_test_len_origS$test_len
pca10_test_len_origS$tab_out                                                         # works using St.Dev. of features but no identity with diag elements S


## test 2: correlation between features
pca10_test_corcostest_orig <- corcosTest(t(biplot_test_pca10$feat_coord), Y)
pca10_test_corcostest_orig$coscorr_test


## Properties test 3: correlation between scores and features
pca10_test_featPC <- corPCTest(V = V_pca10, 
                              Z = Z_pca10, 
                              Y = Y, 
                              lambda = lambda_pca10,  
                              S = S_pca10)
pca10_test_featPC$testS
pca10_test_featPC$check_eq19
pca10_test_featPC$corr_feat_PC



## Test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- caveat: in principle I couldn't run this test because neither prcomp nor biplot give me the variance-covariance matrix of Y
pca10_test_mahala <- MahalaTest(A = biplot_test_pca10$obser_coord, Y = Y, S = S_pca10)
pca10_test_mahala$mahalaTest                                                                       # fails when we consider its actual covariance matrix
pca10_test_mahala$dist_tab

## notice: it would have worked if we had used n instad of n-1
all.equal(pca10_test_mahala$dist_tab[,3], pca10_test_mahala$dist_tab[,5]*nrow(Y))



##                                           ####
## ---- B benchmark vs SVD packages          ####
## ---- B.0 Load Packages (for benchmark)    -----------------------------------------------------
library(ggbiplot)
library(pcaMethods)
library(PCAtools)

## libraries from bioconductor
# https://www.bioconductor.org/install/
# -- https://bioconductor.org/packages/release/bioc/html/pcaMethods.html
# -- https://bioconductor.org/packages/release/bioc/html/PCAtools.html
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcaMethods")
# BiocManager::install("PCAtools")
# BiocManager::install("BiocSingular")

## ---- B.1a [BENCHMARK] - pcaMethods, pca()  ####
## in essence, this is just prcomp() operating under the hood: pca() --> svdPca() -->  prcomp()
## NOTICE: when calling prcomp via svdPca the default is: prcomp(Matrix, center = FALSE, scale. = FALSE)

## (A) call pre-built function 
test_pca2 <- pcaMethods::pca(X, nPcs = ncol(X), method = "svd")
Z_pca2 <- test_pca2@scores                                                                    # scores
V_pca2 <- test_pca2@loadings                                                                  # loadings (eigenvectors)
stDev_pca2 <- test_pca2@sDev                                                                  # scores standard deviation 
sv_pca2 <- test_pca2@sDev*sqrt(nrow(Y)-1)


## (B) compare with from scratch SVD PCA (sv not returned)
## -- here I just prove it's identical to prcomp
all.equal(test_pca2@scores, test_pca1$x)
all.equal(test_pca2@loadings, test_pca1$rotation)
all.equal(as.numeric(test_pca2@sDev), test_pca1$sdev)
all.equal(stDev_pca2^2, my_lambda)

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
isEqualUnname(stDev_pca2, sqrt(my_lambda_svd))                                               ## check: Eq. 7
isEqualUnname(apply(Z_pca2, 2, var), stDev_pca2^2)                                           ## check: Eq. 5



## ---- B.1b [BENCHMARK] - pcaMethods, biplot() and slplot() ####

## (D.1) Biplot coordinate, 2 options: 
## -- D.1.1) replica of stats::biplot -- I will omit this since prcomp and biplot have already been covered
## -- see View(pcaMethods:::biplot.pcaRes)
## -- see also https://rdrr.io/bioc/pcaMethods/man/biplot-methods.html
biplot_test_pca2 <- pcaMethods:::biplot(test_pca2)                                             # returns no output except the graphical device

## -- D.1.2)  slplot: a "Side by side scores and loadings plot"
## -- see e.g. https://rdrr.io/bioc/pcaMethods/man/slplot-pcaRes-method.html#heading-3; and see https://rdrr.io/bioc/pcaMethods/src/R/methods-pcaRes.R 
biplot_test_pca2_sl <- pcaMethods::slplot(test_pca2)                                           # returns some output but it's entirely about the graphical device

## (D.2) Replication of biplot coordinates
## -- skipped because scores and loadings already provided, 
## -- no further computations related to the biplot coordinates seems to take place in slplot except for graphical representation


## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only)
## -- omitted:  not a PC biplot 

## (E.2) test 2: correlation between features
## -- omitted:  not a PC biplot

## (E.3) test 3:correlation between scores and features
pca2_test_featPC <- corPCTest(V = test_pca2@loadings, Z = test_pca2@scores, Y = Y, lambda = stDev_pca2^2)
pca2_test_featPC$check_eq19
pca2_test_featPC$corr_feat_PC

## (E.4) test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- omitted: not a PC biplot


## ---- B.2a [BENCHMARK] - PCAtools, pca()    ####

## Like pcaMethods, this is also reminiscent of princomp but does not call it directly
## -- the relevant functions are borrowed from  BiocSingular https://rdrr.io/bioc/BiocSingular/man/runSVD.html
## -- pca() --> runPCA() --> runSVD() from a different package BiocSingular

## Resources
## -- View(PCAtools::pca)  
## -- View(BiocSingular::runExactSVD)
## -- findMethods("runPCA")
## -- findMethods("runSVD")
## -- https://rdrr.io/bioc/BiocSingular/src/R/runPCA.R
## -- https://rdrr.io/bioc/BiocSingular/src/R/runSVD.R

## (A) call pre-built function 
## -- excerpt from runPCA page: "This function simply calls runSVD and converts the results into a format similar to that returned by prcomp. "
test_pca3 <- PCAtools::pca(X, transposed = T)                                 # there is an issue in the code. if !transposed then it...transposes
#test_pca3b <- BiocSingular::runPCA(X, 2)                                     # equivalent. I specify how many PCs to retain (i.e., 2)
Z_pca3 <- as.matrix(test_pca3$rotated)                                        # scores. From findMethods("runPCA") out$x <- sweep(svd.out$u, 2, svd.out$d, "*") and rotated = data.frame(pcaobj$x), 
V_pca3 <- as.matrix(test_pca3$loadings)                                       # loadings or (right) eigenvectors: From findMethods("runPCA") out$rotation <- svd.out$v and loadings = data.frame(pcaobj$rotation)
stDev_pca3 <- test_pca3$sdev                                                  # From findMethods("runPCA") sdev = svd.out$d/sqrt(nrow(x) - 1)) which is our Eq.7
sv_pca3 <- test_pca3$sdev*sqrt(nrow(Y)-1)                                     # available internally but nor returned, hence had to be reproduced

## (B) compare with from scratch SVD PCA (sv not returned)
isEqualUnname(Z_pca3, my_Z_svd)
isEqualUnname(V_pca3, my_V)
isEqualUnname(sv_pca3, my_sv)                                                 ## based on Eq.7 ths is the stDev of scores
isEqualUnname(stDev_pca3^2, my_lambda)

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
isEqualUnname(apply(Z_pca3, 2, var), stDev_pca3^2)                            # verifies Eq.7 i.e., sv/sqrt(n - 1) = sqrt(eigenvalues of S) 
isEqualUnname(stDev_pca3, my_sv/sqrt(nrow(Y)-1))                              # verifies Eq.5 i.e., sv/sqrt(n - 1) = St.Dev of scores 


## ---- B.2b [BENCHMARK] - PCAtools, biplot() ####

## (D.1) Biplot coordinates: the function biplot does not return them
## Generates a scores and loadings, not a biplot but features are scaled "cosmetically"; coordinates not returned 
PCAtools::biplot(test_pca3, showLoadings = TRUE)
PCAtools::biplot(test_pca3, showLoadings = TRUE)$data                         ## Cannot get output but can get what it sees as "input" i.e. the PCA scores

##(D.2) Replication of biplot coordinates (call own function)
## -- First reproduce coordinates for a PC biplot (to show discrepancy)
biplot_test_pca3 <- replica_biplotPCAtools(test_pca3, showLoadings = TRUE)
biplot_observ_coord_pca3 <- as.matrix(biplot_test_pca3$observ_coord)
biplot_feat_coord_pca3 <- as.matrix(biplot_test_pca3$feat_coord)
isEqualUnname(biplot_observ_coord_pca3, my_Z_svd[,1:2])        ## it's a scores and loadings biplot
isEqualUnname(biplot_feat_coord_pca3, my_V[,1:2])              
isEqualUnname(biplot_test_pca3$feat_coord_scales, my_V[,1:2])  ## actually visualised features' coordinates, scaled by some constant for visualisation purposes



## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only)
## -- omitted:  not a PC biplot 

## (E.2) test 2: correlation between features
## -- omitted:  not a PC biplot

## (E.3) test 3:correlation between scores and features
pca3_test_featPC <- corPCTest(V = V_pca3,
                              Z = Z_pca3,      
                              Y = Y, 
                              lambda = stDev_pca3^2
) 
pca3_test_featPC$testS
pca3_test_featPC$check_eq19
pca3_test_featPC$corr_feat_PC


## (E.4) test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- omitted: not a PC biplot


## ---- B.3  [BENCHMARK] - ggbiplot           ####


## PART 1: get_SVD (illustrated with prcomp object. Accepts princomp, dudi.pca and more) 
## (A) call pre-built function
test_pca4 <- ggbiplot::get_SVD(test_pca1)
V_pca4 <- test_pca4$V
D_pca4 <- test_pca4$D                                                   # issue: not singular values -- see later why / see commentary in replica of this function
U_pca4 <- test_pca4$U                                                   # issue: not left eigenvectors -- see later why / see commentary in replica of this function
Z_pca4 <- test_pca4$U %*% test_pca4$D                                   # alternatively we could use the INPUT, but this is correctly deduced from how get_SVD obtains U

## (B) compare with from scratch SVD PCA (sv not returned)
isEqualUnname(Z_pca4, my_Z_svd)                                        # fails -- see later why / see commentary in replica of this function
isEqualUnname(V_pca4, my_V)                                            # success, loadings or right eigenvectors
isEqualUnname(U_pca4, my_U)                                            # fails -- see later why / see commentary in replica of this function
isEqualUnname(D_pca4, my_D)                                            # fails -- see later why / see commentary in replica of this function

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
## -- this might seem a bit too strict but is based on the underlying computational logic of get_SVD
## -- another option is to leave these checks to the PCA function used as input
isEqualUnname(apply(Z_pca4,2,var), my_lambda_svd)                      # Eq.5, fail
isEqualUnname(stDev_pca1, apply(test_pca1$x,2,sd))                     # Eq.7, fail

## (D) Replication to explain issue with "U"
# -- replicate how U is actually computed (incorreclty):
# -- starting from Z as given by prcomp, U = ZD{-1}*sqrt(n) because it knows that D = sqrt(eigenvalues) instead of D = diag(singular values)
U_pca4_replica <- test_pca1$x %*% solve(diag(test_pca1$sdev)*sqrt(nrow(Y)))      
isEqualUnname(U_pca4_replica, U_pca4)



## PART 2): ggbiplot
## (A) call pre-built function
## -- caveat: the output df.v is not quite what is computed in ggbiplot; "plot_env" is a ggplot thing that probably scales things further
test_pca4_biplot <- ggbiplot::ggbiplot(test_pca1, scale = 1, pc.biplot = FALSE)              # pc biplot
test_pca4_biplot_sl <- ggbiplot::ggbiplot(test_pca1, scale = 0, pc.biplot = FALSE)           # scores and loadings

## (B) compare with from scratch SVD PCA - same as get_SVD
isEqualUnname(test_pca4_biplot$plot_env$df.u, my_biplot_observ_coord[,1:2])                  # fails -- same reason as get_SVD
isEqualUnname(test_pca4_biplot$plot_env$df.v[,1:2], t(my_biplot_feat_coord)[,1:2])           # fails -- different from loadings given by get_SVD
isEqualUnname(test_pca4_biplot_sl$plot_env$df.u, my_Z[,1:2])                    # fails -- same reason as get_SVD
isEqualUnname(test_pca4_biplot_sl$plot_env$df.v[,1:2], my_V[,1:2]) 

## (C) Replication of biplot coordinates (call own function)
## run my replica to understand where these coordinates come from
pca4_biplot_replica <- replica_ggbiplot(test_pca1, pcbiplot_alpha = 1, pcbiplot_option = FALSE)
pca4_biplot_replica_sl <- replica_ggbiplot(test_pca1, pcbiplot_alpha = 0, pcbiplot_option = FALSE)

isEqualUnname(pca4_biplot_replica$observ_coord[,1:2], test_pca4_biplot$data)
isEqualUnname(pca4_biplot_replica$observ_coord[,1:2], test_pca4_biplot$plot_env$df.u)
isEqualUnname(pca4_biplot_replica$feat_coord_processed[,1:2], test_pca4_biplot$plot_env$df.v[,1:2])          

## (D) introduce approximation: I cannot fully replicate df.v there is always a constant that I cannot pick from the script
gg_correction_factor <- unique( c(as.matrix(pca4_biplot_replica$feat_coord_processed[,1:2]) / as.matrix(test_pca4_biplot$plot_env$df.v[,1:2])))
gg_correction_factor
isEqualUnname(pca4_biplot_replica$feat_coord_processed[,1:2]/gg_correction_factor, test_pca4_biplot$plot_env$df.v[,1:2])


## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only)
## -- caveat 1: the test requires all the dimensions (in all cases, the function does not produce a covariance matrix for Y)
## -- caveat 2: SVD functions do not compute a covariance matrix for Y
B_trsp_gg <- t(pca4_biplot_replica$feat_coord_processed/gg_correction_factor)
lambda_gg <- D_pca4^2                                              # if prcomp or princomp are used

## -- Test 1: features vector length: fails (biplot feat coord are multiplied by a "mystery" constant)
gg_test_len_orig <- lenTest(B_trsp_gg, Y = Y)              
gg_test_len_orig$testS
gg_test_len_orig$test_len
gg_test_len_orig$tab_out

## -- Test 2: correlation between features
gg_test_corcostest_orig <- corcosTest(B_trsp_gg, Y)
gg_test_corcostest_orig$coscorr_test
gg_test_corcostest_orig$coscorr_export_tab

## -- Test 3: correlation between features and PCs 
gg_test_featPC <- corPCTest(V = V_pca4, Z = Z_pca4, Y = Y, lambda = lambda_gg)
gg_test_featPC$check_eq19
gg_test_featPC$corr_feat_PC


## --- Test 4: Mahalanobis distance
## CAVEAT: the package doesn't provide a covariance matrix. We force the "correct" S
gg_test_mahala_orig <- MahalaTest(A = test_pca4_biplot$plot_env$df.u, Y = Y, S = cov(Y))          
gg_test_mahala_orig$mahalaTest   
gg_test_mahala_orig$dist_tab
  



##                                           ####
## ---- C benchmark vs eigendecomp packages  ####
## ---- C.0 Load Packages (for benchmark)    -----------------------------------------------------
library(ade4)
library(amap)
library(psych)


## ---- C.1 [BENCHMARK] - ade4, dudi.pca()   ####
## View(ade4::as.dudi)   # this is where the call to princomp happens, and where the biplot coordinates are computed
test_pca5 <- ade4::dudi.pca(X, center = TRUE, scale = FALSE, nf = ncol(X), scannf = FALSE)
Z_pca5 <- as.matrix(test_pca5$li)                                            # scores
V_pca5 <- as.matrix(test_pca5$c1)                                            # loadings
lambda_pca5 <- test_pca5$eig                                                 # eigenvalues or scores' variance
stDev_pca5 <- sqrt(lambda_pca5)                                              # scores' St Dev ASSUMED to be sqrt(eigenvalues) ?

## The function also provides biplot coordinates that are passed on to scatter.dudi
## -- PC Biplot coordinates consistently with Gabriel 1971
biplot_feat_coord_pca5 <-  test_pca5$co 
biplot_obsevr_coord_pca5 <- test_pca5$l1

## -- Scores and loadings biplot (Permute = TRUE in scatter.dudi)
biplot_sl_pca5_observ <- test_pca5$li                     # get scores
biplot_sl_pca5_feat <- test_pca5$c1                       # get loadings

## compare with "from scratch" eigendecomposition PCA
isEqualUnname(my_Z, Z_pca5)                          # due to issues similar to princomp (n instead of n-1; symmetric in eigen and weird sample cov matrix)
isEqualUnname(my_V_eig, V_pca5)                      # due to issues to princomp (n instead of n-1; symmetric in eigen and weird sample cov matrix)
isEqualUnname(my_lambda, lambda_pca5)                # false

## sign issue: ok in absolute values
isEqualAbs(my_Z, Z_pca5)
isEqualAbs(my_V_eig, V_pca5)

## Run my replica to retrieve some extra info (should get exact results as the original but allows m)
my_dude <- replica_dude(X, original = T, test_symmetric = T) 
my_dude$test_cov_orig
my_dude$test_myprin                                              # flag: if true, I got the same results with princomp (as replicated by me)
isEqualUnname(my_dude$Z, test_pca5$li)                           # check scores are excatly the same
isEqualUnname(my_dude$V, as.matrix(test_pca5$c1))                # check loadings are exaclty the same
isEqualUnname(my_dude$ev, test_pca5$eig)                         # check eigenvalues are exactly the same

## get the covariance matrix and what they use as "singular values"
S_pca5 <- my_dude$internal_cov
sv_pca5 <- my_dude$sv                                ## what they call dval, in reality it's not sv but just sqrt(ev)

## Test Eq. 5 and 7 (should fail)
isEqualUnname(apply(Z_pca5,2,var), lambda_pca5)      ## checking Eq. 5
isEqualUnname(sv_pca5, my_sv / (nrow(Y) - 1))        ## checking Eq. 7

## Now apply "fixes" and try again
my_dude_fix <- replica_dude(X, original = F, test_symmetric = F)
isEqualUnname(my_dude_fix$Z, my_Z)                           # check scores are excatly the same
isEqualUnname(my_dude_fix$V, my_V_eig)                       # check loadings are exaclty the same
isEqualUnname(my_dude_fix$ev, my_lambda)                     # check eigenvalues are exactly the same
isEqualUnname(apply(my_dude_fix$Z,2,var), my_dude_fix$ev)    # checking Eq. 5
isEqualUnname(my_dude_fix$sv, my_sv)                         # checking Eq. 7



## ---- C.1a[BENCHMARK] - ade4, scatter()    ####
## View(ade4:::scatter.dudi); assume Permute = FALSE (defaul) in scatter.dudi
## The PC biplot coordinates are incorrect because they use the wrong sv (also scores and loadings are different)
# compare with my grid's results: fail
isEqualUnname(as.matrix(biplot_obsevr_coord_pca5), my_biplot_observ_coord)                                                    
isEqualUnname(as.matrix(biplot_feat_coord_pca5), t(my_biplot_feat_coord))

# compare with princomp: also fails (ade4 doesn't use the singular values but just the square root of the eigenvalues)
isEqualUnname(biplot_test_pca10$obser_coord, as.matrix(biplot_obsevr_coord_pca5))
all.equal(unname(biplot_test_pca10$feat_coord),
          unname(as.matrix(biplot_feat_coord_pca5)))

## Now the biplot (no computations, visualisation only)
## -- caveat: permute = FALSE gives the right PC biplot BUT the observations are represented as arrows and the features as points...(wrong)
## -- caveat: permute = TRUE put the arrows correctly to represent features but is not a PC biplot
test_scatter_dudi_PCbiplot <- ade4:::scatter.dudi(test_pca5)                   # uses c0 and l1: PC biplot  ...but arrows are wrongly displayed       
test_scatter_dudi_sl <- ade4:::scatter.dudi(test_pca5, permute = TRUE )        # uses li and c1 i.e., scores and loadings





## test 1: features' vector length
## -- caveat: this test requires all the dimensions
B_trsp_pca5 <- t(as.matrix(biplot_feat_coord_pca5))
pca5_test_len_orig <- lenTest(B_trsp_pca5, Y, S = S_pca5)         # fails??
pca5_test_len_orig$testS
pca5_test_len_orig$test_len
pca5_test_len_orig$tab_out


## the following quantities should be equal
sqrt(diag(S_pca5))
apply(Y,2,sd) 
pca5_test_len_orig$length_B/sqrt(nrow(Y)-1)  # pca5_test_len_orig$length_B/sqrt(nrow(Y)-1)


## Properties test 2: correlation between features (pass)
pca5_test_corcostest_orig <- corcosTest(t(biplot_feat_coord_pca5), Y)
pca5_test_corcostest_orig$coscorr_test


## Properties test 3: correlation between scores and features (fail?)
pca5_test_featPC_1 <- corPCTest(V = V_pca5, 
                                 Z = Z_pca5, 
                                 Y = Y, 
                                 lambda = lambda_pca5,  
                                 S = S_pca5)
pca5_test_featPC_1$testS
pca5_test_featPC_1$check_eq19


## Test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- caveat euclidean and mahalanobis dist are identical instead of proportional (i.e. not as in Eq.20)
pca5_test_mahala <- MahalaTest(A = biplot_obsevr_coord_pca5, Y = Y, S = S_pca5)
pca5_test_mahala$mahalaTest                                                                       # fails when we consider its actual covariance matrix
pca5_test_mahala$dist_tab




## ---- C.2 [BENCHMARK] - amap, acp()        ####

## for details: amap:::acp
## -- Uses eigendecomposition
## -- from the documentaiton: "offer a variant of princomp and prcomp functions, with a slightly different graphic representation"
test_pca6 <-  amap:::acp(X, reduce = FALSE)

Z_pca6 <- test_pca6$scores                                  # scores
V_pca6 <- test_pca6$loadings                                # loadings or eigenvecotrs
lambda_pca6 <- test_pca6$eig^2                              # the name "eig" is misleading:these these the sqrt of the eigenvalues
stDev_pca6 <- test_pca6$sdev                                # scores St.Dev. computed directly from Z: if Eq.5 holds, these should be identical to the sqrt(eigenvlues) - but aren't

## compare with "from scratch" eigendecomposition PCA
isEqualUnname(my_Z, Z_pca6)                                 # fails due to covariance matrix
isEqualUnname(my_V_eig, V_pca6)                             # fails due to covariance matrix
isEqualAbs(my_Z, Z_pca6)                                    ## sign issue: ok in abs value
isEqualAbs(my_V_eig, V_pca6)                                ## sign issue: ok in abs value

## verify Eq. 5 and 7
isEqualUnname(apply(Z_pca6,2,var), lambda_pca6)             # Eq 5 doesn't work
isEqualUnname(apply(Z_pca6,2,sd), stDev_pca6)               # St Dev of scores computed directly...
isEqualUnname(test_pca6$eig, stDev_pca6)                    # ...but differs from eig = sqrt(eigenvalues)
isEqualUnname(test_pca6$eig, my_sv/sqrt(nrow(Y)-1))         # Eq.7 fails - RHS sqrt (eigenvalues); LHS features the singular values from scratch
isEqualUnname(test_pca6$eig, my_sv)                         # Eq.7 works without the denominator meaning Sqrt(eigenvectors) actually equals the singular values. This is due to choice of cov matrix

## Run my replcia to get the covariance matrix
pca6_replica <- replica_acp(X, original = T)
S_pca6 <- pca6_replica$cov_mat

## Check with "fixes"
pca6_replica_fix <- replica_acp(X, original = F)
isEqualUnname(my_Z, pca6_replica_fix$Z)
isEqualUnname(my_V_eig, pca6_replica_fix$V)
isEqualUnname(apply(pca6_replica_fix$Z,2,var), pca6_replica_fix$lambda)     # Eq 5 works

          




## ---- C.2a[BENCHMARK] - amap, plot()       ####
## scores and loadings bipot, not PC biplot. Does not return computations
amap:::plot.acp(test_pca6)                                                      # calls amap:::biplot.acp under the hood to overlay a plot of features

## reconstruction based on View(amap:::biplot.acp) and View(amap:::plot.acp)
biplot_observ_coord_pca6_sl <- test_pca6$scores                              # called U in the original script. Plotted as arrows
biplot_feat_coord_pca6_sl <- test_pca6$loadings                              # ...also called U in the original script. Plotted as arrows

## check with from scratch
isEqualAbs(biplot_observ_coord_pca6_sl, my_Z)
isEqualAbs(biplot_feat_coord_pca6_sl, my_V_eig)

## Properties test 1: features' vector length
## -- caveat: won't work because it's not a PC biplot
pca6_test_len <- lenTest(t(as.matrix(biplot_feat_coord_pca6_sl)), Y, S = S_pca6)         
pca6_test_len$testS
pca6_test_len$test_len
pca6_test_len$tab_out


## Properties test 2: correlation between features 
## -- caveat: won't work because it's not a PC biplot
pca6_test_corcostest <- corcosTest(t(biplot_feat_coord_pca6_sl), Y)
pca6_test_corcostest$coscorr_test
pca6_test_corcostest$coscorr_export_tab


## Properties test 3: correlation between scores and features (fail)
pca6_test_featPC <- corPCTest(V = V_pca6, Z = Z_pca6, Y = Y, lambda = lambda_pca6, S = S_pca6 )
pca6_test_featPC$testS
pca6_test_featPC$check_eq19
pca6_test_featPC$corr_feat_PC


## Properties test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- caveat: won't work because it's not a PC biplot
pca6_test_mahala_orig <- MahalaTest(A = biplot_observ_coord_pca6_sl, Y = Y, S = S_pca6)          
pca6_test_mahala_orig$mahalaTest                                                                    






## ---- C.3 [BENCHMARK] - psych, principal() ####
## it computes an eigen decomposition of the covariance matrix, but WON'T centre the data so input must be Y
test_pca7 <- psych::principal(Y, nfactors = 0, cor = "cov", rotate = "none")

## Caveat: the values returned are normalised by default (the process cannot be undone)
## -- From the help for this function: "...The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings"
Z_pca7_norm <- test_pca7$scores[,]                             # scores BUT normalised "as if" they were biplot's observ coords [Eq.11 in paper except uses the feature's sd instead of singular values]
V_pca7_norm <- test_pca7$loadings[,]                           # loading BUT normalised "as if" they were biplot's feat coords [Eq.12 in paper, but with sd instead of singular values]
lambda_pca7 <- test_pca7$values                                # eigenvalues - this is OK

## quick fix to reverse built-in normalisation process 
Z_pca7 <- sweep(Z_pca7_norm[,], 2, sqrt(lambda_pca7), "*")      # still not the same as Z due to some signs. Compare [Eq.11 in paper] 
V_pca7 <- sweep(V_pca7_norm[,], 2, sqrt(lambda_pca7), "/")      # still not the same as V_eig due to some signs. [Eq.12 in paper] but with score's sd instead of singular values
isEqualUnname(Z_pca7, my_Z)
isEqualUnname(V_pca7, my_V_eig)
isEqualUnname(my_lambda, lambda_pca7)
isEqualAbs(Z_pca7, my_Z)
isEqualAbs(V_pca7, my_V_eig)

## Test Eq.5 & 7
isEqualUnname(apply(Z_pca7_norm,2,var), lambda_pca7)           # Eq.5 fails as is
isEqualUnname(apply(Z_pca7,2,var), lambda_pca7)                # Eq.5 ok but requires DE-WEIGHTING to work, otherwise it won't
isEqualUnname(sqrt(lambda_pca7), my_sv/sqrt(nrow(Y)-1))        # Eq.7 ok but only if we use MY sv. This function uses sqrt(ev) as a proxy for sv but not explicitly


## Run my replica to (1) get covariance matrix; (2) grap issues with signs
## -- to understand how I replicate results and signs, see my function script for details
replica_pca7 <- replica_psych(Y, original = T, signs_thing = T)
isEqualUnname(replica_pca7$scores, test_pca7$scores[,])
isEqualUnname(replica_pca7$loadings, test_pca7$loadings[,])                


## -- get covariance matrix
S_pca7 <- replica_pca7$cov_mat                                     # retrieve covariance matrix for next section



## ---- C.3a[BENCHMARK] - psych, biplot()    ####
## based on the scrip View(psych::biplot.psych)

## don't be fooled by the names: these are scaled just like you would do in a biplot
biplot_test_pca7_observ <- test_pca7$scores[,]
biplot_test_pca7_feat <- test_pca7$loadings[,] 

## won't match our results: recall this is a strange "hybrid" set of coordiantes, rather unconventional
isEqualAbs(biplot_test_pca7_observ, my_biplot_observ_coord_eig)          
isEqualAbs(biplot_test_pca7_feat, my_biplot_feat_coord_eig)

## quick fix to show difference from a PC biplot
biplot_test_pca7_observ_fixed <- Z_pca7_norm * 1/sqrt(nrow(Y) - 1)
biplot_test_pca7_feat_fixed <- V_pca7_norm * sqrt(nrow(Y) - 1)
isEqualAbs(biplot_test_pca7_observ_fixed, my_biplot_observ_coord_eig)
isEqualAbs(biplot_test_pca7_feat_fixed, t(my_biplot_feat_coord_eig))



## Test geometrical properties (feature properties as in PC biplot)
B_trsp_pca7 <- t(biplot_test_pca7_feat)                                  # useful for testing

## -- Test 1: feautre's vector length: won't work as is
pca7_test_len <- lenTest(B_trsp_pca7, Y, S = S_pca7)         
pca7_test_len$testS
pca7_test_len$test_len
pca7_test_len$tab_out


## -- Test 2: correlations between features
pca7_test_corcostest <- corcosTest(B_trsp_pca7, Y)           # works even with incorrect length??
pca7_test_corcostest$coscorr_test
pca7_test_corcostest$coscorr_export_tab


## -- Test 3: correlation between features and PCs
## Obviously won't work with the "as is" loadings and scores because of the internal manipulations
pca7_test_featPC <- corPCTest(V = V_pca7_norm,
                              Z = Z_pca7_norm, 
                              Y = Y, 
                              lambda = lambda_pca7,
                              S = S_pca7
) 
pca7_test_featPC$check_eq19
pca7_test_featPC$corr_feat_PC


## Test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## -- fails with the "as is" output
pca7_test_mahala <- MahalaTest(A = biplot_test_pca7_observ, Y = Y, S = S_pca7)          
pca7_test_mahala$mahalaTest                                                                     
pca7_test_mahala$dist_tab
                                                            



##                                           ####
## ---- D benchmark vs generalised SVD       ####
## ---- D.0 load packages (for benchmark)    ####
library(FactoMineR)
library(PCAmixdata)
library(vegan)
## ---- D.1a [BENCHMARK] - FactoMineR's PCA() ####
test_pca8 <- FactoMineR::PCA(X, scale.unit = FALSE)               # View source script: View(FactoMineR::PCA)
sv_pca8 <- test_pca8$svd$v                                        # singular values. However, these are not the usual (see later)
D_pca8 <- diag(as.vector(test_pca8$svd$vs))                       # diagnosalised singular values
lambda_pca8 <- test_pca8$eig[,1]                                  # in the original script this is just vs^2
V_pca8 <- test_pca8$svd$V
U_pca8 <- test_pca8$svd$U
biplot_observ_coord_pca8 <- test_pca8$ind$coord                   # biplot, features coordinates
biplot_feat_coord_pca8 <-test_pca8$var$coord                      # biplot, observations coordiantes
Z_pca8 <- biplot_observ_coord_pca8                                # my addition: no scores as such in the original script

## biplot coordinates (seems and odd "hybrid" betewen a scores and loadings biplot and a PC biplot...)
isEqualUnname(biplot_feat_coord_pca8, t(t(V_pca8) * sv_pca8) )    # from original script. this is B = VD  (I use B transposed = DV')
isEqualUnname(biplot_observ_coord_pca8, t(t(U_pca8) * sv_pca8))   # from original script. equivalent to Z = UD and A = Z (alpha = 1). But scores Z are not provided as such


## Compare SVD output with from scratch: all fail (more than just sign difference)
isEqualUnname(lambda_pca8, my_lambda_svd)
isEqualUnname(sv_pca8, my_sv)
isEqualUnname(U_pca8, my_U)
isEqualUnname(V_pca8, my_V)
isEqualUnname(Z_pca8, my_Z_svd)
isEqualUnname(biplot_observ_coord_pca8, my_biplot_observ_coord)
isEqualUnname(biplot_feat_coord_pca8, my_biplot_feat_coord)

## try absolute values
isEqualAbs(V_pca8, my_V)
isEqualAbs(U_pca8, my_U)
isEqualAbs(Z_pca8, my_Z_svd)
isEqualAbs(biplot_observ_coord_pca8, my_biplot_observ_coord)
isEqualAbs(biplot_feat_coord_pca8, my_biplot_feat_coord)


## test scores variance and eigenvalues vs sv
isEqualUnname(apply(Z_pca8, 2, var), lambda_pca8)                 # Eq. 5 fails if I use the eigenvalues returned by the script
isEqualUnname(apply(Z_pca8, 2, var), my_lambda_svd)               # oddly correct if I use the correctly computed eigenvalues
isEqualUnname(sv_pca8/sqrt(nrow(Y)-1), sqrt(lambda_pca8))         # Eq. 7 fails
isEqualUnname(sv_pca8, sqrt(lambda_pca8))                         # sv in Factominer's svd is sv/sqrt(n-1) in my framework ?


## Here I replicate Factorminer's "generalised svd" with some fixes
## -- All differences ultimately due to (1) svd of Y / sqrt(n) instead of Y; (2) internal tweaking with signs. 
svd_pca8 <- FactoMineR::svd.triplet(Y)                                       # this "generalized" SVD is what the feeds the function PCA()
svd_pca8_replica <- replica_svdTriplet(Y, original = T, sign_leave = F)      # this is just svd(Y / sqrt(n)) see my own function for detail
svd_pca8_replica_fix <- replica_svdTriplet(Y, original = T, sign_leave = T)  # same except I "switch off" the messing around with signs to ease comparison

## re-check U, V and sv now that signs have been fixed
isEqualUnname(my_sv/sqrt(nrow(Y)), svd_pca8_replica_fix$sv)                  # singular values = my singular values divided by sqrt(n)
isEqualUnname(my_U*sqrt(nrow(Y)), svd_pca8_replica_fix$U)                    # left eigenvectors = my U times sqrt(n). Similar to Gabriel 1971 eq.52
isEqualUnname(my_V, svd_pca8_replica_fix$V)                                  # V is the same

## now re-check factominer's biplot coordinates using the previous equivalences
biplot_feat_coord_pca8_2 <- t(t(my_V) * my_sv/sqrt(nrow(Y)))                  # B = VD/sqrt(n)  (I use B transposed = DV')
biplot_observ_coord_pca8_2 <- t(t(my_U*sqrt(nrow(Y))) * my_sv/sqrt(nrow(Y)))  # Z = (U*sqrt(n))*(D/sqrt(n)) = UD
isEqualUnname(biplot_observ_coord_pca8_2, my_Z_svd)                            # observ coord = Z as in scores and loadings biplot
isEqualUnname(biplot_feat_coord_pca8_2, t(my_biplot_feat_coord)/sqrt(nrow(Y))) # feat coord = as in PC biplot but divided by sqrt(n)



## ---- D.1b [BENCHMARK] - FactoMineR biplot  ####
## based on View(FactoMineR::plot.PCA). Biplot coordinates are pased on from the PCA function
A_pca8 <- biplot_observ_coord_pca8                                    # scores, as in a scores and loadings biplot
B_trsp_pca8 <- t(biplot_feat_coord_pca8 )                             # as in a PC biplot (check Eq. 12 in paper)

## Test geometrical properties (feature properties as in PC biplot)
## -- Test 1: feature's vector length
pca8_test_len <- lenTest(B_trsp_pca8, Y, S = NULL)                    # because it's SVD-based there is no covariance matrix produced
pca8_test_len$testS
pca8_test_len$test_len
pca8_test_len$tab_out

## -- Test 2: correlations between features
pca8_test_corcostest <- corcosTest(B_trsp_pca8, Y)
pca8_test_corcostest$coscorr_test
pca8_test_corcostest$coscorr_export_tab

## -- Test 3: correlation between features and PCs
## ----- this is the only function to provides these correlations, but the test fails 
pca8_test_featPC <- corPCTest(V = V_pca8,
                                  Z = Z_pca8, 
                                  Y = Y, 
                                  lambda = lambda_pca8, 
                                  ) 
pca8_test_featPC$check_eq19
isEqualUnname(c(pca8_test_featPC$corr_feat_PC$corr), c(test_pca8$var$cor))    # my test replicates the correlations correctly, but the test fails






## ---- D.2a [BENCHMARK] - vegan              ####
library(vegan)
## The script itself is very difficult to navigate
## -- https://rdrr.io/cran/vegan/src/R/pca.R
## -- https://rdrr.io/cran/vegan/src/R/eigenvals.R
## A very useful resource: https://ourcodingclub.github.io/tutorials/ordination/

## (A) call pre-built function 
test_pca_x1 <-  vegan::rda(X, scale = FALSE)                         # performs svd of  Y/sqrt(nrow(Y) - 1)
scores_pca_x1 <- vegan::scores(test_pca_x1, choices = 1:ncol(Y) )    # the closest we get to PCA scores is the "sites" from the function "scores"

## (B) compare with from scratch SVD PCA (sv not returned)
isEqualUnname(test_pca_x1$CA$v, my_V)
isEqualUnname(test_pca_x1$CA$u, my_U)
isEqualUnname(test_pca_x1$CA$eig, my_lambda_svd)                    # eigenvalues (of S) correct even if this this an SVD


## (C) Issue: no "scores" in Vegan: in fact, sites and scores are PC biplot coordinates  
## (C.1) replicate vegan's svd and scores: 
svd_veg_test <- replica_svd_vegan(X)                                # reproduce svd of  Y/sqrt(nrow(Y) - 1)
vegan_biplot_replica <- replica_vegan_scores_FAIL(test_pca_x1)      # reproduce "scores" as "sites" - but didn't manage

## check if I replicated correctly
isEqualUnname(svd_veg_test$u, test_pca_x1$CA$u)                     # succeeds
isEqualUnname(svd_veg_test$v, test_pca_x1$CA$v)                     # succeeds
isEqualUnname(svd_veg_test$d^2, test_pca_x1$CA$eig)                 # this means d from svd( Y/sqrt(nrow(Y) - 1) ) is sqrt(eigenvalues) 
isEqualUnname(svd_veg_test$d, my_sv/sqrt(nrow(Y)-1))                # this means d from svd are divided by sqrt(n-1) = sqrt(eigenvalues) because of Eq.7
isEqualUnname(vegan_biplot_replica$sites_replica, scores_pca_x1$sites)      # fail to replicate exactly, but my replica bring useful info nevertheless


## (C.2) try theoretical scores from scratch
Z_pca_x1 <- Y %*% test_pca_x1$CA$v                                          # manual computation of scores - how it should be (but isn't)
isEqualUnname(Z_pca_x1, my_Z_svd)                                           # success

# now we have 3 scores - not too far off from each other but not quite the same either
head(Z_pca_x1)
head(vegan_biplot_replica$sites_replica)                                    
head(scores_pca_x1$sites)
vegan_biplot_replica$sites_replica / scores_pca_x1$sites                    # every column  divided by a constant
scores_pca_x1$sites / Z_pca_x1                                              # every column  divided by a constant

## (C.3) try vegan "scores" vs biplot observations coordinates
pca_x1_sv <- svd_veg_test$d                                                 # sv, as computed internally i.e., not really sv: in fact this is sqrt(eigenvalues) but "technically"
A_x1 <- sweep(Z_pca_x1, 2, pca_x1_sv, "/")                                  # Eq.11: A = ZD^(-1) however obviously wrong because D = diag [sqrt(eigenvalue_k)]
isEqualUnname(A_x1, my_biplot_observ_coord)                                 # fail: not our PC biplot coordinates
mystery_scores_constant <- scores_pca_x1$sites / A_x1                       # A = ZD^{-1} * "one mystery constant"
mystery_scores_constant

## (C.4) try vegan "scores" vs our correct observations coordinates
pca_x1_sv_correct <- svd_veg_test$d*sqrt(nrow(Y) - 1)                           # based on Eq.7: sqrt(eigenvalues) = sv/sqrt(nrow(Y) - 1)
A_x1_correct <-  sweep(Z_pca_x1, 2, pca_x1_sv_correct, "/")                     # A = ZD^{-1} 
isEqualUnname(A_x1_correct, my_biplot_observ_coord)                             # success - same as  our PC biplot coordinates
isEqualUnname(pca_x1_sv_correct, my_sv)                                         # correct

scores_pca_x1$sites / A_x1_correct                                              # vegan's score = our correct matrix A * mystery constant 
t(my_biplot_feat_coord)/scores_pca_x1$species                                   #  Same goes for  "species"

sqrt(sqrt((nrow(Y) - 1) * sum(test_pca_x1$CA$eig)))                             #... mystery solved: this is  "const" in vegan's scores script
vegan_biplot_replica$costant                                                    # here it is, I saved it: see function's script


## chck identities
isEqualUnname(t(my_biplot_feat_coord)/vegan_biplot_replica$costant, scores_pca_x1$species)
isEqualUnname(my_biplot_observ_coord*vegan_biplot_replica$costant, scores_pca_x1$sites)

## (E) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
isEqualUnname(apply(scores_pca_x1$sites, 2, var), test_pca_x1$CA$eig)         # fails -- in fact these are not "scores" as we know them see above
isEqualUnname(apply(Z_pca_x1, 2, var), test_pca_x1$CA$eig)                    # works if I use MY scores. this is because eig is computed correctly from a "special" SVD, so this is not cheating
isEqualUnname(sqrt(test_pca_x1$CA$eig), pca_x1_sv/sqrt(nrow(Y) - 1))          # fails
isEqualUnname(sqrt(test_pca_x1$CA$eig), my_sv/sqrt(nrow(Y) - 1))              # works if I use MY sv



## ---- D.2b [BENCHMARK] - vegan biplot       ####

## (D.1) Biplot coordinates: the function biplot does not return them
## The coordinates of a biplot for package vegan come from the applying the function "scores()" to the output of "rda()"
## -- View(vegan:::scores)
## -- View(vegan:::biplot.rda)
test_vegan_biplot <- vegan:::biplot.rda(test_pca_x1)

biplot_features_coord_pca_x1 <- scores_pca_x1$species
biplot_observations_coord_pca_x1 <-  scores_pca_x1$sites

## relationship with "correct" biplot coordinates
## The result is a PC biplot in principle, but scaled: A = ZD^{-1} * constant; B' = DV' / constant
isEqualUnname(biplot_observations_coord_pca_x1, my_biplot_observ_coord*vegan_biplot_replica$costant)
isEqualUnname(biplot_features_coord_pca_x1,t(my_biplot_feat_coord)/vegan_biplot_replica$costant)


## (E) Geometrical properties tests
## (E.1) test 1: features' vector length (PC biplot only). fails (biplot feat coord are multiplied by a "mystery" constant)
pca_x1_test_len_orig <- lenTest(t(biplot_features_coord_pca_x1), Y)          # vegan doesn't provide a covariance matrix         
pca_x1_test_len_orig$testS
pca_x1_test_len_orig$test_len
pca_x1_test_len_orig$tab_out

## (E.2) test 2: correlation between features
pca_x1_test_corcostest_orig <- corcosTest(t(biplot_features_coord_pca_x1), Y)
pca_x1_test_corcostest_orig$coscorr_test
pca_x1_test_corcostest_orig$coscorr_export_tab

## (E.3) test 3: correlation between scores and features
## CAVEAT: vegan does not provide scores as we know them
pca_x1_test_featPC <- corPCTest(V = test_pca_x1$CA$v,
                                Z = scores_pca_x1$sites,      # not quite the scores as we know them
                                Y = Y, 
                                lambda = test_pca_x1$CA$eig 
) 
pca_x1_test_featPC$check_eq19
pca_x1_test_featPC$corr_feat_PC


## (E.4) test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
## CAVEAT: vegan is svd based so it doesn't provide a covariance matrix. We force the "correct" S
pca_x1_test_mahala_orig <- MahalaTest(A = biplot_observations_coord_pca_x1, Y = Y, S = cov(Y))          
pca_x1_test_mahala_orig$mahalaTest                                                                     
pca_x1_test_mahala_orig$dist_tab
                                                                  





## ---- D.4  [BENCHMARK] - PCAmix             ####

## Notes
## -- like package FactoMiner and vegan PCAmix performs an svd on a specially weighted Y i.e. Y*(sqrt(n/n-1))
## -- unlike other packages it forces normalisation besides centering 

## (A) call pre-built function 
test_pca9 <- PCAmixdata::PCAmix(X.quanti = X,
                                ndim = ncol(X),
                                )
Y_weight_pca9 <- test_pca9$Z[,]                         # the matrix being decompoosed, as returned by the function
V_pca9 <- test_pca9$V                                   # loadings are not what you think (if you read the documentaiton $coef reads like loadings but it isn't)
lambda_pca9 <- test_pca9$eig[,1]                        # just like the vegan package, sqrt(eigenvalue) should comply with Eq.7
Z_pca9 <- test_pca9$scores                              # Scores are not what you think:


## (B) Run my replica (PCAmix is tricker tha others so replication is done early)
pca9_replica <- replica_svd_PCAmix(X, original = T)                            # retrieve the weighted matrix returned by Z in the original function
pca9_replica_fixes <- replica_svd_PCAmix(X, original = F, sign_leave = T)      # switch off the built-in normalisation, and swith off the messing about with signs. generates output comparable with from scratch results
sv_pca9 <- pca9_replica$d

## -- show that replication worked
isEqualUnname(Y_weight_pca9, pca9_replica$Y_weigh)
isEqualUnname(V_pca9, pca9_replica$v)
isEqualUnname(lambda_pca9, pca9_replica$d^2)
isEqualUnname(Z_pca9, pca9_replica$Y_weigh %*% pca9_replica$v)

## -- show that weighted matrix being decomposed is Y* (= centre and NORMALISE X) times sqrt(nrow(X)/(nrow(X) - 1)),  see annotations in my function
Y_star <- scale(X, center = T, scale = T)[,]
isEqualUnname(Y_star*sqrt(nrow(Y)/(nrow(Y) - 1)), test_pca9$Z[,]) 


## (C) compare with from scratch SVD PCA 
## -- show that PCAmix can NOT be compared directly with from scratch results, unlike other packages
isEqualUnname(V_pca9, my_V)
isEqualUnname(lambda_pca9, my_lambda_svd)
isEqualUnname(Z_pca9, my_Z_svd)

## -- repeat comparison after switching off normalisation via our replica
isEqualUnname(pca9_replica_fixes$u, my_U*sqrt(nrow(Y)))           # SVD returns U multiplied by sqrt(n) - notice signs are ok too, because I acted on svd.triplet
isEqualUnname(pca9_replica_fixes$v, my_V)                         # V is almost unaffected - notice signs are ok too, because I acted on svd.triplet
isEqualUnname(pca9_replica_fixes$d, my_sv/sqrt(nrow(Y)-1))        # SVD of weighted data matrix returns sv divided by sqrt(n-1) = sqrt(eigenvalues)

eigenvalues_pca9_fix <- pca9_replica_fixes$d^2
sv_pca9_fix <- pca9_replica_fixes$d*sqrt(nrow(Y) - 1)             # manual computation of singular values Eq. 7 in my paper. Will be useful later           


## (D) test eq.5 and 7
isEqualUnname(apply(Z_pca9,2,var), lambda_pca9)
isEqualUnname(sqrt(lambda_pca9), sv_pca9/sqrt(nrow(Y) - 1))


## Biplot function: View(PCAmixdata:::plot.PCAmix)
## biplot coordinates: it's NOT a score and loadings biplot because Z and V are "scaled" already in a way similar to Gabriel's PC biplot decomposition
biplot_observ_coord_pca9 <- test_pca9$ind$coord                 # matrix A = Z in a PCA biolot 
biplot_feat_coord_pca9 <-  test_pca9$quanti$coord               # matrix B = VD in PCA biplot

## same hybrid approach as FactoMiner (similar to HJ biplots)
isEqualUnname(biplot_observ_coord_pca9, test_pca9$scores)                              # A = Z like in a scores and loadings biplot
isEqualUnname(t(biplot_feat_coord_pca9),  diag(sqrt(lambda_pca9)) %*% t(test_pca9$V))  # B' = DV' like in a PCA biplot

isEqualUnname(biplot_observ_coord_pca9, my_biplot_observ_coord)
isEqualUnname(biplot_feat_coord_pca9, my_biplot_feat_coord)


## Use fixes to understand biplot coordinates assuming we don't standardise Y
fixed_Z_pca9 <-  pca9_replica_fixes$Y_weigh %*%  pca9_replica_fixes$v 
fixed_biplot_observ_coord_pca9 <- fixed_Z_pca9
fixed_biplot_feat_coord_pca9 <-  diag(pca9_replica_fixes$d) %*% t(pca9_replica_fixes$v)

correction_factor <- sqrt(nrow(Y))/sqrt(nrow(Y)-1)

## now compare with grid
isEqualUnname(fixed_Z_pca9, my_Z_svd*correction_factor)
isEqualUnname(fixed_biplot_observ_coord_pca9, my_Z_svd*correction_factor)
isEqualUnname(fixed_biplot_feat_coord_pca9, my_biplot_feat_coord / sqrt(nrow(Y)-1))


## Geometric properties test
## -- Test 1: features vector length: fails (biplot feat coord are multiplied by a "mystery" constant)
pca9_test_len_orig <- lenTest(t(biplot_feat_coord_pca9), Y)          # PCA mix doesn't provide a covariance matrix         
pca9_test_len_orig$testS
pca9_test_len_orig$test_len
pca9_test_len_orig$tab_out

## -- Test 2: correlation between features
pca9_test_corcostest_orig <- corcosTest(t(biplot_feat_coord_pca9), Y)
pca9_test_corcostest_orig$coscorr_test
pca9_test_corcostest_orig$coscorr_export_tab

## -- Test 3: correlation between features and PCs
## CAVEAT: vegan does not provide scores as we know them
pca9_test_featPC <- corPCTest(V = test_pca9$V,
                              Z = test_pca9$scores,      
                              Y = Y, 
                              lambda = test_pca9$eig[,1] 
) 
pca9_test_featPC$testS
pca9_test_featPC$check_eq19
pca9_test_featPC$corr_feat_PC

## --- Test 4: Mahalanobis distance
## CAVEAT: the package doesn't provide a covariance matrix. We force the "correct" S
pca9_test_mahala_orig <- MahalaTest(A = biplot_observ_coord_pca9, Y = Y, S = cov(Y))          
pca9_test_mahala_orig$mahalaTest   
pca9_test_mahala_orig$dist_tab

##                                           ####
## ---- R benchmark vs latest biplot fun     ####
## ---- E.0 load packages                   ####
library(bpca)
library(biplotEZ)
library(BiplotGUI)
library(MultBiplotR)
## ---- E.1 [BENCHMARK] - bpca              #### 

## (A) call pre-built function 
test_bpca <- bpca(X, scale = F, method = "gh", var.rb = T)           # the "GH" method is equivalent to a PC biplot (alpha = 0 in Gabriel's factorisation)
test_bpca_sl <- bpca(X, scale = F, method = "jk", var.rb = T)        # the "JK" method is equivalent to a scores and loadings biplot (alpha = 1 in Gabriel's factorisation)
Z_bpca <- test_bpca_sl$coord$objects                                 ## Verify scores: this package does not generate PCA scores per se, but we can "force it" by choosing "JK" as a method (i.e. alpha = 1 in Gabriel's factorisation yields Z)
V_bpca <- test_bpca_sl$coord$variables
lambda_pca <- test_bpca$eigenvalues                                  # notice the incorrect terminology (these are singular values)
corr_feat_bpca <- test_bpca$var.rb                                   # correlation coefficients for all variables

## (B) compare with from scratch SVD PCA (sv not returned)
isEqualUnname(V_bpca, test_bpca$eigenvectors)                        # success
isEqualUnname(test_bpca$eigenvalues, my_lambda_svd)                  # fail. what they call eigenvalues is in fact sv
isEqualUnname(test_bpca$eigenvalues, my_sv)                          # success: what they call eigenvalues are in fact singular values
isEqualUnname(test_bpca$eigenvectors, my_V)                          # success

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
## -- wrong because unfortunately from the documentation "eigenvalues" is incorreclty defined as "A vector of the eigenvalues."
isEqualUnname(apply(Z_bpca,2,var), lambda_pca)                       # Eq 5
eigenvalues_correct <- (test_bpca$eigenvalues/(sqrt(nrow(Y)-1)))^2   # Eq 7 at is should be: what they call eigenvalutes are in fact sv 
isEqualUnname(apply(Z_bpca,2,var), eigenvalues_correct )             # Eq 5 holds with fixes

## (D.1) Biplot coordinates
biplot_coord_observ_bpca <- test_bpca$coord$objects
biplot_coord_feat_bpca <- test_bpca$coord$variables
isEqualUnname(test_bpca$coord$objects, my_biplot_observ_coord)       # fail because of Gabriel Eq.52 hard-wired additional scaling
isEqualUnname(test_bpca$coord$variables, my_biplot_feat_coord)       # fail because of Gabriel Eq.52 hard-wired additional scaling

## (D.2) Replication of biplot coordinates (call own function)
## -- see comment within function to understand results
my_bpca <- replica_bpca(X)

## Verify geometrical properties (for PC biplot version only)
## -- Test 1: features vector length: fails (biplot feat coord are multiplied by a "mystery" constant)
pca_bpca_test_len_orig <- lenTest(t(test_bpca$coord$variables), Y)          # vegan doesn't provide a covariance matrix         
pca_bpca_test_len_orig$testS
pca_bpca_test_len_orig$test_len
pca_bpca_test_len_orig$tab_out

## -- Test 2 
pca_bpca_cocostest <- corcosTest(t(test_bpca$coord$variables), Y)
pca_bpca_cocostest$coscorr_test

## -- Test 3 correlation between scores and features
pca_bpca_test_featPC <- corPCTest(V = test_bpca$eigenvectors, Z = Z_bpca, Y = Y, lambda = test_bpca$eigenvalues)
pca_bpca_test_featPC$check_eq19
pca_bpca_test_featPC$corr_feat_PC

## --- Test 4: Mahalanobis distance
pca_bpca_test_mahala_orig <- MahalaTest(A = test_bpca$coord$objects, Y = Y, S = cov(Y))          
pca_bpca_test_mahala_orig$mahalaTest 
pca_bpca_test_mahala_orig$dist_tab




## ---- E.2 [BENCHMARK] - EZ                ####
## Calibrated axis biplot as in Gower et al 2011, package written by some of Gower's co-authors
## -- underlying math: https://cran.r-project.org/web/packages/biplotEZ/vignettes/biplotEZ.html
## -- based on View(biplotEZ:::biplot) and  View(biplotEZ:::PCA.biplot)              

## (A) call pre-built function 
## A.1) display only (typical call - https://cran.r-project.org/web/packages/biplotEZ/vignettes/biplotEZ.html)
biplotEZ::biplot(X, scaled = FALSE) |>                                                      #  a bit misleading: centres and scales only 
  PCA(correlation.biplot = TRUE) |>                                                         # everything happens here
  plot()

## A.2) call with outputs (2 steps)
biplot_EZ_output <- biplotEZ::biplot(X, center = T, scaled = F)                             # this function centres and scales
biplot_EZ_pca_PCbiplot <- biplotEZ::PCA(biplot_EZ_output,                                   # PC biplot:  in principle A = U; B' = DV' but turns out ot be WRONG
                                correlation.biplot = TRUE,
                                e.vects = 1:ncol(X))
biplot_EZ_pca_sl <- biplotEZ::PCA(biplot_EZ_output,                                         # Scores and loadings biplot: in principle A = Z; B = V
                                  correlation.biplot = FALSE,
                                  e.vects = 1:ncol(X))
Y_EZ <- biplot_EZ_output$X[,]                                                               # centred and scaled data matrix. the only useful numerical output returned by biplot fun.
V_EZ <- biplot_EZ_pca_sl$Lmat	                                                              # actually feat coord but workaround to get V with all dimensions. Called "the matrix for transformation to the principal components"
Z_EZ_est <- (Y_EZ %*% V_EZ)                                                                 # scores not returned as such. ONLY IF PC BIPLOT, returned as Z but that's a coincidence
lambda_EZ <- biplot_EZ_pca_sl$eigenvalues                                                   # eigenvalues returned                                                  

## (B) compare with from scratch SVD PCA (sv not returned)
## -- first, default setting: scores and loadings biplot (observ. coord are just the PCA scores)
isEqualUnname(Y_EZ, my_Y)     
isEqualUnname(V_EZ, my_V)
isEqualUnname(Z_EZ_est, my_Z_svd)                                                           # internally computed only if a PC biplot is generated 
isEqualUnname(lambda_EZ, my_lambda_svd)                                                     # EZ eigenvalues = sv^2, NOT from S = cov(Y)
isEqualUnname(lambda_EZ, my_sv^2)                                                           # what they call eigenvalues  is just sv^2 from svd(Y)  

## (C) Check Eq.5 and 7 (own eigen-decomposition as benchmark)
isEqualUnname(apply(Z_EZ_est, 2, var), lambda_EZ)                                        
isEqualUnname(sqrt(lambda_EZ), sqrt(lambda_EZ)/(nrow(Y_EZ) - 1))                            # Eq.7 fails because what they call eigenvalues = sv^2 
isEqualUnname(sqrt(my_lambda_svd), sqrt(lambda_EZ)/(nrow(Y_EZ) - 1))                        # Eq.7 won't even work with my eigenvalues and their sv


## (D) Biplot coordinates - returned with only 2 dimensions by default (Cannot be changed) 
## (D.1) scores and loadings biplot coordinates 
EZ_observ_coord_sl_r2 <- biplot_EZ_pca_sl$Z                                                 # biplot observ coord matrix A = Z, first two columns of scores used as coordinates for observations (points)
EZ_feat_coord_sl_r2 <-   biplot_EZ_pca_sl$Vr                                                # biplot feat coord matrix B = V, first two columns of loadings used as coordinates for features 
EZ_observ_coord_sl <- Z_EZ_est                                                              # workaround: all dimensions
EZ_feat_coord_sl <- V_EZ	                                                                  # workaround: all dimensions

isEqualUnname(EZ_observ_coord_sl, my_Z_svd)
isEqualUnname(EZ_feat_coord_sl, my_V)                                                         
isEqualUnname(biplot_EZ_pca_PCbiplot$Z, my_biplot_observ_coord[,1:2])                       # fails: yet this should be A = ZD^{-1} but it isn't
isEqualUnname(biplot_EZ_pca_PCbiplot$Lmat, my_biplot_feat_coord[, 1:2])



# (D.2) Replication of biplot coordinates (call own function)
## -- First reproduce coordinates for a PC biplot (to show discrepancy), see annotations in my function
EZ_replica <- replica_biplotEZ(X, correlation.biplot = T, original = T)
S_EZ <- EZ_replica$S_EZ_replica                                                                  # get covariance matrix computed internally
isEqualUnname(EZ_replica$biplot_obscoord[,1:2], biplot_EZ_pca_PCbiplot$Z)                            # successful replication of what seems an incorrect result
isEqualUnname(EZ_replica$biplot_featcoord, biplot_EZ_pca_PCbiplot$Lmat)                              # successful replication of what seems an incorrect result

## -- Re-run replica, this time applying "fixes" to amend the source of discrepancies (see function annotations for details)
EZ_replica_fixes <- replica_biplotEZ(X, correlation.biplot = T, original = F)
isEqualUnname(EZ_replica_fixes$biplot_obscoord, my_biplot_observ_coord)                            
isEqualUnname(EZ_replica_fixes$biplot_featcoord, my_biplot_feat_coord)     


## (E) Geometrical properties tests
B_trsp_EZ <- t(EZ_replica$biplot_featcoord)
colnames(B_trsp_EZ) <- colnames(ext_PCA$biplot$B_trsp)
rownames(B_trsp_EZ) <- rownames(ext_PCA$biplot$B_trsp)

## (E.1) test 1: features' vector length (PC biplot only)
## -- caveat: biplot feat coord are multiplied by a "mystery" constant
pca_EZ_test_len_orig <- lenTest(B_trsp_EZ, 
                                Y = Y, 
                                #S = S_EZ                                      # covariance matrix provided internally
                                S = NULL
                                )                                              
pca_EZ_test_len_orig$testS
pca_EZ_test_len_orig$test_len                                                  # fail regardless of S
pca_EZ_test_len_orig$tab_out

## (E.2) test 2: correlation between features
pca_EZ_test_corcostest_orig <- corcosTest(B_trsp_EZ, Y)
pca_EZ_test_corcostest_orig$coscorr_test
pca_EZ_test_corcostest_orig$coscorr_export_tab

## (E.3) test 3: correlation between scores and features
pca_EZ_test_featPC <- corPCTest(V = V_EZ, Z = Z_EZ_est, Y = Y, lambda = lambda_EZ, 
                                S = S_EZ
                                #S = NULL
                                )
pca_EZ_test_featPC$check_eq19
pca_EZ_test_featPC$corr_feat_PC                                               # works if S = Y'Y still doesn't correspond to sd y                

## (E.4) test 4: euclid dist between observations on biplot = Mahalnaobis distance between observations in data matrix
pca_EZ_test_mahala_orig <- MahalaTest(A = EZ_replica$biplot_obscoord, Y = Y, S = S_EZ)          
pca_EZ_test_mahala_orig$mahalaTest  
pca_EZ_test_mahala_orig$dist_tab                                             # notice: this would work if the eucl dist sq was multiplied by the correct constant



## ---- E 3 [BENCHMARK] - BiplotGUI         ####
## probably a precursor of biplotEZ, same coauthors. Has only one function

test_bipGUI <- BiplotGUI::Biplots(X)

## View(BiplotGUI::Biplots)
## -- only exports via User Interface. 
## -- telling overal with biplot_EZ with some significant differences

## call my replica
myGUY <- replica_biplotGUI(X) 

## similarly to EZ biplot there are many discrepancies with our basic model
isEqualUnname(myGUY$Z_guy, my_Z)                # fail
isEqualUnname(myGUY$V_guy, my_V_eig)            # fail
isEqualAbs(myGUY$Z_guy, my_Z)                   # success in abs
isEqualAbs(myGUY$V_guy, my_V_eig)               # success in abs

isEqualUnname(myGUY$eigenval, my_lambda)        # fail
isEqualUnname(myGUY$sv_guy_wrong, my_sv)        # fail

## Eq.5 and 7
isEqualUnname(apply(myGUY$Z_guy,2,var), myGUY$eigenval)
isEqualUnname(sqrt(myGUY$eigenval), myGUY$sv_guy_wrong/sqrt(nrow(Y)-1) )

## PC Biplot coordinats
isEqualUnname(myGUY$biplot_obscoord_pc, my_biplot_observ_coord_eig)
isEqualUnname(myGUY$biplot_featcoord_pc, my_biplot_feat_coord_eig)



## Geometric properties test
B_trstp_guy <- t(myGUY$biplot_featcoord_pc) 
S_guy <- myGUY$S_guy

## -- Test 1: features vector length: fails (biplot feat coord are multiplied by a "mystery" constant)
pcaguy_test_len_orig <- lenTest(B_trstp_guy, Y = Y, S = S_guy)          # PCA mix doesn't provide a covariance matrix         
pcaguy_test_len_orig$testS
pcaguy_test_len_orig$test_len
pcaguy_test_len_orig$tab_out

## -- Test 2: correlation between features
pcaguy_test_corcostest_orig <- corcosTest(B_trstp_guy, Y)
pcaguy_test_corcostest_orig$coscorr_test
pcaguy_test_corcostest_orig$coscorr_export_tab

## -- Test 3: correlation between features and PCs
## CAVEAT: vegan does not provide scores as we know them
pcaguy_test_featPC <- corPCTest(V = myGUY$V_guy,
                              Z = myGUY$Z_guy,      
                              Y = Y, 
                              lambda = myGUY$eigenval,
                              S = S_guy
) 
pcaguy_test_featPC$testS
pcaguy_test_featPC$check_eq19
pcaguy_test_featPC$corr_feat_PC

## --- Test 4: Mahalanobis distance
## CAVEAT: the package doesn't provide a covariance matrix. We force the "correct" S
pcaguy_test_mahala_orig <- MahalaTest(A = myGUY$biplot_obscoord_pc, Y = Y, S = S_guy)          
pcaguy_test_mahala_orig$mahalaTest   
pcaguy_test_mahala_orig$dist_tab




## ---- E 23 Multibiplot                     ##
library(MultBiplotR)
test_mbr <- MultBiplotR::PCA.Biplot(X, 
                                    alpha = 0,                   # 0 = GH biplot i.e. PC biplot alternatively 1 = scores and loadings
                                    Scaling = 4)                 # specificed in function "InitialTransform". 4=  "Column centering":

isEqualUnname(test_mbr$EigenValues/(nrow(X)-1), my_lambda_svd)   # similar issue to EZ plot and others 
isEqualUnname(test_mbr$EV, my_V[,1:2])
isEqualUnname(test_mbr$RowCoordinates, my_biplot_observ_coord[,1:2])
isEqualUnname(test_mbr$ColCoordinates, t(my_biplot_feat_coord)[,1:2])
test_mbr$Structure





## ---- E 4 [BENCHMARK] - MultiBiplot       ####

pca_multi <- MultBiplotR::PCA.Analysis(X, 
                          dimension = ncol(X),
                          Scaling = 4)

## scores and loadings biplot
pca_multi_biplot_sl <- MultBiplotR::PCA.Biplot(X,
                                            dimension = ncol(X),
                                            Scaling = 4,
                                            alpha = 1)

## PC biplot
pca_multi_biplot_pcb <- MultBiplotR::PCA.Biplot(X,
                                               dimension = ncol(X),
                                               Scaling = 4,
                                               alpha = 0)

## check scale data is just Y
isEqualUnname(pca_multi$Scaled_Data, Y)

## Extract
lambda_multi <- pca_multi$EigenValues
V_multi <- pca_multi$EV
Z_multi_forced <- pca_multi_biplot_sl$RowCoordinates   # Scoring: no scoring matrix as such. Coorindates of sl biplot instead. Turns out this is wrong
multi_biplot_coord_observ_sl <- pca_multi_biplot_sl$RowCoordinates
multi_biplot_coord_feat_sl <- pca_multi_biplot_sl$ColCoordinates
multi_biplot_coord_observ_pcb <- pca_multi_biplot_pcb$RowCoordinates
multi_biplot_coord_feat_pcb <- pca_multi_biplot_pcb$ColCoordinates

## compare with from scratch
isEqualUnname(pca_multi$EV, my_V)                      # eigenvectors ok
isEqualUnname(pca_multi$EigenValues, my_lambda)        # issue with "eigenvalues" misnomer (= sv^2). Same for both functions
isEqualUnname(pca_multi$EigenValues, my_sv^2)
isEqualUnname(Z_multi_forced, my_Z_svd)
isEqualUnname(my_biplot_observ_coord, multi_biplot_coord_observ_pcb)
isEqualUnname(t(my_biplot_feat_coord), multi_biplot_coord_feat_pcb)
isEqualUnname(my_V, multi_biplot_coord_feat_sl)                     # not what you expect
isEqualUnname(my_Z, multi_biplot_coord_observ_sl) 

## Deep dive: run my replica
## -- Scores and biplot computaitons are correct to a point in the script, but further procesing causes coordinates to be inorrect even for the sl biplot
## -- see my replica functionb for details
my_multi_sl <- replica_multBiplotR(Y, alpha = 1)
my_multi_pcb <- replica_multBiplotR(Y, alpha = 0)

## check replication worked
isEqualUnname(my_multi_pcb$RowCoordinates, multi_biplot_coord_observ_pcb)
isEqualUnname(my_multi_pcb$ColCoordinates, multi_biplot_coord_feat_pcb)
isEqualUnname(my_multi_sl$RowCoordinates, multi_biplot_coord_observ_sl)
isEqualUnname(my_multi_sl$ColCoordinates, multi_biplot_coord_feat_sl)

## here is the fix: proving I isolated the issue
isEqualUnname(my_multi_pcb$correct_coord_observ, my_biplot_observ_coord) 
isEqualUnname(my_multi_pcb$correct_coord_feat, t(my_biplot_feat_coord))


## Testing Eq.5 & 7 with "as is" output
isEqualUnname(apply(Z_multi_forced,2,var), lambda_multi)
isEqualUnname(sqrt(lambda_multi), my_multi_sl$sv/sqrt(nrow(Y)-1))

## Gemoetrical properties
B_trsp_multi <- t(multi_biplot_coord_feat_pcb)

## -- Test 1: features vector length: fails (biplot feat coord are multiplied by a "mystery" constant)
multi_test_len_orig <- lenTest(B_trsp_multi, Y = Y)              
multi_test_len_orig$testS
multi_test_len_orig$test_len
multi_test_len_orig$tab_out

## -- Test 2: correlation between features
multi_test_corcostest_orig <- corcosTest(B_trsp_multi, Y)
multi_test_corcostest_orig$coscorr_test
multi_test_corcostest_orig$coscorr_export_tab

## -- Test 3: correlation between features and PCs
## CAVEAT: vegan does not provide scores as we know them
multi_test_featPC <- corPCTest(V = V_multi,
                                Z = Z_multi_forced,      
                                Y = Y, 
                                lambda = lambda_multi
) 
multi_test_featPC$testS
multi_test_featPC$check_eq19
multi_test_featPC$corr_feat_PC

pca_multi_biplot_pcb$Structure                            # the function returns correlations
isEqualUnname(c(pca_multi_biplot_pcb$Structure), c(multi_test_featPC$corr_feat_PC$corr))

## --- Test 4: Mahalanobis distance
## CAVEAT: the package doesn't provide a covariance matrix. We force the "correct" S
multi_test_mahala_orig <- MahalaTest(A = multi_biplot_coord_observ_pcb, Y = Y, S = cov(Y))          
multi_test_mahala_orig$mahalaTest   
multi_test_mahala_orig$dist_tab


##                                           ####
## Tables for paper (LaTeX export) ####
## basic matrices
clinic_Tab_X <- round(unitise_cols(X_real),2)
clinic_Tab_Y <- round(my_Y,2)
clinic_Tab_Zsvd <- round(my_Z_svd,2)
clinic_Tab_Z <- round(my_Z_svd,2)
clinic_Tab_A <- round(my_biplot_observ_coord,2)
clinic_Tab_B <- round(my_biplot_feat_coord,2)
clinic_Tab_A_r2 <- round(my_biplot_observ_coord,2)[,1:2]
clinic_Tab_B_r2 <- round(t(my_biplot_feat_coord),2)[,1:2]
clinic_Tab_V <- round(my_V,2)
clinic_Tab_lambda <- round(my_lambda,2)
clinic_Tab_sv <- round(my_sv, 2)

## properties
clinic_Tab_B_len <- round(ext_PCA$biplot$b_len_sd,2)
clinic_Tab_coscorr <- ext_PCA$biplot$coscorr[,c(1:2,4:7)]
clinic_Tab_coorPCFeat <- matrix(ext_corPC$corr_feat_PC$corr, nrow = ncol(V), byrow = T)
rownames(clinic_Tab_coorPCFeat) <- rownames(my_biplot_feat_coord)
colnames(clinic_Tab_coorPCFeat) <- colnames(my_biplot_feat_coord)

## Clinical example Table 1 
clinic_Tab_extd_obs_top <- cbind(clinic_Tab_X,
                         clinic_Tab_Y,
                         clinic_Tab_Zsvd,
                         clinic_Tab_A_r2)
clinic_Tab_extd_obs_col_means <-  apply(clinic_Tab_extd_obs,2, mean)                                                          # mean of raw data
clinic_Tab_extd_obs_col_s.var  <- apply(clinic_Tab_extd_obs,2,var)                                                            # sample variance of raw data, see section 2.1 in the paper
clinic_Tab_extd_obs_top_bottom <- rbind(clinic_Tab_extd_obs_col_means, clinic_Tab_extd_obs_col_s.var)
rownames(clinic_Tab_extd_obs_top_bottom) <- c("mean", "s.var")
clinic_Tab_extd_obs <- rbind(clinic_Tab_extd_obs_top, clinic_Tab_extd_obs_top_bottom)

## Clinical example Table 2
clinic_Tab_extd_feat <- cbind(clinic_Tab_V,
                                  clinic_Tab_B_r2,
                                  clinic_Tab_B_len)

clinic_Tab_extd_svev <- rbind(clinic_Tab_lambda,
                                     clinic_Tab_sv)





library(xtable)
print(xtable(clinic_Tab_extd_obs, type = "latex"), file = "export_clinical_Tab_obs.tex")
print(xtable(clinic_Tab_extd_feat, type = "latex"), file = "export_clinical_Tab_feat.tex")
print(xtable(clinic_Tab_extd_svev, type = "latex"), file = "export_clinical_Tab_feat_svev.tex")
print(xtable(clinic_Tab_coscorr, type = "latex"), file = "export_clinical_Tab_feat_corr.tex")
print(xtable(t(clinic_Tab_coorPCFeat), type = "latex"), file = "export_clinical_Tab_featPCcorr.tex")

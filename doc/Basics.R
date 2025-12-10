## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 6,comment = "#",
                      collapse = TRUE,results = "hold",
					  fig.align = "center")
library(fashr)
library(ggplot2)

## ----fig.height=8, fig.width=7------------------------------------------------
set.seed(1)
N <- 100
propA <- 0.5; propB <- 0.3; propC <- 0.2
sigma_vec <- c(0.05, 0.1, 0.2)

sizeA <- N * propA
data_sim_list_A <- lapply(1:sizeA, function(i) simulate_process(sd_poly = 0.2, type = "nondynamic", sd = sigma_vec, normalize = TRUE))

sizeB <- N * propB
if(sizeB > 0){
data_sim_list_B <- lapply(1:sizeB, function(i) simulate_process(sd_poly = 1, type = "linear", sd = sigma_vec, normalize = TRUE))

}else{
  data_sim_list_B <- list()
}

sizeC <- N * propC
data_sim_list_C <- lapply(1:sizeC, function(i) simulate_process(sd_poly = 0, type = "nonlinear", sd = sigma_vec, sd_fun = 1, p = 1, normalize = TRUE))

datasets <- c(data_sim_list_A, data_sim_list_B, data_sim_list_C)
labels <- c(rep("A", sizeA), rep("B", sizeB), rep("C", sizeC))
indices_A <- 1:sizeA
indices_B <- (sizeA + 1):(sizeA + sizeB)
indices_C <- (sizeA + sizeB + 1):(sizeA + sizeB + sizeC)

dataset_labels <- rep(as.character(NA),100)
dataset_labels[indices_A] <- paste0("A",seq(1,length(indices_A)))
dataset_labels[indices_B] <- paste0("B",seq(1,length(indices_B)))
dataset_labels[indices_C] <- paste0("C",seq(1,length(indices_C)))
names(datasets) <- dataset_labels

par(mfrow = c(3, 3))
for(i in indices_A[1:3]){
  plot(datasets[[i]]$x, datasets[[i]]$y, type = "p", col = "black", lwd = 2, xlab = "time", ylab = "effect", ylim = c(-1.5, 1.5), main = paste("Category A: ", i))
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lwd = 1)
}

for(i in indices_B[1:3]){
  plot(datasets[[i]]$x, datasets[[i]]$y, type = "p", col = "black", lwd = 2, xlab = "time", ylab = "effect", ylim = c(-1.5, 1.5), main = paste("Category B: ", i))
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lwd = 1)
}

for(i in indices_C[1:3]){
  plot(datasets[[i]]$x, datasets[[i]]$y, type = "p", col = "black", lwd = 2, xlab = "time", ylab = "effect", ylim = c(-1.5, 1.5), main = paste("Category C: ", i))
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lwd = 1)
}

## -----------------------------------------------------------------------------
length(datasets)

## -----------------------------------------------------------------------------
datasets[[1]]

## -----------------------------------------------------------------------------
table(labels)

## ----results="hide"-----------------------------------------------------------
fash_fit1 <- fash(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, order = 1)
fash_fit1

## -----------------------------------------------------------------------------
fash_fit1

## -----------------------------------------------------------------------------
fash_fit1$prior_weights

## ----fig.height=3.5, fig.width=4----------------------------------------------
psd_colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00", "#ffff33", "#a65628")
visualize_fash_prior(fash_fit1, constraints = "initial") +
  scale_color_manual(values = psd_colors)
visualize_fash_prior(fash_fit1, constraints = "orthogonal") +
  scale_color_manual(values = psd_colors)

## -----------------------------------------------------------------------------
fash_fit1_adj <- BF_update(fash_fit1)
fash_fit1_adj$prior_weights

## ----fig.height=3.5, fig.width=4----------------------------------------------
visualize_fash_prior(fash_fit1_adj, constraints = "initial") +
  scale_color_manual(values = psd_colors)
visualize_fash_prior(fash_fit1_adj, constraints = "orthogonal") +
  scale_color_manual(values = psd_colors)

## ----fig.height=4.5, fig.width=3.5--------------------------------------------
plot(fash_fit1_adj, plot_type = "heatmap",
     selected_indices = c(paste0("A",1:5),
	                      paste0("B",1:5),
						  paste0("C",1:5)))

## ----fig.height=4, fig.width=5------------------------------------------------
plot(fash_fit1_adj, plot_type = "structure", discrete = TRUE)

## ----fig.height=4, fig.width=4------------------------------------------------
fdr_result1_adj <- fdr_control(fash_fit1_adj, alpha = 0.1, plot = TRUE)

## -----------------------------------------------------------------------------
detected_indices1 <- which(fash_fit1_adj$lfdr < 0.01)
length(detected_indices1)

## -----------------------------------------------------------------------------
sum(labels[detected_indices1] != "A")/(sizeC + sizeB)

## -----------------------------------------------------------------------------
sum(labels[detected_indices1] == "A")/length(detected_indices1)

## -----------------------------------------------------------------------------
fitted_beta <- predict(fash_fit1_adj, index = detected_indices1[1])
fitted_beta

## -----------------------------------------------------------------------------
fitted_beta_new <- predict(fash_fit1_adj, index = detected_indices1[1],
                           smooth_var = seq(0, 16, length.out = 100))
head(fitted_beta_new)

## -----------------------------------------------------------------------------
fitted_beta_samples <- predict(fash_fit1_adj, index = detected_indices1[1], 
                               smooth_var = seq(0, 16, length.out = 100), 
                               only.samples = TRUE, M = 50)
str(fitted_beta_samples)

## -----------------------------------------------------------------------------
plot(fitted_beta_new$x,fitted_beta_new$mean,type = "l",lwd = 2,
     xlab = "condition",ylab = "function value")

## ----fig.height=4, fig.width=4------------------------------------------------
plot(fash_fit1_adj, selected_unit = detected_indices1[1],
     plot_type = "function")
lines(datasets[[detected_indices1[1]]]$x, datasets[[detected_indices1[1]]]$truef,
      col = "black", lwd = 1, lty = "dashed")

## ----results="hide"-----------------------------------------------------------
evaluations <- seq(1, 16, length.out = 100)
functional_F <- function(x){
  max_first_half <- max(x[evaluations <= 8])
  max_second_half <- max(x[evaluations > 8])
  return(max_first_half - max_second_half)
}
lfsr_result <- testing_functional(
  functional = functional_F,
  fash = fash_fit1_adj,
  indices = 1:100,
  smooth_var = evaluations
)

## ----fig.height=2.5, fig.width=7----------------------------------------------
par(mfrow = c(1,3))
lfsr_result <- lfsr_result[order(lfsr_result$indices),]
hist(lfsr_result[indices_A,"lfsr"],n = 16,main = "A",xlab = "lfsr")
hist(lfsr_result[indices_B,"lfsr"],n = 16,main = "B",xlab = "lfsr")
hist(lfsr_result[indices_C,"lfsr"],n = 16,main = "C",xlab = "lfsr")

## ----results="hide"-----------------------------------------------------------
functional_switch <- function(x){
  x_pos <- x[x > 0]
  x_neg <- x[x < 0]
  if(length(x_pos) == 0 || length(x_neg) == 0){
    return(0)
  }
  min(max(abs(x_pos)), max(abs(x_neg))) - 0.25
}
lfsr_switch <- testing_functional(
  functional = functional_switch,
  fash = fash_fit1_adj,
  indices = 1:100,
  lfsr_cal = function(x) mean(x <= 0),
  smooth_var = evaluations
)

## ----fig.height=2.5, fig.width=7----------------------------------------------
par(mfrow = c(1,3))
lfsr_switch <- lfsr_switch[order(lfsr_switch$indices),]
hist(lfsr_switch[indices_A,"lfsr"],n = 16,main = "A",xlab = "lfsr")
hist(lfsr_switch[indices_B,"lfsr"],n = 16,main = "B",xlab = "lfsr")
hist(lfsr_switch[indices_C,"lfsr"],n = 16,main = "C",xlab = "lfsr")

## ----fig.height=4, fig.width=4------------------------------------------------
par(mfrow = c(1,1))
detected_switch <- lfsr_switch$indices[lfsr_switch$lfsr < 0.01]
plot(fash_fit1_adj, selected_unit = detected_switch[1],
     plot_type = "function")

## ----results="hide"-----------------------------------------------------------
fash_fit2 <- fash(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, order = 2)

## -----------------------------------------------------------------------------
fash_fit2_adj <- BF_update(fash_fit2)

## ----fig.height=3.5, fig.width=4----------------------------------------------
visualize_fash_prior(fash_fit2, constraints = "initial")
visualize_fash_prior(fash_fit2_adj, constraints = "initial")

## -----------------------------------------------------------------------------
detected_indices2 <- which(fash_fit2_adj$lfdr < 0.01)

## -----------------------------------------------------------------------------
sum(labels[detected_indices2] == "C")/sizeC

## -----------------------------------------------------------------------------
sum(labels[detected_indices2] != "C")/length(detected_indices2)

## ----fig.height=4, fig.width=4------------------------------------------------
selected_index <- detected_indices2[1]
plot(fash_fit2_adj, selected_unit = selected_index, plot_type = "function")
lines(datasets[[selected_index]]$x, datasets[[selected_index]]$truef,
      col = "black", lwd = 1, lty = "dashed")


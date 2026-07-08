pkgload::load_all(quiet = TRUE)
n_envopt <- 10000L
# ... minimal setup same as compare script ...
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group1 <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight1 <- c(ctl, trt)
n1 <- length(weight1)
weight <- c(weight1, weight1)
p_setup <- Prior_Setup(weight ~ group, data = data.frame(weight = weight, group = group_stacked <- factor(c(group1, group1), levels = levels(group1))), family = gaussian())
y <- p_setup$y; x_block <- p_setup$x; l1 <- ncol(x_block)
block <- factor(c(rep("B1", n1), rep("B2", n1)), levels = c("B1", "B2"))
x1 <- x_block[seq_len(n1), , drop = FALSE]; zeros <- matrix(0, n1, l1)
x_old <- rbind(cbind(x1, zeros), cbind(zeros, x1))
mu1 <- p_setup$mu; Sigma1 <- p_setup$Sigma
prior_list_old <- list(mu = c(mu1, mu1), Sigma = rbind(cbind(Sigma1, matrix(0,l1,l1)), cbind(matrix(0,l1,l1), Sigma1)), dispersion = p_setup$dispersion, shape = p_setup$shape, rate = p_setup$rate, Precision = solve(rbind(cbind(Sigma1, matrix(0,l1,l1)), cbind(matrix(0,l1,l1), Sigma1))), max_disp_perc = 0.99)
prior_list_block <- list(mu = mu1, Sigma = Sigma1, shape = p_setup$shape, rate = p_setup$rate, max_disp_perc = 0.99)
sim_env <- rindepNormalGamma_reg_with_envelope(n=1L,y=y,x=x_old,prior_list=prior_list_old,n_envopt=n_envopt,Gridtype=3L,use_parallel=FALSE,progbar=FALSE)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(n=1L,y=y,x=x_block,block=block,prior_list=prior_list_block,n_envopt=n_envopt,Gridtype=3L,use_parallel=FALSE,progbar=FALSE,verbose=FALSE,offset=rep(0,length(y)),wt=rep(1,length(y)),p_re=-1L,n_rss_iter=10L,RSS_ML=NA_real_,use_opencl=FALSE,group_levels=character(0),re_names=character(0))
cbL <- sim_env$Envelope$cbars; cbB1 <- sim_block$build_out$block_envelopes[[1]]$cbars; cbB2 <- sim_block$build_out$block_envelopes[[2]]$cbars
lgL <- sim_env$UB_list$lg_prob_factor; lgB1 <- sim_block$build_out$dispersion_envelope$block_dispersion[[1]]$lg_prob_factor
cat('cbL dim', paste(dim(cbL), collapse='x'), ' cbB1', paste(dim(cbB1), collapse='x'), '\n')
cat('Sample cbars legacy row 1:', paste(round(cbL[1,],4), collapse=' '), '\n')
cat('Sample cbars B1 row 1:     ', paste(round(cbB1[1,],4), collapse=' '), '\n')
# nearest legacy face to each B1 cbars in cols 1:2
for (j in 1:3) {
  d <- apply(cbL, 1, function(r) max(abs(r[1:2]-cbB1[j,])))
  cat('B1 face', j, 'best legacy match max|diff| cols1:2 =', min(d), 'at legacy face', which.min(d), '\n')
}
# compare lg when pairing by GridIndex if available
giL <- sim_env$Envelope$GridIndex
giB <- sim_block$build_out$block_envelopes[[1]]$GridIndex
cat('GridIndex legacy length', length(giL), ' block', length(giB), '\n')
if (!is.null(giL) && !is.null(giB)) cat('First 5 legacy GI:', paste(head(giL,5), collapse=' '), '\n')

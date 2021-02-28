# run 
options(warn=-1)
# simulation of the data
source('./model/lib.R')
source('./model/param.R')
source('./model/simu.R')

eff_matrix_mtx <- se_matrix_mtx <- p_matrix_mtx <- eff_matrix_det_mtx <- se_matrix_det_mtx <- p_matrix_det_mtx <- eff_matrix_sto_mtx  <- se_matrix_sto_mtx <- p_matrix_sto_mtx <- NULL

ci_blower_det <- ci_bupper_det <- ci_blower_sto <- ci_bupper_sto <- NULL
ci_blower_det_d <- ci_bupper_det_d <- ci_blower_sto_d <- ci_bupper_sto_d <- NULL
ci_asy_lower_det <- ci_asy_upper_det <- ci_asy_lower_sto <- ci_asy_upper_sto <- NULL
length_list <- within_list<- list()

length_q_b <- length_es_b <- length_q_db <- length_es_db <- length_q <- length_es <- NULL
within_q_b <- within_es_b <- within_q_db <- within_es_db <- within_q <- within_es <- NULL
set.seed(12345)
for(s in 1:nsim){
    o_1 <- A_1 <- o_2 <- A_2 <- y_1 <- y_2 <- epsilon <- rep(NA, n)
###############################################################
for(i in 1:n){
    y_1[i] <- 0
    A_1[i] <- 2 * rbinom(1, 1, p) - 1
    o_1[i] <- 2 * rbinom(1, 1, p) - 1
    A_2[i] <- 2 * rbinom(1, 1, p) - 1
    p_v <- expit(o_1[i], A_1[i], delta_1, delta_2)
    o_2[i] <- 2 * rbinom(1, 1, p_v) - 1
    epsilon[i] <- rnorm(1, 0, 1)
    y_2[i] <- gamma_1 + gamma_2 * o_1[i] + gamma_3 * A_1[i] + gamma_4 * o_1[i] * A_1[i] +
    gamma_5 * A_2[i] + gamma_6 * o_2[i] * A_2[i] + gamma_7 * A_1[i] * A_2[i] + epsilon[i]
}
y_2 <- round(y_2, 4)
usedata <- as.data.frame(cbind(A_1, o_1, A_2, o_2, y_1, y_2))
colnames(usedata) <- c("A1", "O1", "A2", "O2", "Y1", "Y2")

qfit <- lm(Y2 ~ O1 + A1 + O1:A1 + A2 + O2:A2 + A1:A2, data = usedata)

# inference at the second stage
eff_matrix <-  as.matrix(round(summary(qfit)$coef, 3 )[,1])
se_matrix <-  as.matrix(round(summary(qfit)$coef, 3 )[,2])
p_matrix <-  as.matrix(round(summary(qfit)$coef, 3 )[,4])



pre_data <- cbind(rep(1, dim(usedata)[1]), usedata$O1, usedata$A1, rep(1, length(usedata$A2)), usedata$O1 * usedata$A1, usedata$O2, usedata$A1)

eff_A1 <- as.matrix(eff_matrix[c(1:3, 5), ])
eff_A2 <- as.matrix(eff_matrix[c(4, 6:7), ])

pre_A1 <- pre_data[, c(1:3, 5)]
pre_A2 <- pre_data[, c(4, 6:7)]

first_part <- pre_A1 %*% eff_A1
second_part <- pre_A2 %*% eff_A2
mean_return <- cbind(first_part - 1 * second_part, first_part + 1 * second_part)

prob_matrix <- softmax(c(-1, 1), second_part) 
#  Pseudo outcomes
q_det_pseudo <- c(first_part + abs(second_part))
q_sto_pseudo <- apply(mean_return * prob_matrix, 1, sum)
usedata$q_det_pseudo <- q_det_pseudo
usedata$q_sto_pseudo <- q_sto_pseudo


# first stage
qfit_det <- lm(q_det_pseudo ~ O1 + A1 + O1:A1, data = usedata)
qfit_sto <- lm(q_sto_pseudo ~ O1 + A1 + O1:A1, data = usedata)

qfit_det_re <- as.matrix(round(summary(qfit_det)$coef, 3)[,1])
qfit_sto_re <- as.matrix(round(summary(qfit_sto)$coef, 3)[,1])

# bootstraps
qfit_det_b1_re_mtx <- qfit_sto_b1_re_mtx <- NULL
u_list_1p_mtx <- u_list_2p_mtx <- NULL
for(b1 in 1:B1){
    index <- sample(1:dim(usedata)[1], nb1, replace = FALSE)
    usedata_b1 <- usedata[index, ]
    qfit_det_b1 <- lm(q_det_pseudo ~ O1 + A1 + O1:A1, data = usedata_b1)
    qfit_sto_b1 <- lm(q_sto_pseudo ~ O1 + A1 + O1:A1, data = usedata_b1)

    qfit_det_b1_re <- as.matrix(round(summary(qfit_det_b1)$coef, 3)[,1])
    qfit_sto_b1_re <- as.matrix(round(summary(qfit_sto_b1)$coef, 3)[,1])

    qfit_det_b1_re_mtx <- cbind(qfit_det_b1_re_mtx, qfit_det_b1_re)
    qfit_sto_b1_re_mtx <- cbind(qfit_sto_b1_re_mtx, qfit_sto_b1_re)
    u_list_1 <- u_list_2 <- NULL
    for(b2 in 1:B2){
        index_b2 <- sample(index, nb2, replace = FALSE)
        usedata_b2 <- usedata[index_b2, ]

        qfit_det_b2 <- lm(q_det_pseudo ~ O1 + A1 + O1:A1, data = usedata_b2)
        qfit_sto_b2 <- lm(q_sto_pseudo ~ O1 + A1 + O1:A1, data = usedata_b2)

        qfit_det_b2_re <- as.matrix(round(summary(qfit_det_b2)$coef, 3)[,1])
        qfit_sto_b2_re <- as.matrix(round(summary(qfit_sto_b2)$coef, 3)[,1])
         
        u_list_part_1 <- (c(qfit_det_b2_re - qfit_det_re) <= 0) * 1
        u_list_part_2 <- (c(qfit_sto_b2_re - qfit_sto_re) <= 0) * 1

        u_list_1 <- cbind(u_list_1, u_list_part_1)
        u_list_2 <- cbind(u_list_2, u_list_part_2)
    }
    u_list_1p <- apply(u_list_1, 1, mean)
    u_list_2p <- apply(u_list_2, 1, mean)
    u_list_1p_mtx <- cbind(u_list_1p_mtx, u_list_1p)
    u_list_2p_mtx <- cbind(u_list_2p_mtx, u_list_2p)
}
qlist1l <- as.vector(apply(u_list_1p_mtx, 1, function(x) quantile(x, a1)))
qlist1u <- as.vector(apply(u_list_1p_mtx, 1, function(x) quantile(x, a2)))

# without interpolation    
qlist1l <- sapply(1:length(qlist1l), function(m) min(u_list_1p_mtx[m, u_list_1p_mtx[m, ] >= qlist1l[m]]))
qlist1u <- sapply(1:length(qlist1u), function(m) min(u_list_1p_mtx[m, u_list_1p_mtx[m, ] >= qlist1u[m]]))


qlist2l <- as.vector(apply(u_list_2p_mtx, 1, function(x) quantile(x, a1)))
qlist2u <- as.vector(apply(u_list_2p_mtx, 1, function(x) quantile(x, a2)))

# without interpolation
qlist2l <- sapply(1:length(qlist2l), function(m) min(u_list_2p_mtx[m, u_list_2p_mtx[m, ] >= qlist2l[m]]))
qlist2u <- sapply(1:length(qlist2u), function(m) min(u_list_2p_mtx[m, u_list_2p_mtx[m, ] >= qlist2u[m]]))


## Bootstrap
# deterministically 
blower_det <- as.vector(apply(qfit_det_b1_re_mtx, 1, function(x) quantile(x, a1)))
bupper_det <- as.vector(apply(qfit_det_b1_re_mtx, 1, function(x) quantile(x, a2)))

blower_det <- sapply(1:length(blower_det), function(m) min(qfit_det_b1_re_mtx[m, qfit_det_b1_re_mtx[m, ] >= blower_det[m]]))
bupper_det <- sapply(1:length(bupper_det), function(m) min(qfit_det_b1_re_mtx[m, qfit_det_b1_re_mtx[m, ] >= bupper_det[m]]))

# stochastically 
blower_sto <- as.vector(apply(qfit_sto_b1_re_mtx, 1, function(x) quantile(x, a1)))
bupper_sto <- as.vector(apply(qfit_sto_b1_re_mtx, 1, function(x) quantile(x, a2)))

blower_sto <- sapply(1:length(blower_sto), function(m) min(qfit_sto_b1_re_mtx[m, qfit_sto_b1_re_mtx[m, ] >= blower_sto[m]]))
bupper_sto <- sapply(1:length(bupper_sto), function(m) min(qfit_sto_b1_re_mtx[m, qfit_sto_b1_re_mtx[m, ] >= bupper_sto[m]]))



## Double bootstrap
# deterministically 
blower_det_d <- as.vector(sapply(1:dim(qfit_det_b1_re_mtx)[1], function(m) quantile(qfit_det_b1_re_mtx[m, ], qlist1l[m])))
bupper_det_d <- as.vector(sapply(1:dim(qfit_det_b1_re_mtx)[1], function(m) quantile(qfit_det_b1_re_mtx[m, ], qlist1u[m])))

blower_det_d <- sapply(1:length(blower_det_d), function(m) min(qfit_det_b1_re_mtx[m, qfit_det_b1_re_mtx[m, ] >= blower_det_d[m]]))
bupper_det_d <- sapply(1:length(bupper_det_d), function(m) min(qfit_det_b1_re_mtx[m, qfit_det_b1_re_mtx[m, ] >= bupper_det_d[m]]))

# stochastically 
blower_sto_d <- as.vector(sapply(1:dim(qfit_sto_b1_re_mtx)[1], function(m) quantile(qfit_sto_b1_re_mtx[m, ], qlist2l[m])))
bupper_sto_d <- as.vector(sapply(1:dim(qfit_sto_b1_re_mtx)[1], function(m) quantile(qfit_sto_b1_re_mtx[m, ], qlist2u[m])))

blower_sto_d <- sapply(1:length(blower_sto_d), function(m) min(qfit_sto_b1_re_mtx[m, qfit_sto_b1_re_mtx[m, ] >= blower_sto_d[m]]))
bupper_sto_d <- sapply(1:length(bupper_sto_d), function(m) min(qfit_sto_b1_re_mtx[m, qfit_sto_b1_re_mtx[m, ] >= bupper_sto_d[m]]))



eff_matrix_det <-  as.matrix(round(summary(qfit_det)$coef, 3 )[,1])
se_matrix_det <-  as.matrix(round(summary(qfit_det)$coef, 3 )[,2])
p_matrix_det <-  as.matrix(round(summary(qfit_det)$coef, 3 )[,4])

eff_matrix_sto <-  as.matrix(round(summary(qfit_sto)$coef, 3 )[,1])
se_matrix_sto <-  as.matrix(round(summary(qfit_sto)$coef, 3 )[,2])
p_matrix_sto <-  as.matrix(round(summary(qfit_sto)$coef, 3 )[,4])
    # eff_matrix_mtx <- cbind(eff_matrix_mtx, result$eff_matrix)
    # eff_matrix_det_mtx <- cbind(eff_matrix_det_mtx, result$eff_matrix_det)
    # eff_matrix_sto_mtx <- cbind(eff_matrix_sto_mtx, result$eff_matrix_sto)

    # se_matrix_mtx <- cbind(se_matrix_mtx, result$se_matrix)
    # se_matrix_det_mtx <- cbind(se_matrix_det_mtx, result$se_matrix_det)
    # se_matrix_sto_mtx <- cbind(se_matrix_sto_mtx, result$se_matrix_sto)

    # p_matrix_mtx <- cbind(p_matrix_mtx, result$p_matrix)
    # p_matrix_det_mtx <- cbind(p_matrix_det_mtx, result$p_matrix_det)
    # p_matrix_sto_mtx <- cbind(p_matrix_sto_mtx, result$p_matrix_sto)

    # ci_blower_det <- cbind(ci_blower_det, result$blower_det)
    # ci_bupper_det <- cbind(ci_bupper_det, result$bupper_det)
    # ci_blower_sto <- cbind(ci_blower_sto, result$blower_sto)
    # ci_bupper_sto <- cbind(ci_bupper_sto, result$bupper_sto)

    # ci_blower_det_d <- cbind(ci_blower_det_d, result$blower_det_d)
    # ci_bupper_det_d <- cbind(ci_bupper_det_d, result$bupper_det_d)
    # ci_blower_sto_d <- cbind(ci_blower_sto_d, result$blower_sto_d)
    # ci_bupper_sto_d <- cbind(ci_bupper_sto_d, result$bupper_sto_d)  

    # length_bootstrap <- cbind(result$bupper_det - result$blower_det, result$bupper_sto - result$blower_sto, result$bupper_det_d - result$blower_det_d, result$bupper_sto_d - result$blower_sto_d)
    # length_all <- cbind(length_bootstrap, result$se_matrix_det * 1.96 * 2, result$se_matrix_sto * 1.96 * 2)

    # length_q_b <- cbind(length_q_b, result$bupper_det - result$blower_det)
    # length_es_b <- cbind(length_es_b, result$bupper_sto - result$blower_sto)
    # length_q_db <- cbind(length_q_db, result$bupper_det_d - result$blower_det_d)
    # length_es_db <- cbind(length_es_db, result$bupper_sto_d - result$blower_sto_d)
    
    # within_bootstrap <- cbind((result$bupper_det >= result$eff_matrix_det) * (result$blower_det <= result$eff_matrix_det) , 
    #                       (result$bupper_sto >= result$eff_matrix_sto) * (result$blower_sto <= result$eff_matrix_sto) , 
    #                       (result$bupper_det_d >= result$eff_matrix_det) * (result$blower_det_d <= result$eff_matrix_det) , 
    #                       (result$bupper_sto_d >= result$eff_matrix_sto) * (result$blower_sto_d <= result$eff_matrix_sto))
    
    # within_q_b <- cbind(within_q_b, (result$bupper_det >= result$eff_matrix_det) * (result$blower_det <= result$eff_matrix_det))
    # within_es_b <- cbind(within_es_b, (result$bupper_sto >= result$eff_matrix_sto) * (result$blower_sto <= result$eff_matrix_sto))
    # within_q_db <- cbind(within_q_db, (result$bupper_det_d >= result$eff_matrix_det) * (result$blower_det_d <= result$eff_matrix_det))
    # within_es_db <- cbind(within_es_db, (result$bupper_sto_d >= result$eff_matrix_sto) * (result$blower_sto_d <= result$eff_matrix_sto))

    # colnames(length_all) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap", "Q", "ESARSA")
    # colnames(within_bootstrap) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap")
    # length_list[[i]] <- length_all
    # within_list[[i]] <- within_bootstrap 
    print(s)
}
# second-stage
emp_se <- round( apply(eff_matrix_mtx, 1, sd), 4)
mean_se <- round(apply(se_matrix_mtx, 1, mean), 4)
p_dist <- round(rowMeans((p_matrix_mtx < 0.05)*1), 4)

# first-stage deteriministically 
emp_se_det <- round( apply(eff_matrix_det_mtx, 1, sd), 4)
mean_se_det <- round(apply(se_matrix_det_mtx, 1, mean), 4)
p_dist_det <- round(rowMeans((p_matrix_det_mtx < 0.05)*1), 4)

# first-stage stochastically
emp_se_sto <- round( apply(eff_matrix_sto_mtx, 1, sd), 4)
mean_se_sto <- round(apply(se_matrix_sto_mtx, 1, mean), 4)
p_dist_sto <- round(rowMeans((p_matrix_sto_mtx < 0.05)*1), 4)

sec_result <- rbind(emp_se, mean_se, p_dist)
first_result <- rbind(emp_se_det, mean_se_det, p_dist_det, emp_se_sto, mean_se_sto, p_dist_sto)


# length
avg_q_b <- apply(length_q_b, 1, mean)
avg_es_b <- apply(length_es_b, 1, mean)
avg_q_db <- apply(length_q_db, 1, mean)
avg_es_db <- apply(length_es_db, 1, mean)
length_output <- cbind(avg_q_b, avg_es_b, avg_q_db, avg_es_db, 2 * 1.96 * mean_se_det, 2 * 1.96 * mean_se_sto)
colnames(length_output) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap", "Q", "ESARSA")


# within
within_output <- cbind(apply(within_q_b, 1, mean), 
                       apply(within_es_b, 1, mean), 
                       apply(within_q_db, 1, mean), 
                       apply(within_es_db, 1, mean))
colnames(within_output) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap")                     

####### save results ##################################################
save(length_list, file = "./output/length_list.RData")
save(within_list, file = "./output/within_list.RData")
write.csv(length_output, file = "./output/length_output.csv")
write.csv(within_output, file = "./output/within_output.csv")
#######################################################################
save(eff_matrix_det_mtx, file = "./output/eff_matrix_det_mtx.RData")
save(eff_matrix_sto_mtx, file = "./output/eff_matrix_sto_mtx.RData")
save(se_matrix_det_mtx , file = "./output/se_matrix_det_mtx.RData")
save(se_matrix_sto_mtx, file = "./output/se_matrix_sto_mtx.RData")
save(ci_blower_det, file = "./output/ci_blower_det.RData")
save(ci_bupper_det, file = "./output/i_bupper_det.RData")
save(ci_blower_sto , file = "./output/ci_blower_sto.RData")
save(ci_bupper_sto , file = "./output/ci_bupper_sto.RData")

save(ci_blower_det_d , file = "./output/ci_blower_det_d.RData")
save(ci_bupper_det_d, file = "./output/ci_bupper_det_d.RData")
save(ci_blower_sto_d , file = "./output/ci_blower_sto_d.RData")
save(ci_bupper_sto_d , file = "./output/ci_bupper_sto_d.RData")

write.csv(sec_result, file = "./output/sec_result.csv")
write.csv(first_result, file = "./output/first_result.csv")



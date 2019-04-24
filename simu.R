# simulation of the data
source('lib.R')
source('param.R')



simu <- function(){
# simulation study simulate the data
###############################################################
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


return(list(eff_matrix = eff_matrix, se_matrix = se_matrix, p_matrix = p_matrix, eff_matrix_det = eff_matrix_det, 
            se_matrix_det= se_matrix_det, p_matrix_det= p_matrix_det, eff_matrix_sto = eff_matrix_sto, se_matrix_sto = se_matrix_sto, 
            p_matrix_sto= p_matrix_sto, blower_det = blower_det, bupper_det = bupper_det, blower_sto = blower_sto, bupper_sto = bupper_sto, 
            blower_det_d = blower_det_d, bupper_det_d = bupper_det_d, blower_sto_d = blower_sto_d, bupper_sto_d = bupper_sto_d))
}

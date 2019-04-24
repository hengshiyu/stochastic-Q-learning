# run 
options(warn=-1)
# simulation of the data
source('lib.R')
source('param.R')
source('simu.R')

eff_matrix_mtx <- se_matrix_mtx <- p_matrix_mtx <- eff_matrix_det_mtx <- se_matrix_det_mtx <- p_matrix_det_mtx <- eff_matrix_sto_mtx  <- se_matrix_sto_mtx <- p_matrix_sto_mtx <- NULL

ci_blower_det <- ci_bupper_det <- ci_blower_sto <- ci_bupper_sto <- NULL
ci_blower_det_d <- ci_bupper_det_d <- ci_blower_sto_d <- ci_bupper_sto_d <- NULL
ci_asy_lower_det <- ci_asy_upper_det <- ci_asy_lower_sto <- ci_asy_upper_sto <- NULL
length_list <- within_list<- list()

length_q_b <- length_es_b <- length_q_db <- length_es_db <- length_q <- length_es <- NULL
within_q_b <- within_es_b <- within_q_db <- within_es_db <- within_q <- within_es <- NULL
set.seed(12345)
for(s in 1:nsim){
    result <- simu()
    eff_matrix_mtx <- cbind(eff_matrix_mtx, result$eff_matrix)
    eff_matrix_det_mtx <- cbind(eff_matrix_det_mtx, result$eff_matrix_det)
    eff_matrix_sto_mtx <- cbind(eff_matrix_sto_mtx, result$eff_matrix_sto)

    se_matrix_mtx <- cbind(se_matrix_mtx, result$se_matrix)
    se_matrix_det_mtx <- cbind(se_matrix_det_mtx, result$se_matrix_det)
    se_matrix_sto_mtx <- cbind(se_matrix_sto_mtx, result$se_matrix_sto)

    p_matrix_mtx <- cbind(p_matrix_mtx, result$p_matrix)
    p_matrix_det_mtx <- cbind(p_matrix_det_mtx, result$p_matrix_det)
    p_matrix_sto_mtx <- cbind(p_matrix_sto_mtx, result$p_matrix_sto)

    ci_blower_det <- cbind(ci_blower_det, result$blower_det)
    ci_bupper_det <- cbind(ci_bupper_det, result$bupper_det)
    ci_blower_sto <- cbind(ci_blower_sto, result$blower_sto)
    ci_bupper_sto <- cbind(ci_bupper_sto, result$bupper_sto)

    ci_blower_det_d <- cbind(ci_blower_det_d, result$blower_det_d)
    ci_bupper_det_d <- cbind(ci_bupper_det_d, result$bupper_det_d)
    ci_blower_sto_d <- cbind(ci_blower_sto_d, result$blower_sto_d)
    ci_bupper_sto_d <- cbind(ci_bupper_sto_d, result$bupper_sto_d)  

    length_bootstrap <- cbind(result$bupper_det - result$blower_det, result$bupper_sto - result$blower_sto, result$bupper_det_d - result$blower_det_d, result$bupper_sto_d - result$blower_sto_d)
    length_all <- cbind(length_bootstrap, result$se_matrix_det * 1.96 * 2, result$se_matrix_sto * 1.96 * 2)

    length_q_b <- cbind(length_q_b, result$bupper_det - result$blower_det)
    length_es_b <- cbind(length_es_b, result$bupper_sto - result$blower_sto)
    length_q_db <- cbind(length_q_db, result$bupper_det_d - result$blower_det_d)
    length_es_db <- cbind(length_es_db, result$bupper_sto_d - result$blower_sto_d)
    
    within_bootstrap <- cbind((result$bupper_det >= result$eff_matrix_det) * (result$blower_det <= result$eff_matrix_det) , 
                          (result$bupper_sto >= result$eff_matrix_sto) * (result$blower_sto <= result$eff_matrix_sto) , 
                          (result$bupper_det_d >= result$eff_matrix_det) * (result$blower_det_d <= result$eff_matrix_det) , 
                          (result$bupper_sto_d >= result$eff_matrix_sto) * (result$blower_sto_d <= result$eff_matrix_sto))
    
    within_q_b <- cbind(within_q_b, (result$bupper_det >= result$eff_matrix_det) * (result$blower_det <= result$eff_matrix_det))
    within_es_b <- cbind(within_es_b, (result$bupper_sto >= result$eff_matrix_sto) * (result$blower_sto <= result$eff_matrix_sto))
    within_q_db <- cbind(within_q_db, (result$bupper_det_d >= result$eff_matrix_det) * (result$blower_det_d <= result$eff_matrix_det))
    within_es_db <- cbind(within_es_db, (result$bupper_sto_d >= result$eff_matrix_sto) * (result$blower_sto_d <= result$eff_matrix_sto))

    colnames(length_all) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap", "Q", "ESARSA")
    colnames(within_bootstrap) <- c("Q-bootstrap", "ESARSA-bootstrap", "Q-Dbootstrap", "ESARSA-bootstrap")
    length_list[[s]] <- length_all
    within_list[[s]] <- within_bootstrap 
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



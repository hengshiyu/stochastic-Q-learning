# library of functions


expit <- function(o_v, a_v, delta_1, delta_2){
    exp_v <- exp(-(delta_1 * o_v + delta_2 * a_v))
    p_v <- 1/(1+exp_v)
    return(p_v)
}

softmax <- function(action_vec, sec_part){
    action_matrix <- as.matrix(action_vec)
    sec_matrix <- as.matrix(sec_part)
    output <- sec_matrix %*% t(action_matrix)
    output <- exp(output)
    output0 <- output
    for(i in 1:dim(output)[1]){
        output0[i, ] <- output[i, ]/sum(output[i, ])
    }
    return(output0)
}
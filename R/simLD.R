#' Simulate spontaneous mutation model
#'
#' Simulate a Luria-Delbruck experiment assuming spontaneous mutation of resistance to phages
#' @param n_gens Number of generations in tube
#' @param mut_rate Mutation rate
#' @param ncells_init Initial number of cells in a tube transferred from a flask before growth
#' @param n_sample Number of cells to plate from a tube after growth
#' @param n_plates Number of replicate plates per experiment.
#' @param ncells_res_init Number of mutant cells in a tube before growth representing standing variation
#' @return A named vector of four: mean number of resistant colonies in experiment A (plating on n_plates plates from one tube), mean number of resistant colonies in experiment B (plating on n_plates plates from n_plates tubes), variance in experiment A, and variance in experiment B.
#' @export
simLD_spo <- function(n_gens, mut_rate , ncells_init , n_sample , n_plates,  ncells_res_init = 0){
        tube_count_a <- sim_tube_count(n_gens = n_gens, 
                                       mut_rate = mut_rate, 
                                       ncells_init = ncells_init)
        tubes_count_b <- sapply(1:n_plates, 
                                function(x){
                                        sim_tube_count(n_gens = n_gens, 
                                                       mut_rate = mut_rate, 
                                                       ncells_init = ncells_init, 
                                                       ncells_res_init = ncells_res_init)
                                }
        )
        plates_count_a <- sim_plate_count(n_re = tube_count_a, 
                                          n_wt = ncells_init * 2^n_gens - tube_count_a, 
                                          n_plates = n_plates,
                                          n_sample = n_sample
        )
        plates_count_b <- sapply(tubes_count_b, 
                                 function(x){
                                         sim_plate_count(n_re = x,
                                                         n_wt = ncells_init * 2^n_gens - x,
                                                         n_plates = 1,
                                                         n_sample = n_sample
                                         )
                                 }
        )
        result <- c(mean_a = mean(plates_count_a), 
                    mean_b = mean(plates_count_b), 
                    var_a = var(plates_count_a), 
                    var_b = var(plates_count_b))
        return(result)
}


#' Simulate induced mutation model
#'
#' Simulate a Luria-Delbruck experiment assuming induced mutation by exposure to phages
#' @param n_plates Number of replicate plates per experiment.
#' @param mut_rate Mutation rate
#' @param n_sample Number of cells to plate from a tube after growth
#' @return A named vector of four: mean number of resistant colonies in experiment A (plating on n_plates plates from one tube), mean number of resistant colonies in experiment B (plating on n_plates plates from n_plates tubes), variance in experiment A, and variance in experiment B.
#' @export
simLD_ind <- function(n_plates, n_sample, mut_rate){
        plates_a <- rbinom(n_plates, n_sample, mut_rate)
        plates_b <- rbinom(n_plates, n_sample, mut_rate)
        res <- c(mean_a = mean(plates_a),
                 mean_b = mean(plates_b),
                 var_a = var(plates_a),
                 var_b = var(plates_b)
        )
        return(res)
}





#' Simulate growth
#'
#' Simulate cell growth in a tube
#' @param n_gens Number of generations in tube
#' @param mut_rate Mutation rate
#' @param ncells_init Initial number of cells in a tube transferred from a flask before growth
#' @param ncells_res_init Number of mutant cells in a tube before growth representing standing variation
#' @return A vector of integers for numbers of colonies on plates with phages
#' @export
sim_tube_count <- function(n_gens, mut_rate, ncells_init, ncells_res_init = 0){
        n_wt <- ncells_init - ncells_res_init
        n_re <- ncells_res_init
        for(t in 1:n_gens){
                n_re <- 2 * n_re + rbinom(1, 2 * n_wt, mut_rate)
                n_wt <- ncells_init * 2^t - n_re
        }
        return(n_re)
}



#' Simulate plating
#'
#' Simulate plating cell culture from a tube containing phage-resistant and wild type (phage-sensitive) cells to plates with phages
#' @param n_re Number of resistant cells in the tube
#' @param n_wt Number of wild type cells in the tube
#' @param n_plates Number of replicate plates per experiment
#' @param n_sample Number of cells to plate from a tube 
#' @return A vector of integers for numbers of colonies on plates with phages
#' @export
sim_plate_count <- function(n_re, n_wt, n_plates, n_sample){
        # Number of resistant cells in plates
        plates_count <- c()
        if(length(n_sample) == 1){
                n_sample <- rep(n_sample, n_plates)
        }
        # Loop over n_plates plates
        for(i in 1:n_plates){
                plates_count <- c(plates_count, rhyper(1, n_re, n_wt, n_sample[i]))
                n_re <- n_re - plates_count[i]
                n_wt <- n_wt - (n_sample[i] - plates_count[i])
        }
        return(plates_count)
}


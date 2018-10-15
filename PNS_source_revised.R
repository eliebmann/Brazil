####Defining functions for analysis and presentation of data###########
Svy_fast <- function(x, obj, y) {
  require(dplyr)
  #calculating the proportion/mean by MH dx after account for design effects and weights
  res <- svyby(as.formula(paste0("~", x)), by= as.formula(paste0("~",y)), obj, svymean, ci = TRUE, vartype = "ci")
  
  #Now, trying to get the DF returned from svby into a nicer data frane
  #Forß there is a 1:3 factor level:column ratio for the SVBY, mean + CI output.
  # The first column is for MH_dx. The point estimates are displayed in order, followed by the CIs.  
  if(is.factor(dat_final[,x])) {
    #print(res)
    len <- length(levels(dat_final[,x]))
    colnames(res)[c(2:(len+1))] <- paste0("%", levels(dat_final[,x]))
    #colnames(res)[c(2:(len+1))] <- levels(dat_final[,x])
    colnames(res)[c(len+2):((len*3)+1)] <- paste(levels(dat_final[,x]), "CI", sep = " ")
    
    #Now, reordering columns in effort to brack CIs and reduce the returned DF.
    res <- res[, order(names(res))]
    #Making sure that the MH_dx comes at the end of the DF.
    #res<- res %>%select(-y, y)
    res <- res[,c(y, names(res)[-which(names(res)==y)])]
    #Re-naming columns for easier manipulation
   # print(res)
    CI_names <- names(res)[grep("CI", names(res))][order(names(res)[grep("CI", names(res))])]
  # print(CI_names)
    CI_names1 <- CI_names[!grepl("\\.", CI_names)]
    #print(CI_names1)
    index_df <- as.data.frame(sapply(CI_names1, function(x) which(grepl(x, names(res)))))
    
    ##Multiplying the point estimates by 100.
    Mult_names <- names(res)[!grepl("CI", names(res)) & !names(res) %in% y]
    res[, Mult_names] <- lapply(res[, Mult_names], function(x) round(x*100, 3))  
    #Function that re-scales and brackets the CI for each level of the variable.
    bracket_sub <- function(x) {
      brack <- data.frame(CI = rep(NA, nrow(res)))
      for(i in 1:nrow(res)) {
        brack[i, ] <- paste0(paste(paste0(paste0("[", 100*round(res[i,x[1]], 3)), ","), 100*round(res[i, x[2]],3),sep= " "), "]")
      }
      return(brack)
    }
    brack_res <- lapply(index_df, bracket_sub)
    #print(brack_res)
    #Takes the list output and collapses it into a single dataframe
    brack_res <- do.call(cbind, lapply(brack_res, function(x) x))
    colnames(brack_res) <- paste(CI_names1, "95%CI", sep = " ")
    #Eliminating the original CI columns from the svyby output to make way for the bracketed one. 
    res_df <- res[,-which(names(res) %in% CI_names)]
    res_df <- cbind(res_df, brack_res)
    res_df <- res_df[, order(names(res_df))]
    #Ensuring the MH_dx variable comes first. 
    #res_df <- res_df %>% select(y, everything())
    res_df <- res_df[,c(y, names(res_df)[-which(names(res_df)==y)])]
    
  } else {
    #This is code for the condition in which the var is numeric. 
    colnames(res)[2] <- x
    colnames(res)[3] <- paste(x, "CI LL", sep = " ")
    colnames(res)[4] <- paste(x, "CI UL", sep = " ")
    
    brack <- data.frame("CI" = rep(NA, nrow(res)))
    for(i in 1:nrow(res)) {
      brack[i, ] <- paste0(paste(paste0(paste0("[", round(res[i,3], 2)), ","), round(res[i, 4],2),sep= " "), "]")
    }
    res_df <- res[, -which(grepl("CI", names(res)))]
    res_df <- cbind(res_df, brack)
  #  res_df <- res_df %>% select(!! y_q, everything())
    #Now,
    
  }
  return(res_df)
  
}

#This function is help clean up the output from the first lapply function for table 1.
DF_reshuffle <- function(x) {
  require(dplyr)
  if(ncol(x) > 3) {
    require(reshape2)
    res1 <- as.data.frame(t(x))
    res1 <- res1[-1, ]
    res1$id <- gsub('^[%]*| .*\\S', '', rownames(res1))
    res_m <- melt(res1, id.vars = "id")
    res_m$log <- grepl("\\[", res_m$value)
    res_cast <- dcast(res_m, value + id ~ variable + log, drop = T)[,-1]
    res_cast <-res_cast %>%
      group_by(id) %>%
      summarise_each(funs(first(.[!is.na(.)]))) 
    colnames(res_cast) <- c("Var", "No DX %", "95% CI", "DX %", "95% CI")
  } else {
    res1 <- as.data.frame(t(x))
    res1 <- res1[-1, ]
    res1$id <-rownames(res1)
    res_m <- melt(res1, id.vars = "id")
    res_m$log <- grepl("\\[", res_m$value)
    res_cast <- dcast(res_m, value + id ~ variable + log, drop = T)[,-1]
    res_cast$ID <- rep(colnames(x)[2], nrow(res_cast))
    res_cast <-res_cast %>%
      group_by(ID) %>%
      summarise_each(funs(first(.[!is.na(.)]))) 
    res_cast <- res_cast[,-2]
    colnames(res_cast) <- c("Var", "No DX %", "95% CI", "DX %", "95% CI")
  }
  return(res_cast)
}

#Function for chi-sq and t tests#######
Svy_test <- function(x,y,obj) {
  if(!is.numeric(dat_final[,y])) {
    chi_res <- svychisq(as.formula(paste0(paste0(paste0("~", x),"+"), y)), obj)
    test_df <- data.frame("Pval" = chi_res$p.value)
    
    if(test_df[1, "Pval"] < .05 & test_df[1, "Pval"] >= .01) {
      test_df$Ast[1] <- "*"
    } else if (test_df[1, "Pval"] < .01 & test_df[1,"Pval"] >= .001) {
      test_df$Ast[1] <- "**"
    } else if(test_df[1, "Pval"] < .001) {
      test_df$Ast[1] <- "***"
    } else {
      test_df$Ast[1] <- ""
    }
    test_df$Pval <- round(test_df$Pval, 3)
    
  } else {
    t_res <- svyttest(as.formula(paste0(paste0(x,"~"),y)), obj) 
    test_df <- data.frame("Pval" = t_res$p.value)
    
    if(test_df[1, "Pval"] < .05 & test_df[1, "Pval"] >= .01) {
      test_df$Ast[1] <- "*"
    } else if (test_df[1, "Pval"] < .01 & test_df[1,"Pval"] >= .001) {
      test_df$Ast[1] <- "**"
    } else if(test_df[1, "Pval"] < .001) {
      test_df$Ast[1] <- "***"
    } else {
      test_df$Ast[1] <- ""
    }
    test_df$Pval <- round(test_df$Pval, 3)
  }
  return(test_df)
}


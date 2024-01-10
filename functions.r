# Fucntions to be used in all scripts

# paired Wilcoxon test
# k is a column in dat
pairedWilcox = function(k, dat)
{
  # create matrix patients vs time points  
  library(Matrix)
  mat = sparseMatrix(i = as.numeric(as.factor(dat$patientID)), 
                     j = as.numeric(as.factor(dat$timePoint)), x = 1, 
                     dimnames = list(levels(as.factor(dat$patientID)),
                                     levels(as.factor(dat$timePoint))))
  mat = as.data.frame(as.matrix(mat))
  
  # running paired test
  # Baseline vs C1D1
  p1 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,2]), paired = T)$p.value
  #  C1D1 vs  C2D1
  p2 = wilcox.test(as.numeric(mat[,2]), as.numeric(mat[,3]), paired = T)$p.value
  #  Baseline vs  C2D1
  p3 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,3]), paired = T)$p.value
  
  # mean of delta between time points
  d1 = mean(as.numeric(mat[,2]) -  as.numeric(mat[,1]), na.rm = T)
  d2 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,2]), na.rm = T)
  d3 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,1]), na.rm = T)
  
  
  return(c(pValue_BvsT1 = p1, 'mean(T1-B)' = d1, pValue_T1vsT21 = p2, 
           'mean(T2-T1)' = d2, pValue_BvsT2 = p3,'mean(T2-B)' = d3))
}

# run tests
#Wilcoxon test for each time point comparing responders vs non-responders. 
# And paired Wilcoxon p-values for baseline vs post run-in and 
# post run-in vs week 8.
#Then take differences between all time points and 
#compare responders vs non-responders to see 
#if there is significant difference in changes after treatment
runWilcox = function(k, dat, dta_pts)
  {
  mat <- data.frame(matrix(nrow = length(unique(dat$patientID)),
                           ncol = 3, dimnames = list(unique(dat$patientID),
                                                     levels(factor(dat$timePoint)))), check.names = F)
  for(i in levels(factor(dat$timePoint)))
  {
    d1 = dat %>% filter(timePoint == i)
    mat[d1$patientID, i] = d1[,k]
  }
  # add response
  mat$RECIST = dta_pts[rownames(mat),"RECIST"]
  # create treatment effect 
  # substract baseline from post run-in and week 8 
  mat$delta_T1_B =  mat[,2]-mat[,1]
  mat$delta_T2_B =  mat[,3]-mat[,1]
  # and week 8 from post run-in 
  mat$delta_T2_T1 =  mat[,3]-mat[,2]
  
  # running paired test
  # Baseline vs post run-in
  p1 = wilcox.test(mat[,1], mat[,2], paired = T)$p.value
  # post run-in vs week 8
  p2 = wilcox.test(mat[,2], mat[,3], paired = T)$p.value
  # baseline vs week 8
  p2_2 = wilcox.test(mat[,1], mat[,3], paired = T)$p.value
  
  # test by response
  # for baseline 
  p3 = wilcox.test(Baseline ~ RECIST, data = mat)$p.value
  # for post run-in 
  p4 = wilcox.test(mat[,2] ~ mat[,"RECIST"], data = mat)$p.value
  # for week 8
  p5 = wilcox.test(mat[,3] ~ mat[,"RECIST"], data = mat)$p.value
  # compare differences in treatment effect in responders vs non-responders
  p6 = wilcox.test(delta_T1_B ~ RECIST, data = mat)$p.value
  p7 = wilcox.test(delta_T2_B ~ RECIST, data = mat)$p.value
  p8 = wilcox.test(delta_T2_T1 ~ RECIST, data = mat)$p.value
  
  res = c(BvsT1 = p1,T1vsT2 = p2, BvsT2 = p2_2,
          B_RvsNR = p3, T1_RvsNR = p4, T2_RvsNR = p5, 
           d_T1_B_RvsNR = p6, d_T2_B_RvsNR = p7, d_T2_T1_RvsNR = p8,
          mean_d_T1_B_R = mean(mat %>% filter(RECIST == "Responder") %>% select(delta_T1_B) %>% unlist(), na.rm = T),
          mean_d_T1_B_NR = mean(mat %>% filter(RECIST == "Non-responder") %>% select(delta_T1_B) %>% unlist(), na.rm = T),
          mean_d_T2_B_R = mean(mat %>% filter(RECIST == "Responder") %>% select(delta_T2_B) %>% unlist(), na.rm = T),
          mean_d_T2_B_NR = mean(mat %>% filter(RECIST == "Non-responder") %>% select(delta_T2_B) %>% unlist(), na.rm = T),
          mean_d_T2_T1_R = mean(mat %>% filter(RECIST == "Responder") %>% select(delta_T2_T1) %>% unlist(), na.rm = T),
          mean_d_T2_T1_NR = mean(mat %>% filter(RECIST == "Non-responder") %>% select(delta_T2_T1) %>% unlist(), na.rm = T))
  return(round(res,3))
}

#===================================
##Create All Facet for plotting
CreateAllFacet <- function(df, col){
  df$facet <- df[[col]]
  temp <- df
  temp$facet <- "All"
  merged <-rbind(temp, df)
  
  # ensure the facet value is a factor
  merged[[col]] <- as.factor(merged[[col]])
  
  return(merged)
}

#==================================
# print the number of pairs between time points
# the input table should have patientID and timePoint columns
printPairs = function(dat)
{
  library(Matrix)
  times = levels(as.factor(dat$timePoint))
  mat = sparseMatrix(i = as.numeric(as.factor(dat$patientID)), 
                     j = as.numeric(as.factor(dat$timePoint )), x = 1, 
                     dimnames = list(levels(as.factor(dat$patientID)),
                                     times))
  mat = as.data.frame(as.matrix(mat))
  
  cat("The number of paired samples between", times[1],"and", times[2],":",
      mat[which(mat[,times[1]] > 0 & mat[,times[2]] > 0),] %>% nrow,"\n")
  cat("The number of paired samples between", times[1],"and", times[3],":",
      mat[which(mat[,times[1]] > 0 & mat[,times[3]] > 0),] %>% nrow,"\n")
  cat("The number of paired samples between", times[3],"and", times[2],":",
      mat[which(mat[,times[3]] > 0 & mat[,times[2]] > 0),] %>% nrow)
  
  return(mat)
  
}
# function from CG to calculate ORR
get_freq_table <- function(dta_pts, v_name) {
  f_tbl <- function(dta) {
    rst <- dta %>%
      group_by(Subtype, Resp) %>%
      summarize(N = n()) %>%
      left_join(dta %>%
                  group_by(Subtype) %>%
                  summarize(Total = n())) %>%
      mutate(Freq = 100 * N / Total) %>%
      select(Subtype, Resp, N, Total, Freq) %>%
      data.frame()
    
    ci <- NULL
    for (i in seq_len(nrow(rst))) {
      cur_ci <- binom.test(rst[i, "N"], rst[i, "Total"])$conf.int
      ci     <- round(rbind(ci, 100 * cur_ci),2)
    }
    
    colnames(ci) <- c("CI_Lb", "CI_Ub")
    cbind(rst, data.frame(ci))
  }
  
  dta_pts$Resp <- dta_pts[[v_name]]
  dta_pts_2    <- dta_pts %>%
    mutate(Subtype = "Overall")
  
  rst <- rbind(f_tbl(dta_pts), f_tbl(dta_pts_2))
}

#=======================
# it finds the minimum Y per patient
# what is Y I'm not sure
get_tumor_best <- function(dat_tumor) {
  dta_tumor %>%
    filter(0 != Days) %>%
    group_by(ID) %>%
    filter(Y == min(Y)) %>%
    mutate(Y = Y * 100) %>%
    arrange(desc(Y))
}

##========================
# modified functions from CG
##========================
get_median_txt <- function(rst_survival, unit = "Months ", digits = 1) {
  rst_sum <- rbind(summary(rst_survival)$table)
  print(rst_sum)
  grps    <- sapply(rownames(rst_sum),
                    function(x) gsub("Subtype=|Group=", "", x))
  # median survival by group
  med     <- round(rst_sum[, "median"], digits)
  # lower 95% CI
  ci_l <- round(rst_sum[, "0.95LCL"], digits)
  # upper 95% CI
  ci_u <- round(rst_sum[, "0.95UCL"], digits)
  
  gm<- apply(cbind(grps, med, ci_l, ci_u), 1,
             function(x) paste0(x[1], " = ", x[2], ' (95% CI = [', x[3],', ', x[4],'])'))
  
  median_txt <- paste("Median Survival in ", unit, " \n ",
                      paste(gm, collapse = "\n"), sep = "")
  
  median_txt
}

# returns survival rate at specified time
get_surv_txt <- function(rst_survival, time = 6, digits = 0) {
  rst_sum <- summary(rst_survival, times = time)
  grps    <- gsub("Subtype=|Group=", "",rst_sum$strata)
  # survival by group
  surv_rate     <- round(rst_sum$surv*100, digits)
  # lower 95% CI
  ci_l <- round(rst_sum$lower*100, digits)
  # upper 95% CI
  ci_u <- round(rst_sum$upper*100, digits)
  
  # gm<- apply(cbind(grps, surv_rate, ci_l, ci_u), 1,
  #           function(x) paste0(x[1], " = ", x[2], ' (95% CI = [', x[3],', ', x[4],'])'))
  gm<- apply(cbind(surv_rate, ci_l, ci_u), 1,
             function(x) paste0( x[1], '% [', x[2],'%, ', x[3],'%]'))
  gm
}

## plot survival curves
plot_surv <- function(dat, type_surv = "PFS", title = "",
                      add_overall = TRUE, xlim = 70, xtext = 50,
                      x_6m = 9, x_12 = 18, add6m = FALSE) {
  
  # p-value for HR+ vs TNBC
  dat$time  <- dat[[paste0(type_surv, "_Time")]]
  dat$event <- dat[[paste0(type_surv, "_Event")]]
  rst_surv_two_groups <- survfit(Surv(time, event) ~ Subtype, data = dat)
  #extract p-value
  p_val = round(surv_pvalue(rst_surv_two_groups, data = dat)$pval,3)
  #     print(p_val)
  
  # add All as a subtype
  dat       <- CreateAllFacet(dat, "Subtype")
  dat$time  <- dat[[paste(type_surv, "_Time",  sep = "")]]
  dat$event <- dat[[paste(type_surv, "_Event", sep = "")]]
  
  rst_survival <- survfit(Surv(time, event) ~ Subtype, data = dat)
  # get text for median survival
  median_txt   <- get_median_txt(rst_survival)
  # get text for 6 months
  text_6m = get_surv_txt(rst_survival)
  #get text for 12 months
  text_12m = get_surv_txt(rst_survival, 12)
  
  attr(rst_survival$strata, "names") <-
    sapply(attr(rst_survival$strata, "names"),
           function(x) gsub("Subtype=", "", x))
  
  rst_plot     <- ggsurvplot(rst_survival,
                             data = dat,
                             risk.table = T,
                             break.time.by = 6,
                             legend = "top",
                             xlab = "Time (in months)",
                             ylab = "Survival Probability", conf.int = FALSE,
                             title = title,
                             palette =c(subtypeCol, "Overall" = 'black'),
                             ylim = c(0, 1), xlim = c(0, xlim),
                             legend.title = "",
                             tables.y.text = FALSE)
  
  #    rst_plot$table <- rst_plot$table
  
  rst_plot$plot <- rst_plot$plot +
    # add text with median survival
    ggplot2::annotate("text", x = xtext, y = 0.65,
                      size = 3.5, label = median_txt) + 
    # add horizontal line at 50% survival rate
    geom_hline(yintercept = .5, linetype="dashed", color = 'grey25',)
  
  if(add6m) rst_plot$plot <- rst_plot$plot +
    # add vertical line at 6 months
    geom_vline(xintercept = 6, linetype="dashed", color = 'grey25')+
    # add survival values for groups at 6 months
    annotate(geom="text", x=x_6m, y=1, label= text_6m[1],color = subtypeCol["HR+"], fill= 'white', )+
    annotate(geom="text", x=x_6m, y=.96, label= text_6m[2],color = subtypeCol["TNBC"])+
    annotate(geom="text", x=x_6m, y=.92, label= text_6m[3],color = 'black')+
    # add vertical line at 12 months
    geom_vline(xintercept = 12, linetype="dashed", color = 'grey25')+
    # add survival values for groups at 6 months
    annotate(geom="text", x=x_12, y=.86, label= text_12m[1],color = subtypeCol["HR+"])+
    annotate(geom="text", x=x_12, y=.82, label= text_12m[2],color = subtypeCol["TNBC"])+
    annotate(geom="text", x=x_12, y=.78, label= text_12m[3],color = 'black')
  
  
  # print p-value
  print(paste("The p-value for HR+ vs TNBC is", p_val))
  
  print(rst_plot)
}

pval_surv <- function(dat, type_surv = "PFS") {
  dat$time  <- dat[[paste(type_surv, "_Time",  sep = "")]]
  dat$event <- dat[[paste(type_surv, "_Event", sep = "")]]
  
  surv_dat <- dat %>%
    dplyr::filter(Subtype != "Overall") %>%
    data.frame()
  
  rst <- survfit(Surv(time, event) ~ Subtype,
                 data = surv_dat)
  surv_pvalue(rst, data = surv_dat)
}


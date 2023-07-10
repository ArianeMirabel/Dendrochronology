#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("Rgraphviz"))

#urlRF <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz"
#install.packages(urlRF, repos=NULL, type="source") 


#####
invisible(lapply(c("data.table","randomForestSRC", "ggplot2", "gridExtra", "ggrepel"),
                 function(pk){if(!pk %in% row.names(installed.packages())){install.packages(pk)} 
                   library(pk,character.only=T)}))

load("Data_RandomForest")

RF_data <- RF_data_allT[Pr...t.. <= 0.1,][,species := as.factor(species)][,
                          .(species, Moy_ba, Moy_age, Elevation, Slope, Orientation, MAT, MAP, SummerSMI, t.value)]
RF_data <- RF_data[!is.na(t.value)]

mtry <- tuneRF(RF_data[,.(species, Moy_ba, Moy_age, Elevation, Slope, Orientation, MAT, MAP, SummerSMI)],
               RF_data$t.value, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best_mtry <- max(3, mtry[mtry[, 2] == min(mtry[, 2]), 1])

rf <- rfsrc(t.value ~ .,data = RF_data, importance = T, mtry = best_mtry, nodesize=1)
 
### partial plots
Vars <- c("Moy_ba", "Moy_age", "Elevation", "MAT", "MAP", "SummerSMI")
pPlots <- plot.variable(rf, Vars, partial = T, show.plots = F)

save(pPlots, file = "RF_PartialPlots")

Sp <- unique(RF_data$species)

Ylab <- data.table(Unit= c("cm²", "year", "m", "°C", "mm", "mm"), Vars = Vars, 
                   Varnames = c("Mean BA", "Mean age", "Elevation", "MAT", "MAP", "Summer SMI"))
pPlots_reg <- lapply(1:length(Vars), function(i){
  
    pP <- pPlots[["pData"]][[i]]
    pP <- data.table(Var = pP$xvar.names, X = pP$x.uniq, Y = pP$yhat, sdY = pP$yhat.se)[,Ymin := Y-sdY][
            ,Ymax := Y+sdY]
    
    reg <- lm(pP$Y~pP$X)
    
    p <- ggplot(pP, aes(x=X, y=Y, ymin = Ymin, ymax = Ymax)) +
      geom_line(colour = "firebrick", lwd = 1, lty = 2) + 
      geom_point( colour = "firebrick") + 
      geom_ribbon( alpha=0.1, colour = "firebrick", linetype = 2) +
      geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], col = "black", lwd = 0.3, lty = 1) +
      ylab("") + theme_classic() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0,0.3,0,0.6), 'cm'), axis.title.y = element_text(vjust = 6, size = 7.5),
            plot.title = element_text(size =9), axis.text=element_text(size=8), plot.subtitle=element_text(size=9)) +
      ggtitle(Ylab[Vars == unique(pP$Var), Varnames], 
              subtitle = paste("x=", format(reg$coefficients[1], digits = 2), "+",
                               format(reg$coefficients[2], digits= 2, scientific = T))) 
    
    p_dens <- ggplot(data = data.frame(samples = pPlots[["pData"]][[i]]$x)) +
      geom_density(aes(x=samples, y=..scaled..), colour = "darkgrey", lwd= 1) + theme_classic() +
      theme(plot.margin = unit(c(0,0.3,0,0.6), 'cm'), axis.title.y = element_text(vjust = 8, size = 7.5)) +
      scale_y_continuous(breaks = c(0,1)) + xlab(Ylab[Vars == unique(pP$Var), Unit]) 
    ifelse (i %in% c(1,4), p_dens <- p_dens + ylab("Samples\ndensity") + 
              theme(plot.subtitle=element_text(size=9)) , p_dens <- p_dens + ylab("")) 
    
    gb2 <- ggplot_build(p_dens)
    lims <- data.table(X = gb2$data[[1]]$x, dens = gb2$data[[1]]$scaled)[, sumdens := cumsum(dens)][
      sumdens <= max(sumdens)*0.95 & sumdens >= max(sumdens)*0.05]
    lims <- data.frame(Xmin = lims[sumdens == min(sumdens),X], Xmax = lims[sumdens == max(sumdens),X])
    
    p <- p + geom_rect(aes(xmin=lims[["Xmin"]], xmax=lims[["Xmax"]], ymin=-Inf, ymax=Inf), 
                       colour = "grey50", alpha = 0.005)
    p_dens <- p_dens + geom_rect(aes(xmin=lims[["Xmin"]], xmax=lims[["Xmax"]]), 
                       color = "grey50", ymin=-Inf, ymax=Inf, alpha = 0.005)
    
    gb2 <- ggplot_build(p_dens)
    gb1 <- ggplot_build(p)
    n1 <- length(gb1$panel$ranges[[1]]$y.labels)
    n2 <- length(gb2$panel$ranges[[1]]$y.labels)
    gA <- ggplot_gtable(gb1)
    gB <- ggplot_gtable(gb2)
    g <- gtable:::rbind_gtable(gA, gB, "last")
    panels <- g$layout$t[grep("panel", g$layout$name)]
    g$heights[panels[1]] <- unit(5, "null"); g$heights[panels[2]] <- unit(1,"null")
 
    return(g)
  })

Arr <- do.call("arrangeGrob", c(pPlots_reg, ncol = 3)); grid.arrange(Arr)

#Variables importance
min_depth_frame <- min_depth_distribution(rf)

importance_frame <- measure_importance(rf)












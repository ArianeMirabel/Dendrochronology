invisible(lapply(c("data.table","mgcv", "stringr","sf","sp", "lme4", "rgeos","nlme",
                   "ggplot2", "reshape2", "gridExtra"),function(pk){
                     if(!pk %in% row.names(installed.packages())){install.packages(pk)}
                     library(pk,character.only=T)}))

#Load dendrochronology series
tmp.data <- as.data.table(read.csv("SitesData_VPDresponse.csv", sep = ";"))
species.run <- unique(tmp.data$species)

# Run and save species-specific models, 2-level lists 'species level/ site level)
VPD_climate <- lapply(species.run, function(sp){
  print(sp)
  
  R2.all <- data.table()
  ptable.all <- data.table()
  apres.detrend <- data.table()
  
  tmp.data.sp <- tmp.data[species == sp,]
  sites.run.sp <- unique(tmp.data.sp$uid_site)
  
  # Site-level
  Ret <- lapply(sites.run.sp, function(si){
    tmp.data.si <- tmp.data.sp[uid_site == si,][,uid_tree := factor(uid_tree)]
    
    # Sites statistics
    tmp.site <- tmp.data.si[, list(Length.tree = .N, Year_from = min(Year), Year_to = max(Year), age_tree = max(Age)), 
                            by = .(species, uid_site, uid_tree)][, Ntrees := .N, by = .(species, uid_site)]
    
    # initialize checks for models run
    model.sel <- 0
    stopp.rc <- 0
    
    # First detrend tree-ring series, removing the age effect
    detrend <- tryCatch(
      {
        if (unique(tmp.site$Ntrees) > 1)
          model.detrend.tree_UID <- gamm(LBai.05 ~ log(BA.t_1) + s(Age, bs='cr'),
                                       random = list(uid_tree=~1),
                                       correlation = corCAR1(value = 0.5, form = ~ Year | uid_tree), 
                                       method = 'REML', data = tmp.data.si)
        
        else 
          model.detrend.tree_UID <- gamm(LBai.05 ~ log(BA.t_1) + s(Age, bs='cr') + 
                                           SMIyear, data = tmp.data.si)
      },
      
      error = function(e) e
    )
  
  tmp.data.si$Resid_Lbai_SMI <- residuals(detrend$lme)
  
  # Fit the GAMM growth-VPD models on the detrended data
  result <- tryCatch(
    {
      if (unique(tmp.site$Ntrees) > 1)
        model.vpd.tree_UID <- gamm(Resid_Lbai_SMI ~ VPDyear + VPDyear_1,
                                     random = list(uid_tree=~1),
                                     method = 'REML', data = tmp.data.si)
      
      else 
        model.vpd.tree_UID <- gamm(Resid_Lbai_SMI ~ VPDyear + VPDyear_1 + s(Age, bs='cr'), 
                                     control = lmeControl(opt = "optim"), data = tmp.data.si)
    },
    
    error = function(e) e
  )
  

    if (!is.null(result$message)) stopp.rc <- 1 else {
      if (summary(model.vpd.tree_UID$gam)$r.sq < 0) stopp.rc <- 1}
    
    if (stopp.rc == 0) {
      
      model.sel <- 1
      mi10 <- model.vpd.tree_UID
      rm(model.vpd.tree_UID)
      
    } else { model.sel <- 0 }
    
    rm(stopp.rc)
    
    # When convergent, store all models details: 
    #predicted values, normal, pearson and normalized residuals, MSE, baskerville estimates
    if (model.sel > 0){
      
      tmp.data.si$pred.Lbai.05.SMI = predict(mi10$lme)
      tmp.data.si$res.Lbai.05.SMI = residuals(mi10$lme)
      tmp.data.si$res.Lbai.pearson.SMI = residuals(mi10$lme, type = "pearson")
      tmp.data.si$res.Lbai.normalized.SMI = residuals(mi10$lme, type = "normalized")
      
      mse = mean(tmp.data.si$res.Lbai.05.SMI^2)
      tmp.data.si$CF = exp(mse/2)
      
      tmp.data.si[,pred.Bai.bv := CF*exp(pred.Lbai.05.SMI) - 0.05]
      tmp.data.si[,res.Bai.bv := Bai-pred.Bai.bv]
      
      tmp.data.si[, ratio := Bai / pred.Bai.bv]
      
      # Store models parameters, summerize in a common table
      ptable <- data.table(Parameter = row.names(summary(mi10$gam)$p.table), summary(mi10$gam)$p.table)
      ptable$uid_site <- si
      ptable.all <- rbind(ptable.all, ptable)
      
      R2 <- data.table(tmp.site, R2 = summary(mi10$gam)$r.sq, model.sel = model.sel)
      R2.all <- rbind(R2.all,R2)
      
      apres.detrend <- rbindlist(list(apres.detrend,tmp.data.si), fill=TRUE)
      
      rm(mi10)
      rm(tmp.data.si, tmp.site, ptable, R2)
      
    } else {
      
      R2 <- data.table(tmp.site, R2 = NA_real_, model.sel = NA_integer_)
      R2.all <- rbind(R2.all,R2)
      
      ptable <- data.table("Parameter" = NA_character_, "Estimate"= NA_real_, "Std. Error"= NA_real_, "t value"=NA_real_, "Pr(>|t|)"=NA_real_)
      ptable$uid_site <- si
      ptable.all <- rbindlist(list(ptable.all, ptable), fill = TRUE)
      
      apres.detrend <- rbindlist(list(apres.detrend,tmp.data.si), fill=TRUE)
      
      rm(tmp.data.si, tmp.site, ptable, R2)
      
    }
    
    return(list(R2.all,ptable.all,apres.detrend))
    
  })
  write.csv(do.call(rbind,lapply(Ret, function(si) return(si[[1]]))), file = paste0("R2_VPDresponse_Normal",sp,".csv"))
  write.csv(do.call(rbind,lapply(Ret, function(si) return(si[[2]]))), file = paste0("ptable_VPDresponse",sp,".csv"))
  return(Ret)
})


# Get fitted models statistics: r², models prameters significance and frequency
#####
R2 <- do.call(rbind,lapply(species.run, function(sp){
  
  ptable <- fread(paste0("VPD_response_Species/ptable_VPDresponse",sp,".csv"))
  
  return(fread(paste0("VPD_response_Species/R2_VPDresponse",sp,".csv"))[,V1:=NULL][
    uid_site %in% ptable[Parameter %in% grep("VPDyear", ptable$Parameter, value = T) & Pr...t.. <= 0.1, uid_site]])
}))
R2.stat <- unique(R2[, meanR2 := round(mean(R2),2)][, meanR2.sp := round(mean(R2),2), by = "species"][
  , sdR2 := round(sd(R2), 2)][,.(species, meanR2.sp, meanR2, sdR2)])

# Parameters occurence
Sp_occurence <- do.call(rbind,lapply(species.run, function(sp){
  
  ptable <- fread(paste0("VPD_response_Species/ptable_VPDresponse",sp,".csv"))
  ptable <- ptable[grep("VPDyear", ptable$Parameter),]
  ptable[, Yr := gsub("year", "", Parameter)]
  
  Prop <- data.table(
    Tot = ptable[Pr...t.. <= 0.05, uniqueN(uid_site)],
    VPD = ptable[Yr == "VPD" & Pr...t.. <= 0.05, uniqueN(uid_site)],
    VPD_Pos = ptable[Pr...t.. <= 0.05 & Yr == "VPD" & Estimate >= 0, uniqueN(uid_site)],
    VPD_Neg = ptable[Pr...t.. <= 0.05 & Yr == "VPD" & Estimate <= 0, uniqueN(uid_site)],
    VPD_1 = ptable[Pr...t.. <= 0.05 & Yr == "VPD_1", uniqueN(uid_site)],
    VPD_1_Pos = ptable[Pr...t.. <= 0.05 & Yr == "VPD_1" & Estimate >= 0, uniqueN(uid_site)],
    VPD_1_Neg = ptable[Pr...t.. <= 0.05 & Yr == "VPD_1" & Estimate <= 0, uniqueN(uid_site)])[
      , species := sp]
  return(Prop)
}))[,.(species, Tot, VPD_1, VPD_1_Pos, VPD_1_Neg, VPD, VPD_Pos, VPD_Neg)]
#####

#Figures: Pie chart, maps and statistics from ArcGis output
#####
tmp.data <- as.data.table(read.csv("SitesData_VPDresponse.csv", sep = ";"))
species.run <- unique(tmp.data$species)

# Pie chart
Nsites.sp <- unique(unique(tmp.data[,.(uid_site, species)])[
  , Nsite.sp := .N, by = species][,.(species,Nsite.sp)])

VPDresp_data <- as.data.table(do.call(rbind,lapply(species.run, function(sp){
  ptable <- as.data.table(read.csv(paste0("VPD_response_Species/ptable_VPDresponse",sp,".csv")))
  ptable <- ptable[grep("VPDyear", ptable$Parameter),]
  ptable[, Yr := gsub("year", "", Parameter)]
  
  ptable <-ptable[,.(uid_site, Yr, t.value, Pr...t..)]
  setnames(ptable, "Yr", "Parameter")
  ptable$species <- sp
  return(ptable)
})))

Pie_data <- VPDresp_data[, value := paste0(sign(t.value),Parameter)]
Pie_data[Pr...t.. >= 0.05, value := "Non Significant"]
Pie_data[value == "-1VPD_1", value := "Previous VPD -"][
  value == "1VPD_1", value := "Previous VPD +"][value =="-1VPD", value := "VPD -"][
  value == "1VPD", value := "VPD +"][value == "NAVPD", value := "Non Significant"]

Pie_data <- unique(Pie_data[, Eff := .N, by = c("species", "value")][,.(species, value, Eff)])[,
                  EffTot := sum(Eff), by = species][, EffProp := Eff/EffTot, by = species]

Pie_data$value <- factor(Pie_data$value, levels= c("VPD +","Previous VPD +", "VPD -","Previous VPD -", 
                                                   "SMI +","SMI -", "T +", "T _","Non Significant"))
Pie_data <- Pie_data[order(species, value),]


sp.labs <- merge(data.frame(species = c("ABIELAS", "PICEENG" , "PICEGLA" , "PICEMAR", "PINUBAN",
                                        "PINUCON", "PINURES", "PSEUMEN","POPUTRE"),
         Name = c("Abies\nlasiocarpa\nN = ", "Picea\nengelmanii\nN = ", "Picea\nglauca\nN = ", 
             "Picea\nmariana\nN = ", "Pinus\nbanksiana\nN = ", "Pinus\ncontorta\nN = ", 
             "Pinus\nresinosa\nN =", "Pseudotsuga\nMenziesii\nN =",  "Populus\ntremuloides\nN = ")), Nsites.sp)
sp.labs <- setNames(paste0(sp.labs$Name, sp.labs$Nsite.sp, sp.labs$Sign), sp.labs$species)


ggplot(Pie_data, aes(x = "", y=EffProp, fill = value, alpha = value)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~species, labeller = labeller(species = sp.labs), nrow = 1) +
  scale_fill_manual(values = c("Previous VPD +" = "royalblue3",  "VPD +" = "royalblue3", 
                               "Previous VPD -" = "darkred", "VPD -" = "darkred",
                               "Non Significant" = "white")) +
  scale_alpha_manual(values = c("Previous VPD +" = 0.5,  "VPD +" = 1, 
                                "Previous VPD -" = 0.5, "VPD -" = 1, "Non Significant" = 1)) + 
  theme_void()+ theme(legend.position = "left") + guides(fill = guide_legend(title=""), alpha = guide_legend(title=""))


#Maps for ArcGis and figures from ArcGIs output
Sites <- fread("DB_VPD/tbl_Sites.csv")[,.(uid_site, latitude, longitude)]          
Sites <- Sites[which(complete.cases(Sites)),]

Sp_map_signif <- lapply(species.run, function(sp){
  ptable <- as.data.table(read.csv(paste0("VPD_response_Species/ptable_VPDresponse",sp,".csv")))
  
  ptable <- ptable[grep("VPDyear", ptable$Parameter),]
  ptable$Yr <- NA
  ptable$Yr[grep("VPDyear_1", ptable$Parameter)] <- "VPDyear_1"; ptable$Yr[grep("VPDyear:", ptable$Parameter)] <- "VPDyear"
  ptable[which(is.na(ptable$Yr))][Parameter == "VPDyear","Yr"] <- "VPDyear"
  
  ptable <- merge(ptable, Sites)[,.(uid_site, Yr, Estimate, Std..Error, t.value, Pr...t.., latitude, longitude)]
  setnames(ptable, "Yr", "Parameter")
  
  assign(paste0("VPDresponse_",sp), st_as_sf(ptable[Parameter == "VPDyear"], coords = c("longitude", "latitude"), 
                                             crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  assign(paste0("VPD_1response_",sp), st_as_sf(ptable[Parameter == "VPDyear_1"], coords = c("longitude", "latitude"), 
                                               crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  
  st_write(get(paste0("VPDresponse_",sp)), paste0("Maps/VPDresponse_Species/",sp,".shp"), 
           driver = "ESRI Shapefile")
  st_write(get(paste0("VPD_1response_",sp)), paste0("Maps/VPDresponse_Species/",sp,"_1.shp"), 
           driver = "ESRI Shapefile")
})

Whole_map_signif <- do.call(rbind,lapply(species.run, function(sp){
  ptable <- as.data.table(read.csv(paste0("VPD_response_Species/ptable_VPDresponse",sp,".csv")))
  ptable <- ptable[grep("VPDyear", ptable$Parameter),]
  ptable$Yr <- NA
  ptable$Yr[grep("VPDyear_1", ptable$Parameter)] <- "VPDyear_1"; ptable$Yr[grep("VPDyear:", ptable$Parameter)] <- "VPDyear"
  ptable[which(is.na(ptable$Yr))][Parameter == "VPDyear","Yr"] <- "VPDyear"
  
  ptable <- merge(ptable, Sites)[,.(uid_site, Yr, Estimate, Std..Error, t.value, Pr...t.., latitude, longitude)]
  setnames(ptable, "Yr", "Parameter")
  return(ptable)
}))

Whole_VPDresponse_signif <- st_as_sf(Whole_map_signif[Parameter == "VPDyear"], coords = c("longitude", "latitude"), 
                              crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
Whole_VPDresponse_1_signif <- st_as_sf(Whole_map_signif[Parameter == "VPDyear_1"], coords = c("longitude", "latitude"), 
                                crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
st_write(Whole_VPDresponse_signif, "Maps/Whole_Allt.shp", driver = "ESRI Shapefile")
st_write(Whole_VPDresponse_1_signif, "Maps/Whole_1_Allt.shp", driver = "ESRI Shapefile")


## Ecozone means barplots
tval_em <- as.data.table(read.csv("Maps/Tvalue_VPD.csv", sep = ";"))[,.(ZONE_NAME, MEAN)]
tval_em[ZONE_NAME == "Boreal PLain", ZONE_NAME := "Boreal Plain"]
setnames(tval_em, c("ZONE_NAME","MEAN"), c("Ecozone","Tval_mean"))
tval_em$Tval_mean <- round(as.numeric(gsub(",",".", tval_em$Tval_mean)),2)
tval_em$Sign <- as.factor(sign(tval_em$Tval_mean))
tval_em$Sign <- gsub("-1","Negative t-value", tval_em$Sign)
tval_em$Sign <- gsub("1","Positive t-value", tval_em$Sign)

tval_1_em <- as.data.table(read.csv("Maps/Tvalue_1_VPD.csv", sep = ";"))[,.(ZONE_NAME, MEAN)]
tval_1_em[ZONE_NAME == "Boreal PLain", ZONE_NAME := "Boreal Plain"]
setnames(tval_1_em, c("ZONE_NAME","MEAN"), c("Ecozone","Tval_mean"))
tval_1_em$Tval_mean <- round(as.numeric(gsub(",",".", tval_1_em$Tval_mean)),2)
tval_1_em$Sign <- as.factor(sign(tval_1_em$Tval_mean))
tval_1_em$Sign <- gsub("-1","Negative t-value", tval_1_em$Sign)
tval_1_em$Sign <- gsub("1","Positive t-value", tval_1_em$Sign)

pt <- ggplot(data=tval_em, aes(x=Ecozone, y=Tval_mean, fill = Sign)) +
  geom_bar(stat="identity") + ylim(-7,7) +
  theme_classic() + theme(axis.text.x=element_text(angle = 60, hjust = 1), legend.title = element_blank(), 
                          legend.position = "left", plot.margin=unit(c(0,0,0,0), 'line')) +
  labs(title = "", x = "", y = "")  +
  scale_fill_manual(values = c("Positive t-value"="royalblue3", "Negative t-value"="darkred")) 

pt_1 <- ggplot(data=tval_1_em, aes(x=Ecozone, y=Tval_mean, fill = Sign)) +
  geom_bar(stat="identity") + ylim(-7,7) +
  theme_classic() + theme(axis.text.x=element_text(angle = 60, hjust = 1), legend.position = "none",
                          plot.margin=unit(c(0,1,0,6), 'line')) +
  labs(title = "", x = "", y = "Mean t-value")  +
  scale_fill_manual(values = c("Positive t-value"="royalblue3", "Negative t-value"="darkred")) 

grid.arrange(pt_1,pt, ncol = 2)


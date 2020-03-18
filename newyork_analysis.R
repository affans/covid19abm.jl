rm(list = ls())

library(data.table)
library(ggplot2)
library(reshape)
library(dplyr)


## plot all the incidence prevalence
df_lat_prev = data.table()
for (isos in c("f00", "f01", "f02", "f03", "f04", "f05", "f06", "f07", "f08", "f09")){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_all.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_all.dat")    
  }
  
  df = fread(fp)
  df_lat_prev[, paste0("iso", isos) := df$lat_prev]
}
cnames = colnames(df_lat_prev)
df_lat_prev$time = 1:500
dfm = melt(df_lat_prev, id.vars = "time", measure.vars = cnames)
colors = c("#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#ffffff", "#e0e0e0", "#bababa", "#878787", "#4d4d4d", "#1a1a1a")
bvals = c("isof00", "isof01", "isof02", "isof03", "isof04", "isof05", "isof06", "isof07", "isof08", "isof09")
lvals = c("no isolation", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")
gg = ggplot(data=dfm)
gg = gg + geom_line(aes(x = time, y = value, color=variable), group=("variable"), size=1.2)
gg = gg + theme_bw()
gg = gg + ggtitle("Infection Prevalence (New York City)")
gg = gg + theme(text = element_text(size=15)) 
gg = gg + scale_color_manual(values=colors, name="Isolation Level", breaks=bvals, labels=lvals)
gg 

## plot hospitalization
df_hos_prev_all = data.table()
df_hos_prev_ag1 = data.table()
df_hos_prev_ag2 = data.table()
df_hos_prev_ag3 = data.table()
df_hos_prev_ag4 = data.table()
df_hos_prev_ag5 = data.table()
isolevels = c("f00", "f01", "f02", "f03", "f04", "f05", "f06", "f07", "f08", "f09")
## load up all, ag1, ag2, ag3, ... separately
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_all.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_all.dat")    
  }
  df = fread(fp)
  df_hos_prev_all[, paste0("iso", isos) := df$icu_prev]
}
## AG 1
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_ag1.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_ag1.dat")    
  }
  df = fread(fp)
  df_hos_prev_ag1[, paste0("iso", isos) := df$icu_prev]
}
## AG 2
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_ag2.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_ag2.dat")    
  }
  df = fread(fp)
  df_hos_prev_ag2[, paste0("iso", isos) := df$icu_prev]
}
## AG 3
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_ag3.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_ag3.dat")    
  }
  df = fread(fp)
  df_hos_prev_ag3[, paste0("iso", isos) := df$icu_prev]
}
## AG 4
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_ag4.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_ag4.dat")    
  }
  df = fread(fp)
  df_hos_prev_ag4[, paste0("iso", isos) := df$icu_prev]
}
## AG 5
for (isos in isolevels){
  if (isos == "f00"){
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau0_", isos, "/timelevel_ag5.dat")
  } else {
    fp = paste0("/data/covid19abm/simresults/newyork/b00485/tau1_", isos, "/timelevel_ag5.dat")    
  }
  df = fread(fp)
  df_hos_prev_ag5[, paste0("iso", isos) := df$icu_prev]
}
cnames = colnames(df_hos_prev_all)  ## column names are the same
df_hos_prev_all$time = 1:500
df_hos_prev_ag1$time = 1:500
df_hos_prev_ag2$time = 1:500
df_hos_prev_ag3$time = 1:500
df_hos_prev_ag4$time = 1:500
df_hos_prev_ag5$time = 1:500

df_hos_prev_all$styp = "all"
df_hos_prev_ag1$styp = "ag1"
df_hos_prev_ag2$styp = "ag2"
df_hos_prev_ag3$styp = "ag3"
df_hos_prev_ag4$styp = "ag4"
df_hos_prev_ag5$styp = "ag5"

df_hos_prev = rbind(df_hos_prev_all, df_hos_prev_ag1, df_hos_prev_ag2, 
                    df_hos_prev_ag3, df_hos_prev_ag4, df_hos_prev_ag5)

dfm = melt(df_hos_prev, id.vars = c("time", "styp"), measure.vars = cnames)
colors = c("#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#ffffff", "#e0e0e0", "#bababa", "#878787", "#4d4d4d", "#1a1a1a")
bvals = c("isof00", "isof01", "isof02", "isof03", "isof04", "isof05", "isof06", "isof07", "isof08", "isof09")
lvals = c("no isolation", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")
gg = ggplot(data=dfm)
gg = gg + geom_line(aes(x = time, y = value, color=variable), group=("variable"), size=1.2)
gg = gg + theme_bw()
gg = gg + ggtitle("ICU Prevalence (New York City)")
gg = gg + theme(text = element_text(size=15)) 
gg = gg + scale_color_manual(values=colors, name="Isolation Level", breaks=bvals, labels=lvals)
gg = gg + facet_wrap( ~ styp, ncol=3, scales="free")
gg 


### for scenarios only
s1 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_all.dat")
s2 = fread("/data/covid19abm/simresults/newyork/b00485/tau1_f05_q00/timelevel_all.dat")
s3 = fread("/data/covid19abm/simresults/newyork/b00485/tau1_f05_q09/timelevel_all.dat")
s4 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q09/timelevel_all.dat")

colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
lat_prev = data.table(t=1:500, s1=s1$lat_prev, s2=s2$lat_prev, s3=s3$lat_prev, s4=s4$lat_prev)
lpm = melt(lat_prev, id.vars = "t", measure.vars = c("s1", "s2", "s3", "s4"))
gg1 = ggplot(lpm)
gg1 = gg1 + geom_line(aes(x = t, y = value, color=variable), size=1.2)
gg1 = gg1 + theme_bw()
gg1 = gg1 + ggtitle("Total infections")
gg1 = gg1 + theme(text = element_text(size=15)) 
gg1 = gg1 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg1 = gg1 + ylab("Prevalence") + xlab("Time (in days)")
gg1 

hos_prev = data.table(t=1:500, s1=s1$hos_prev, s2=s2$hos_prev, s3=s3$hos_prev, s4=s4$hos_prev)
lpm = melt(hos_prev, id.vars = "t", measure.vars = c("s1", "s2", "s3", "s4"))
gg2 = ggplot(lpm)
gg2 = gg2 + geom_line(aes(x = t, y = value, color=variable), size=1.2)
gg2 = gg2 + theme_bw()
gg2 = gg2 + ggtitle("Hospitalization")
gg2 = gg2 + theme(text = element_text(size=15)) 
gg2 = gg2 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg2 = gg2 + ylab("Prevalence") + xlab("Time (in days)")
gg2 

icu_prev = data.table(t=1:500, s1=s1$icu_prev, s2=s2$icu_prev, s3=s3$icu_prev, s4=s4$icu_prev)
lpm = melt(icu_prev, id.vars = "t", measure.vars = c("s1", "s2", "s3", "s4"))
gg3 = ggplot(lpm)
gg3 = gg3 + geom_line(aes(x = t, y = value, color=variable), size=1.2)
gg3 = gg3 + theme_bw()
gg3 = gg3 + ggtitle("ICU")
gg3 = gg3 + theme(text = element_text(size=15)) 
gg3 = gg3 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg3 = gg3 + ylab("Prevalence") + xlab("Time (in days)")
gg3 

ded_prev = data.table(t=1:500, s1=s1$ded_prev, s2=s2$ded_prev, s3=s3$ded_prev, s4=s4$ded_prev)
lpm = melt(ded_prev, id.vars = "t", measure.vars = c("s1", "s2", "s3", "s4"))
gg4 = ggplot(lpm)
gg4 = gg4 + geom_line(aes(x = t, y = value, color=variable), size=1.2)
gg4 = gg4 + theme_bw()
gg4 = gg4 + ggtitle("Total Deaths")
gg4 = gg4 + theme(text = element_text(size=15)) 
gg4 = gg4 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg4 = gg4 + ylab("Total Deaths") + xlab("Time (in days)")
gg4 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(gg1)

gg1 <- gg1 + theme(legend.position="none")
gg2 <- gg2 + theme(legend.position="none")
gg3 <- gg3 + theme(legend.position="none")
gg4 <- gg4 + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(gg1, gg2, legend, gg3, gg4, nrow = 2, ncol=3)

max(s2$hos_prev)
sum(s1$ded_inc)/sum(s1$lat_inc)

sum(s1$ded_inc)/(sum(s1$inf_inc) + sum(s1$iiso_inc))


sum(s1$ded_inc)/(sum(s1$hos_inc)



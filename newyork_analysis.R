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

scenario_data <- function(x){
  s1 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_all.dat")
  s2 = fread("/data/covid19abm/simresults/newyork/b00485/tau1_f05_q00/timelevel_all.dat")
  s3 = fread("/data/covid19abm/simresults/newyork/b00485/tau1_f05_q09/timelevel_all.dat")
  s4 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q09/timelevel_all.dat")
  dt = data.table(s1=s1[, get(x)], s2=s2[, get(x)], s3=s3[, get(x)], s4=s4[, get(x)])
  dt = dt*800 ## 8 mil NYC pop
  dt$time = 1:500
  lpm = melt(dt, id.vars = "time", measure.vars = c("s1", "s2", "s3", "s4"))
  return(lpm)
}

colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

ti = scenario_data("lat_prev")
gg1 = ggplot(ti)
gg1 = gg1 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
gg1 = gg1 + theme_bw()
gg1 = gg1 + ggtitle("Total infections")
gg1 = gg1 + theme(text = element_text(size=15)) 
gg1 = gg1 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg1 = gg1 + ylab("Prevalence") + xlab("Time (in days)")
gg1 = gg1 + scale_y_continuous(labels = comma)
gg1 

hos = scenario_data("hos_prev")
gg2 = ggplot(hos)
gg2 = gg2 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
gg2 = gg2 + theme_bw()
gg2 = gg2 + ggtitle("Hospitalization")
gg2 = gg2 + theme(text = element_text(size=15)) 
gg2 = gg2 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg2 = gg2 + ylab("Prevalence") + xlab("Time (in days)")
gg2 = gg2 + scale_y_continuous(labels = comma)
gg2 


icu = scenario_data("icu_prev")
gg3 = ggplot(icu)
gg3 = gg3 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
gg3 = gg3 + theme_bw()
gg3 = gg3 + ggtitle("ICU")
gg3 = gg3 + theme(text = element_text(size=15)) 
gg3 = gg3 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg3 = gg3 + ylab("Prevalence") + xlab("Time (in days)")
gg3 = gg3 + scale_y_continuous(labels = comma)
gg3 

ded = scenario_data("ded_prev")
gg4 = ggplot(ded)
gg4 = gg4 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
gg4 = gg4 + theme_bw()
gg4 = gg4 + ggtitle("Total Deaths")
gg4 = gg4 + theme(text = element_text(size=15)) 
gg4 = gg4 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s4", "s2", "s3"), labels=c("No control", "Shelter in Place 60+ (SiP60+)", "Self-Isolation Symptomatic (SiS)", "SIP60+ (SiP60+ & SiS)"))
gg4 = gg4 + ylab("Total Deaths") + xlab("Time (in days)")
gg4 = gg4 + scale_y_continuous(labels = comma)
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

## Age Specific Case Fatality Rate
## Look at no isolation case
ag1 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_ag1.dat")
ag2 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_ag2.dat")
ag3 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_ag3.dat")
ag4 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_ag4.dat")
ag5 = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_ag5.dat")
sum(ag1$ded_inc)/sum(ag1$inf_inc)*100
sum(ag2$ded_inc)/sum(ag2$inf_inc)*100
sum(ag3$ded_inc)/sum(ag3$inf_inc)*100
sum(ag4$ded_inc)/sum(ag4$inf_inc)*100
sum(ag5$ded_inc)/sum(ag5$inf_inc)*100

## Overall CFR (base case no isolation)
all = fread("/data/covid19abm/simresults/newyork/b00485/tau0_f00_q00/timelevel_all.dat")
sum(all$ded_inc)/sum(all$lat_in)

## Plots for splitting age specific deaths
getagdata <- function(x) {
  ag1 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag1.dat"))
  ag2 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag2.dat"))
  ag3 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag3.dat"))
  ag4 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag4.dat"))
  ag5 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag5.dat"))
  dt = data.table(ag1=ag1$ded_prev, ag2=ag2$ded_prev, ag3=ag3$ded_prev, ag4=ag4$ded_prev, ag5=ag5$ded_prev)
  dt = dt * 800 ## 8 mil NYC pop
  dt$time = 1:500
  dt = dt[1:350, ]
  dtm = melt(dt, id.vars = "time", measure.vars = c("ag1", "ag2", "ag3", "ag4", "ag5"))
  return(dtm)
}

getplot <- function(x, scen){
  gg4 = ggplot(data=x)
  gg4 = gg4 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
  gg4 = gg4 + theme_bw()
  gg4 = gg4 + ggtitle(scen)
  gg4 = gg4 + theme(text = element_text(size=15)) 
  gg4 = gg4 + scale_color_manual(values=colors, name="Age Groups", breaks=c("ag1", "ag2", "ag3", "ag4", "ag5"), labels=c("0-4", "5-19", "19-49", "50-64", "65+"))
  gg4 = gg4 + ylab("Total Deaths") + xlab("Time (in days)")
  gg4 
  return(gg4)
}

s1ag = getagdata("tau0_f00_q00")
sg1 = getplot(s1ag, "No control")

s2ag = getagdata("tau1_f05_q00")
sg2 = getplot(s2ag, "Self-Isolation Symptomatic (SiS)")

s3ag = getagdata("tau1_f05_q09")
sg3 = getplot(s3ag, "SIP60+ (SiP60+ & SiS)")

s4ag = getagdata("tau0_f00_q09")
sg4 = getplot(s4ag, "Shelter in Place 60+ (SiP60+)")

# extract the legend and remove from the other plots
sglegend = get_legend(sg1)
sg1 <- sg1 + theme(legend.position="none")
sg2 <- sg2 + theme(legend.position="none")
sg3 <- sg3 + theme(legend.position="none")
sg4 <- sg4 + theme(legend.position="none")
# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(sg1, sg4, sglegend, sg2, sg3,  nrow = 2, ncol=3)


x = "tau0_f00_q00"
all = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_all.dat"))
ag1 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag1.dat"))
ag2 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag2.dat"))
ag3 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag3.dat"))
ag4 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag4.dat"))
ag5 = fread(paste0("/data/covid19abm/simresults/newyork/b00485/", x, "/timelevel_ag5.dat"))

sum(ag1$ded_inc)/sum(ag1$inf_inc)*100
sum(ag2$ded_inc)/sum(ag2$inf_inc)*100
sum(ag3$ded_inc)/sum(ag3$inf_inc)*100
sum(ag4$ded_inc)/sum(ag4$inf_inc)*100
sum(ag5$ded_inc)/sum(ag5$inf_inc)*100


(sum(ag3$ded_inc) + sum(ag4$ded_inc))/(sum(ag3$inf_inc) + sum(ag4$inf_inc))
(sum(ag3$ded_inc) + sum(ag4$ded_inc) + sum(ag5$ded_inc))/(sum(ag3$inf_inc) + sum(ag4$inf_inc) + sum(ag5$inf_inc))

sum(ag5$ded_inc)/(sum(ag5$ded_inc) + sum(ag3$ded_inc) + sum(ag4$ded_inc) + sum(ag1$ded_inc) + sum(ag2$ded_inc))


sum(all$ded_inc)/sum(all$inf_inc)






### school closure 
scenario_data <- function(x){
  s1 = fread("/data/covid19abm/schoolclosure/reg/timelevel_ag5.dat")
  s2 = fread("/data/covid19abm/schoolclosure/clo/timelevel_ag5.dat")
  
  dt = data.table(s1=s1[, get(x)], s2=s2[, get(x)])
  dt$time = 1:500
  lpm = melt(dt[1:200], id.vars = "time", measure.vars = c("s1", "s2"))
  return(lpm)
}

ti = scenario_data("lat_inc")
gg1 = ggplot(ti)
gg1 = gg1 + geom_line(aes(x = time, y = value, color=variable), size=1.2)
gg1 = gg1 + theme_bw()
gg1 = gg1 + ggtitle("Total infections")
gg1 = gg1 + theme(text = element_text(size=15)) 
gg1 = gg1 + scale_color_manual(values=colors, name="Scenario", breaks=c("s1", "s2"), labels=c("No control", "School Closure"))
gg1 = gg1 + ylab("Prevalence") + xlab("Time (in days)")
gg1 = gg1 + scale_y_continuous(labels = comma)
gg1 
sum(ti[variable == "s1"]$value)
sum(ti[variable == "s2"]$value)

colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")


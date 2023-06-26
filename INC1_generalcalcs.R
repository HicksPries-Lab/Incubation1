## SOIL INCUBATION CALCULATIONS, SEPT 2022
## Author: Michelle S. Wang, michelle.s.wang.th@dartmouth.edu

# Load packages + functions
library(tidyverse)
#library(SoilR)
library(FME)
library(ggrepel)
library(scales)
library(ggpubr)
library(car)

# Read in data and initialize to how much you want to see
data <- read.csv("IRGA_Measurements.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example


############################################################################
# GENERAL CALCULATIONS #######################################################
# Constants
R <- 82.05746   # [mL*atm/(K*mol)]

# Room Parameters
Pr <- .98 # [atm]
Tr <- 22 + 273  # [K]

# Jar/Soil Parameters
Vjar_P <- 473.176 - 46  # [mL] pint jar - filled sample cup, from Google Sheet 'Incubation Initializations <- Bulk Density' 
Vjar_V <- 473.176 - 57  # [mL] pint jar - filled sample cup, from Google Sheet 'Incubation Initializations <- Bulk Density' 
# Bdensity_P <- .822  # [g/cm^3] bulk density of Palouse soil
# Bdensity_V <- .651 # [g/cm^3] bulk density of Vershire soil 
# Wsoil <- 39       # [g] weight of soil + residue
# Vsoil_P <- Wsoil/Bdensity_P   # [mL = cm^3] Palouse
# Vsoil_V <- Wsoil/Bdensity_P   # [mL = cm^3] Vershire
# Vair_P <- Vjar_P-Vsoil_P        # [mL] Palouse
# Vair_V <- Vjar_V-Vsoil_V        # [mL] Vershire

# n [mol] air inside jar
n_P <- (Pr*Vjar_P) / (R*Tr)     # [mol] Palouse
n_V <- (Pr*Vjar_V) / (R*Tr)     # [mol] Vershire

# Moles/Mass of C inside jar
molmass_C <- 12.011*10^3  # [mg/mol] molar mass of C
data_C <- data %>%
  mutate(moles_C_P = C_ppm*n_P/(10^6)) %>%  # [mol] moles of C in air in jar
  mutate(moles_C_V = C_ppm*n_V/(10^6)) %>%  # [mol] moles of C in air in jar
  mutate(mass_C_P = moles_C_P*molmass_C)  %>%  # [mg] mg of C in air in jar
  mutate(mass_C_V = moles_C_V*molmass_C)  # [mg] mg of C in air in jar

# Removed air inside syringe
Vrem <- 30      # [mL] CO2 rich air removed from jar
Trem <- 25+273  # [K] temp of air removed since in incubator
nrem <- (Pr*Vrem) / (R*Trem)  # [mol] moles of air removed from jar

# CLEAN DATA #######################################################
# Flux
data_all <- data_C %>%
  filter(!Sample %in% c('CO2 FREE', '2008', '2%')) %>% # filters out controls 
  separate(Sample, c("Typ","Num", "Lett"), sep=cumsum(c(1,1,1)), remove = FALSE) %>%  # separates out Palouse/Vershire soil and treatments
  mutate(Date.Time = as.POSIXct(Date.Time, format = '%m/%d/%y %H:%M')) %>% # converts Date.Time from characters to date-time format
  mutate(Date = as.Date(Date)) %>%
  group_by(Sample, Flush) %>%  # group by flush 
  arrange(Date.Time) %>%  # arrange in ascending order
  mutate(time_diff = as.numeric(Date.Time - lag(Date.Time, default = first(Date.Time)), units = 'hours')) %>% # [hours] find time difference in flush groups 
  mutate(mass_diff = ifelse(Typ == 'P', as.numeric(mass_C_P - lag(mass_C_P, default = first(mass_C_P))), 
                            as.numeric(mass_C_V - lag(mass_C_V, default = first(mass_C_V))))) %>% # [mg] find mass_C difference in flush groups for Palouse and Vershire
  mutate(rem_moles_C = C_ppm*nrem/(10^6)) %>%  # [mol] moles of C in removed air
  mutate(rem_mass_C = rem_moles_C*molmass_C)  %>% # [mg] mg of C in removed air
  mutate(adj_mass_diff = as.numeric(ifelse(time_diff != '0', as.numeric(mass_diff + lag(rem_mass_C, default = first(rem_mass_C))), '0'))) %>% # [mg] find adjusted mass by including removal mass_C difference in flush groups
  filter(!time_diff == '0') %>% # delete used values
  mutate(flux = adj_mass_diff/time_diff) %>% #%>% # this depends on prev. line being right)
  filter(Num != '4')  # FILTER OUT CCBP
  #filter_if(~is.numeric(.), all_vars(!is.infinite(.))) # keeps the "last" day of a flux measurement ie. gets rid of the "first" day of each session, that's what we graph

## Save Flux Data
data_flux_P <- data_all %>%
  filter(Typ == 'P') # %>%
  # filter(Num == '1') # just soil

data_flux_V <- data_all %>%
  filter(Typ == 'V')

# EVERY LABEL
# num_labs <- c('Soil Control', 'Corn Stover', 'AD HLFB', 'C-CBP HLFB', 'DASE HLFB')
# names(num_labs) <- c('1', '2', '3', '4','5')

# NEW LABELS - CCBP
num_labs <- c('Soil Control', 'CS1', 'AD1', 'DASE1')
names(num_labs) <- c('1', '2', '3','5')

# just soil 
# num_labs <- c('Soil Control')
# names(num_labs) <- c('1')

# SPECIFY END DATE for the datasheet exports for rcum. CO2 respired values plus the plotting, CHANGE HERE
#end_date <- "2021-12-14 10:27" # !!!! CHANGE THIS TO EXTEND GRAPH !!! "2021-12-14 10:27" corresponds to 43 days
end_date <- "2022-07-28 12:30" # this if the last day for this incubation
end_date135 <- "2022-03-16 17:15" # this is the 135th day date which matches with INC2
lims <- as.POSIXct(strptime(c("2021-11-03 12:45", end_date135), format = "%Y-%m-%d %H:%M")) 
#last_day <- max(data$Day, na.rm = TRUE) # [days] final day of measurement for this datasheet, later when you filter data you can use this to check if it actually filtered it

# Respired
data_resp <- data_all %>% 
  group_by(Sample) %>%
  ungroup(Flush) %>%  # ungroup Flush but keep groups by Sample
  #select(Flush, Sample, Date.Time, flux) %>%  # clean it up
  arrange(Date.Time) %>% # rearrange in ascending order
  mutate(time_days = (Date.Time - lag(Date.Time, k = 1))) %>%     # time difference btwn flux measurements in days
  mutate(time_hours = (Date.Time - lag(Date.Time, k = 1))*24) %>% # time difference btwn flux measurements in hours
  mutate(C_resp = .5*(time_hours)*(flux+lag(flux))) %>% # [mg] trapezoidal area calculation to get C respired
  drop_na(C_resp) %>% # drops rows w/ NAs which arise from the first trapezoid area measurement
  mutate(C_resp_cum = cumsum(as.numeric(C_resp))) %>% # [mg] cumulatively add together trapezoids
  mutate(time = as.numeric(Date.Time - first(Date.Time), units = 'days')) # calculate time difference from first in group in [days]

# write.csv(data_resp, file="respdata.csv", row.names = FALSE)

# Export Resp Data
data_resp_P <- data_resp %>%
  filter(Typ == 'P') %>%
  # filter(Num == '1') %>%
  filter(Date.Time <= end_date135)

stats_resp_P <- data_resp_P %>% # output averages plotted in (lumping together replicates) RESP graphs 
  group_by(Sample, Num) %>%
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  group_by(Num) %>%
  summarise(mean_C_resp_cum = mean(max_C_resp_cum), stdev = sd(max_C_resp_cum)) %>% # [mg]
  mutate(Name = num_labs) %>%
  mutate(Soil = 'P')

stats_resp_P2 <- data_resp_P %>% # output averages plotted in (non-lumped, individual replicates) RESP graphs
  group_by(Typ, Num, Sample) %>%
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  group_by(Typ, Num, Sample) %>%
  summarise(mean_C_resp_cum = mean(max_C_resp_cum)) # [mg]

#print(paste("The incubation period currently spans", last_day, "days!"))
#write.csv(stats_resp_P, file = 'P_summary_cumCresp.csv', row.names = FALSE)

data_resp_V <- data_resp %>%
  filter(Typ == 'V') %>%
  filter(Date.Time <= end_date135)

stats_resp_V <- data_resp_V %>% # output averages plotted in (lumping together replicates) RESP graphs 
  group_by(Sample, Num) %>%
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  group_by(Num) %>%
  summarise(mean_C_resp_cum = mean(max_C_resp_cum), stdev = sd(max_C_resp_cum)) %>% # [mg]
  mutate(Name = num_labs) %>%
  mutate(Soil = 'V')

stats_resp_V2 <- data_resp_V %>% # output averages plotted in (non-lumped, individual replicates) RESP graphs
  group_by(Typ, Num, Sample) %>%
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  group_by(Typ, Num, Sample) %>%
  summarise(mean_C_resp_cum = mean(max_C_resp_cum)) # [mg]

INC1_cumCresp_enddate <- rbind(stats_resp_P, stats_resp_V)
twowayanova_data <- rbind(stats_resp_P2, stats_resp_V2) 

#write.csv(INC1_cumCresp_enddate, file = 'INC1_cumCresp_135enddate.csv', row.names = FALSE)
write.csv(twowayanova_data, file = 'twowayanova_data135.csv', row.names = FALSE)

#print(paste("The incubation period currently spans", last_day, "days!"))
#write.csv(stats_resp_V, file = 'V_summary_cumCresp.csv', row.names = FALSE)

####################################################################
# STATISTICS #######################################################

# 2 WAY ANOVA for INC1, test if treatment and soil type have an effect on mean C resp/fraction of C retained by soil/fraction of C retained by residue by end of incubation

# recode Num to factors
thirteenC_data <- read.csv("INC1_2wayanova_13C.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
twowayanova_data <- merge(twowayanova_data, thirteenC_data, by = 'Sample') 

twowayanova_data$Num <- factor(twowayanova_data$Num,
                               levels = c(1, 2, 3, 4, 5), 
                               labels = c('Soil Control', 'Corn Stover', 'AD HLFB', 'C-CBP HLFB', 'DASE HLFB'))

write.csv(twowayanova_data, file = 'INC1_twowayanova_data.csv', row.names = FALSE)

table(twowayanova_data$Typ, twowayanova_data$Num) # check for balanced design, UNBALANCED SINCE UNQUELA NUM. OF OBSERVATIONS

### MEAN C RESP CUM. ###
ggboxplot(twowayanova_data, x = 'Num', y = 'mean_C_resp_cum', color = 'Typ')  # boxplot shows various treatments and how they compare to each other + soil type

res.aov_meanCresp <- aov(mean_C_resp_cum ~ Typ * Num, data = twowayanova_data) # test interaction btwn Num and Typ
Anova(res.aov_meanCresp, type = 'III') # have to use this bc unbalanced design

# Tukey-Kramer maybe
TukeyHSD(res.aov_meanCresp)  # unclear if this is taking into account unbalanced design, I think it's using Tukey Kramer

### Edit data so it works on fraction data
twowayanova_data2 <- subset(twowayanova_data, Num != 'Soil Control')

### FRACTION OF CARBON REMAINING FROM RESIDUE ###
ggboxplot(twowayanova_data2, x = 'Num', y = 'fr', color = 'Typ')  # boxplot shows various treatments and how they compare to each other + soil type

res.aov_fr <- aov(fr ~ Typ * Num, data = twowayanova_data2) # test interaction btwn Num and Typ
Anova(res.aov_fr, type = 'III') # have to use this bc unbalanced design

# Tukey-Kramer maybe
TukeyHSD(res.aov_fr)  # unclear if this is taking into account unbalanced design, I think it's using Tukey Kramer

### FRACTION OF CARBON REMAINING FROM SOIL ###
ggboxplot(twowayanova_data2, x = 'Num', y = 'fs', color = 'Typ')  # boxplot shows various treatments and how they compare to each other + soil type

res.aov_fs <- aov(fs ~ Typ * Num, data = twowayanova_data2) # test interaction btwn Num and Typ
Anova(res.aov_fs, type = 'III') # have to use this bc unbalanced design

# Tukey-Kramer maybe
TukeyHSD(res.aov_fs)  # unclear if this is taking into account unbalanced design, I think it's using Tukey Kramer

##################################################################
# PLOTTING #######################################################

# Theme and Labels
theme_C <- theme_light() + 
  theme(panel.grid.minor = element_blank(), 
  text = element_text(size = 30), #for facetwrapped plots
  strip.background = element_rect(color="black", fill="#93C5FF", size=1.5, linetype="solid"),
  legend.position = "none",
  plot.title = element_text(hjust = 0.5)
  ) 



# FLUX: Mean and SE Each
p1P<- ggplot(data_flux_P, aes(x=Date.Time, y=flux)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + 
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  theme_C +
  scale_y_continuous(limits=c(0,.25)) +  # sets all plots start at 0 go to .3
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = '135 Day Carbon Flux Evolution in HLFB Amended Palouse Soil') 
p1P

p1V<- ggplot(data_flux_V, aes(x=Date.Time, y=flux)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) +
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  theme_C +
  scale_y_continuous(limits=c(0,.35)) +  # sets all plots start at 0 go to .3
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = '135 Day Carbon Flux Evolution in HLFB Amended Vershire Soil') 
p1V

ggsave("newPflux_135days_scale_mean&se.png", plot = p1P, width = 60, height = 20, units = "cm")  # change this accordingly
ggsave("newVflux_135days_scale_mean&se.png", plot = p1V, width = 60, height = 20, units = "cm")  # change this accordingly


# RESPIRED: Mean and SE Each
p2P <- ggplot(data_resp_P, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + # NON FREE SCALE
  ##facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  ##geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +  # when water was added, comment this out for no lines 
  theme_C +
  scale_y_continuous(limits=c(0,50)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon in 135 Day Incubation of HLFB Amended Palouse Soil') 
p2P

p2V <- ggplot(data_resp_V, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + # NON FREE SCALE
  ##facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  ##geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +  # when water was added, comment this out for no lines 
  theme_C +
  scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon in 135 Day Incubation of HLFB Amended Vershire Soil') 
p2V

ggsave("newPresp_135days_scale_mean&se.png", plot = p2P, width = 60, height = 20, units = "cm")
ggsave("newVresp_135days_scale_mean&se.png", plot = p2V, width = 60, height = 20, units = "cm")


# #lumped figure w/ geom smooth of C respired
theme_lump <- theme_light() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman"), 
        #text = element_text(size = 30), #for facetwrapped plots
        #strip.background = element_rect(color="black", fill="#93C5FF", size=1.5, linetype="solid"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )
# 
# lumped_P <- ggplot(data_resp_P, aes(x=Date.Time, y=C_resp_cum)) +
#   geom_smooth(aes(color = Num), se = TRUE) +
#   scale_x_datetime(limits = lims) +
#   theme_lump +
#   scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
#   scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
#   labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in 267 Day Incubation of HLFB Amended Palouse Soil') 
# lumped_P
# 
# lumped_V <- ggplot(data_resp_V, aes(x=Date.Time, y=C_resp_cum)) +
#   geom_smooth(aes(color = Num), se = TRUE) +
#   scale_x_datetime(limits = lims) +
#   theme_lump +
#   scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
#   scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
#   labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in 267 Day Incubation of HLFB Amended Vershire Soil') 
# lumped_V
# 
# ggsave("Plumped_scale_mean&se.png", plot = lumped_P, width = 60, height = 20, units = "cm")
# ggsave("Vlumped_scale_mean&se.png", plot = lumped_V, width = 60, height = 20, units = "cm")

## C Retention Graphs
# Average C respired per soil type by time
soil_resp <- data_resp %>%
  filter(Num == '1') %>%
  group_by(Typ, Flush) %>%
  summarize(mean_soil_cum = mean(C_resp_cum)) # [mg C] averages out all soil C resp. values per time

# Calculate C retained as percentage of residue C and total treatment C 
datainitC <- read.csv("justinitC.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
data3 <- left_join(data_resp, datainitC, by = 'Sample') 
data3 <- merge(data3, soil_resp, by = c('Typ', 'Flush'))  # matches soil resp. to each measurement at a time point 
data3 <- data3 %>%
  mutate(init_totC = init_C*1000) %>%         # [mg C] initial C (soil+res) in each treatment on average
  mutate(init_resC = init_resC*1000) %>%      # [mg C] initial C (res) in each treatment on average
  mutate(Ctot_ret = 100*(init_totC - C_resp_cum)/init_totC) %>% # [% total C] C retained from total treatment
  mutate(Cres_ret = 100*(init_resC-(C_resp_cum-mean_soil_cum))/init_resC) %>% # [% residue C] C retained from residue in each treatment
  group_by(Num, Typ, Date.Time) %>%
  summarize(meanCtot_ret = mean(Ctot_ret), meanCres_ret = mean(Cres_ret)) %>% # [% total C] average of prev. calculations per treatment on specific days
  ungroup() %>% # necessary to add row after
  add_row(Typ = 'P', Num = c('1', '2', '3', '4', '5'), Date.Time = as.POSIXct('2021-11-01 15:00:00'), meanCres_ret = 100, meanCtot_ret = 100) %>% # add initial anchor point of 100% for all treatments (when incubation began)
  add_row(Typ = 'V', Num = c('1', '2', '3', '4', '5'), Date.Time = as.POSIXct('2021-11-01 15:00:00'), meanCres_ret = 100, meanCtot_ret = 100) 

data4 <- left_join(data_resp, datainitC, by = 'Sample') 
data4 <- data4 %>%
  filter(Typ == 'P') %>%
  filter(Num != '4') %>%
  filter(Num != '1') %>%
  mutate(invC_resp_cum = init_C*1000 - C_resp_cum)  %>%
  mutate(invC_resp_cum_ADJ = case_when(Num == '3' ~ .5*invC_resp_cum,
                                       Num == '2' ~ 1*invC_resp_cum,
                                       Num == '5' ~ .35*invC_resp_cum)) %>%
  mutate(ID = case_when(Num == '2' ~ 'CS1', 
                        Num == '3' ~ 'AD1',
                        Num == '5' ~ 'DASE1'))

write.csv(data4, file = 'INC2_invC_resp.csv', row.names = FALSE)
# OK now you need to export this to the other code and do graph them together

data3P <- data3 %>%
  filter(Typ == 'P') %>%
  mutate(Ctot_Label = round(ifelse(Date.Time == max(Date.Time), meanCtot_ret, NA), 0)) %>%  # add labels to last point of each line
  mutate(Cres_Label = round(ifelse(Date.Time == max(Date.Time), meanCres_ret, NA), 0))

data3V <- data3 %>%
  filter(Typ == 'V') %>%
  mutate(Ctot_Label = round(ifelse(Date.Time == max(Date.Time), meanCtot_ret, NA), 0)) %>%  # add labels to last point of each line
  mutate(Cres_Label = round(ifelse(Date.Time == max(Date.Time), meanCres_ret, NA), 0))

# Graph C retained graphs 
# C retained of only residue graphs
Cret_P <- ggplot(data3P, aes(x=Date.Time, y=meanCres_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(35, 100) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Residue \n [% of Initial Residue C]', title = 'Palouse Soil Incubations') +   # Carbon Retained in Residue in \n 267 Day Incubation of Biofuel Residues in Palouse Soil
  geom_label_repel(aes(label = Cres_Label), min.segment.length = 0, size = 2, force = 2.1, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_P

# C retained of total treatment graphs
Cret_P2 <- ggplot(data3P, aes(x=Date.Time, y=meanCtot_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(60, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment \n [% of Initial Treatment C]', title = 'Palouse Soil Incubations') + # 'Carbon Retained in Treatment in \n 267 Day Incubation of Biofuel Residues in Palouse Soil' 
  geom_label_repel(aes(label = Ctot_Label), min.segment.length = 0 , size = 2, force = .6, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_P2

# C retained of only residue graphs
Cret_V <- ggplot(data3V, aes(x=Date.Time, y=meanCres_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(35, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  #geom_label_repel(aes(), nudge_x = 1, na.rm = TRUE) +  # CHANGE THIS SO THE LABEL WORKS, PAGE IS SAVED IN GOOGLE
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Residue \n [% of Initial Residue C]', title = 'Vershire Soil Incubations') +   # Carbon Retained in Residue in \n 267 Day Incubation of Biofuel Residues in Vershire Soil
  geom_label_repel(aes(label = Cres_Label), min.segment.length = 0, size = 2, force = .5, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_V

# C retained of total treatment graphs
Cret_V2 <- ggplot(data3V, aes(x=Date.Time, y=meanCtot_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') +  
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(60, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment \n [% of Initial Treatment C]', title = 'Vershire Soil Incubations') + # 'Carbon Retained in Treatment in \n 267 Day Incubation of Biofuel Residues in Vershire Soil'
  geom_label_repel(aes(label = Ctot_Label), min.segment.length = 0, size = 2, force = .6, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_V2

ggsave("newCretP_scale_mean&se.png", plot = Cret_P, width = 60, height = 20, units = "cm")
ggsave("newCretV_scale_mean&se.png", plot = Cret_V, width = 60, height = 20, units = "cm")


#Show Initial C
pV_init <- ggplot(data_resp_V, aes(x=Date.Time, y=C_resp_cum)) +
  geom_smooth(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  #stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + # NON FREE SCALE
  #geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  ##geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +  # when water was added, comment this out for no lines 
  theme_C +
  scale_y_continuous(limits=c(0,NA)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon in 267 Day Incubation of HLFB Amended Vershire Soil') 
pV_init



#NOT READY
#################################
# FLUX: Corn Stover Lumped
data_corn_flux <- data_all %>%
  filter(Num != 1)
p_corn_flux <- ggplot(data_corn_flux, aes(x = Date.Time, y = flux, group = Num)) +
  stat_summary(fun.data = "mean_se", aes(group = Num,color = Num)) +
  theme_C +
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux Evolution in Corn Stover \n Treatments Over an 89 Day Incubation') 
p_corn_flux
ggsave("corn_flux.png", plot = p_corn_flux, width = 15, height = 15, units = "cm")

# RESPIRED: Corn Stover Lumped
data_corn_resp <- data_resp %>%
  filter(Num != 1)
p_corn_resp <- ggplot(data_corn_resp, aes(x = Date.Time, y = C_resp_cum, group = Num)) +
  stat_summary(fun.data = "mean_se", aes(group = Num,color = Num)) +
  theme_C +
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in Corn Stover \n Treatments Over an 89 Day Incubation') 
p_corn_resp
ggsave("corn_resp.png", plot = p_corn_resp, width = 15, height = 15, units = "cm")

############################################################################
# SOIL MODELLING #######################################################
# Based off of https://www.bgc-jena.mpg.de/TEE/optimization/2015/12/09/Fractions-Incubations/
# Context from https://escholarship.org/uc/item/9h72f7hk

# Clean data for modelling
data_modP <- data_resp_P %>%
  ungroup(Sample) %>%   # now, not grouped as anything
  select(c('time','Num', 'C_resp_cum'))  %>% # select these columns for ease
  group_by(Num, time) %>%  #  
  summarize(cummCO2 = mean(C_resp_cum)) # sd gives an error for some reason: Stderr = sd(C_resp_cum))   # [mg] amount of carbon respired cumulatively, not in terms of mg C/g soil
  #summarize(cummCO2 = mean(C_resp_cum)/50, Stderr = sd(C_resp_cum/50)) %>%  # /50 so it's in [g C/g soil] since we start w/ ~50g soil, summarizing by all incubations def. loses precision since it's not a rate, it's an absolute amount?, but also it's based off of rate anyways
write.csv(data_mod, file = 'data_modP.csv')

data_modV <- data_resp_V %>%
  ungroup(Sample) %>%   # now, not grouped as anything
  select(c('time','Num', 'C_resp_cum'))  %>% # select these columns for ease
  group_by(Num, time) %>%  #  
  summarize(cummCO2 = mean(C_resp_cum)) # sd gives an error for some reason: Stderr = sd(C_resp_cum))   # [mg] amount of carbon respired cumulatively, not in terms of mg C/g soil
#summarize(cummCO2 = mean(C_resp_cum)/50, Stderr = sd(C_resp_cum/50)) %>%  # /50 so it's in [g C/g soil] since we start w/ ~50g soil, summarizing by all incubations def. loses precision since it's not a rate, it's an absolute amount?, but also it's based off of rate anyways
write.csv(data_mod, file = 'data_modV.csv')

############################################################################
# EXTRA PLOTTING #######################################################

# FLUX: All Plots
plot_all_flux <- ggplot(data = data_all) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Flux [mg/hr]', title = 'Soil Incubation Experiments: C Flux v. Time') +
  facet_wrap(~Sample, nrow = 3) +
  theme_light()
plot_all_flux

# Plot 1 FLUX separately
data_1 <- data_all %>%
  filter(Num == '1')
plot_1_flux <- ggplot(data = data_1) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Carbon Flux v. Time') +
  facet_wrap(~Sample, nrow = 1) +
  theme_light()
plot_1_flux

# Plot 2 and 3 FLUX separately
data_23 <- data_all %>%
  filter(Num != '1')
plot_23_flux <- ggplot(data = data_23) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Carbon Flux v. Time') +
  facet_wrap(~Sample, nrow = 2) +
  theme_light()
plot_23_flux

# Plot LUMPED FLUX v Time [1,2,3 w/ geom_smooth in one image, same scale]
new_labels <- c('1' = 'No Residue Control','2' = '500 um Corn Stover', '3' = '8500 um Corn Stover')
plot_lumped_flux <- ggplot(data_all, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped_flux

# Plot 1 LUMPED FLUX separately
data_lumped1 <- data_all %>%
  filter(Num == '1')
plot_lumped1_flux <- ggplot(data_lumped1, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Cumulative C Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped1_flux

# Plot 2 and 3 LUMPED FLUX separately
data_lumped23 <- data_all %>%
  filter(Num != '1')
plot_lumped23_flux <- ggplot(data_lumped23, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Carbon Flux v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped23_flux

# Plot 1 RESP separately
data_1_resp <- data_resp %>%
  filter(Num == '1')
plot_1_resp <- ggplot(data = data_1_resp) +
  geom_point(aes(x = Date.Time, y = C_resp_cum)) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Incubation of (1) Soil Controls: Cumulative C Respired v. Time') +
  facet_wrap(~Sample, nrow = 1) +
  theme_light()
plot_1_resp

# Plot 2 and 3 RESP separately
data_23_resp <- data_resp %>%
  filter(Num != '1')
plot_23_resp <- ggplot(data = data_23_resp) +
  geom_point(aes(x = Date.Time, y = C_resp_cum)) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Cumulative C Respired v. Time') +
  facet_wrap(~Sample, nrow = 2) +
  theme_light()
plot_23_resp

#Adjusted for controls
ggplot(data_resp, aes(x = Date.Time, y = C_resp_cum, group = Num, col = Num)) +
  geom_point()

# Plot LUMPED RESP. v Time [1,2,3 w/ geom_smooth in one image, same scale]
plot_lumped_resp <- ggplot(data_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Cumulative C Respired v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped_resp

# Plot LUMPED 1 separately
data_lumped1_resp <- data_resp %>%
  filter(Num == '1')
plot_lumped1_resp <- ggplot(data_lumped1_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Cumulative Carbon Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped1_resp

# Plot LUMPED 2 and 3 separately
data_lumped23_resp <- data_resp %>%
  filter(Num != '1')
plot_lumped23_resp <- ggplot(data_lumped23_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Cumulative Carbon Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped23_resp

pacman::p_load(quantmod, ggpmisc, cardidates, SparseM, fANCOVA, signal, foreach, doParallel, grwat, adc, hydroGOF, parallel, doParallel, hydrostats, hydroEvents, tidyr, spatialEco, zoo, npreg, fANCOVA, forecast, FlowScreen, stinepack, imputeTS, missForest, Metrics, dplyr, tstools, tsbox, highfrequency, tidyverse, install = T)

# Load results x=
x1	=	0.5
x2	=	0.9
x3	=	0.0157017579
x4	=	0.5
x5	=	0.105360897
x6	=	0.844400698
x7	=	0.0116280655
x8	=	0.058229537
x9	=	0.346167928
x10	=	0.163222888
x11	=	1
x12	=	0.8
x13	=	0.00875907
x14	=	1
x15	=	2

##################################################################################
## Load
df <- as.data.frame(read.csv("Tule_SC_RB_O_UnT_1.csv"))

## Preallocate
df$Date <- as.Date(df$Date, "%m/%d/%Y")
df$time = 1:nrow(df)
df$strflow = as.numeric(df[,2])

# inspect
sum(is.na(df$strflow))

# assign unaltered
dfr = df

#remove T and pH
dfl = dfr[-c(4:5)]

# impute missing
dfl$strflow = impute.loess(dfr$strflow, s = 0.01, smooth = FALSE)
#dfl$strflow <- ifelse(dfl$strflow <0, 0, df$strflow)
dfl$SC = impute.loess(dfr$SC, s = x13, smooth = FALSE)
#dfl$SC <- ifelse(dfl$SC <0, 0, df$SC)

################################################################################
################################# Eckhardt #####################################
################################################################################

## Step 1: apply Eckhardt RDF and calculate bf and qf
# empty output vector
bf <- rep(NaN, length(dfl$strflow))

# fill in initial value
bf[1] <- dfl$strflow[1]*x2*0.9  # x2 = BFImax

# iterate through remaining values
for (i in 2:length(dfl$strflow)){
  # calculate bf using digital filter
  bf[i] <- (((1-x2)*x1*bf[i-1]) + ((1-x1)*x2*dfl$strflow[i]))/(1-x1*x2) # x7 = recession constant alpha
  
  # make sure 0 <= bf <= strflow
  if (is.na(bf[i]))    bf[i] <- dfl$strflow[i]  # set logical index to prevent crash at BFImax = 1
  if (bf[i]<0|is.na(bf[i]))    bf[i] <- dfl$strflow[i]*x2*0.9  # from Barlow 'Digital Filters' document
  if (bf[i]>dfl$strflow[i]) bf[i] <- dfl$strflow[i]
}

dfl$bf = bf

# calculate qf
dfl$qf = dfl$strflow - bf

################################################################################
## Step 2: calculate BFI and set dry periods to 1
dfl$bfi = dfl$bf/dfl$strflow
dfl$bfi[is.na(dfl$bfi)] <- 1

################################################################################
## Step 3: Get bfi and qfi peaks 
#bfip = bf
bfi_peak_index <- dfl$time[ggpmisc:::find_peaks(dfl$bfi,ignore_threshold = x4)]
bfip <- as.data.frame(dfl$bfi[bfi_peak_index])
bfip$i <- bfi_peak_index
colnames(bfip) <- c('bfi','i')

#qfip = bf
qfi_peak_index <- dfl$time[ggpmisc:::find_peaks(dfl$bfi*-1,ignore_threshold = x5)]
qfip <- as.data.frame(dfl$bfi[qfi_peak_index])
qfip$i <- qfi_peak_index
colnames(qfip) <- c('qfi','i')
  
ggplot() +
  geom_line(data = dfl, aes(x = time, y = bfi), color = "grey", alpha = 0.5, linewidth = 0.5)+
  #geom_line(data = dfl, aes(x = time, y = bfis), color = "blue", alpha = 1, linewidth = 0.5)+
  #geom_vline(xintercept = bfip$i, color = "red", linetype = "dashed", linewidth = 0.3)+
  labs(x = "Time",
       y = "BFI")+
  #xlim(2000, 5000) +
  #ylim(0,50) +
  theme_classic()

################################################################################
################################ SCb bounds ####################################
################################################################################

#SCb Upper Bound (SCy = 0), unrealistic 
SCb_Upper <- as.data.frame(dfl$SC[bfi_peak_index]/bfip$bfi)
SCb_Upper$i <- bfi_peak_index
colnames(SCb_Upper) <- c('SCb', 'i')
SCb_Upper$SCb <- pmax(SCb_Upper$SCb, dfl$SC[bfi_peak_index]) 

#SCb Lower Bound, SCy = SCq
SCb_Lower <- as.data.frame(dfl$SC[bfi_peak_index])
SCb_Lower$i <- bfi_peak_index
colnames(SCb_Lower) <- c('SCb', 'i')
SCb_Lower$SCb <- pmin(SCb_Lower$SCb, dfl$SC[bfi_peak_index],SCb_Upper$SCb)

SCfl <- dfl[c(4)]
#Interpolate SCb lower bound
SCfl$SCb_L = dfl$SC
if (x14 == 1) {
  SCfl$SCb_L = interp1(SCb_Lower$i, SCb_Lower$SCb, dfl$time, "pchip", extrap = NA)
} else if (x14 == 2) {  
  SCfl$SCb_L = interp1(SCb_Lower$i, SCb_Lower$SCb, dfl$time, "linear", extrap = NA)
} else if (x14 == 3) { 
  SCfl$SCb_L = interp1(SCb_Lower$i, SCb_Lower$SCb, dfl$time, "spline", extrap = NA)
}  else { 
  SCfl$SCb_L = interp1(SCb_Lower$i, SCb_Lower$SCb, dfl$time, "cubic", extrap = NA)
}

#Interpolate SCb upper bound
SCfl$SCb_U = dfl$SC
if (x14 == 1) {
  SCfl$SCb_U = interp1(SCb_Upper$i, SCb_Upper$SCb, dfl$time, "pchip", extrap = NA)
} else if (x14 == 2) {  
  SCfl$SCb_U = interp1(SCb_Upper$i, SCb_Upper$SCb, dfl$time, "linear", extrap = NA)
} else if (x14 == 3) { 
  SCfl$SCb_U = interp1(SCb_Upper$i, SCb_Upper$SCb, dfl$time, "spline", extrap = NA)
}  else { 
  SCfl$SCb_U = interp1(SCb_Upper$i, SCb_Upper$SCb, dfl$time, "cubic", extrap = NA)
}

sum(is.na(SCfl$SCb_U[qfi_peak_index]))
sum(is.na(SCfl$SCb_L[qfi_peak_index]))

################################################################################
################################ SCy bounds ####################################
################################################################################

#SCy Upper bound
SCy_Upper <- as.data.frame(((dfl$SC[qfi_peak_index])/(1-dfl$bfi[qfi_peak_index]))
                           -(((SCfl$SCb_L[qfi_peak_index])*((dfl$bf[qfi_peak_index])/dfl$qf[qfi_peak_index]))))
SCy_Upper$i <- qfi_peak_index
colnames(SCy_Upper) <- c('SCy','i')
SCy_Upper[is.na(SCy_Upper)] <- 1
SCy_Upper$SCy <- pmin(SCy_Upper$SCy, dfl$SC[qfi_peak_index])
SCy_Upper[SCy_Upper < 0] <- 0
sum(is.na(SCy_Upper))

#SCy Lower bound
SCy_Lower <- as.data.frame(((dfl$SC[qfi_peak_index])/(1-dfl$bfi[qfi_peak_index]))
                           -(((SCfl$SCb_U[qfi_peak_index])*(dfl$bf[qfi_peak_index]))/dfl$qf[qfi_peak_index]))
SCy_Lower$i <- qfi_peak_index
colnames(SCy_Lower) <- c('SCy','i')
SCy_Lower[is.na(SCy_Lower)] <- 1
SCy_Lower[SCy_Lower < 0] <- 0
SCy_Lower$SCy <- pmin(SCy_Lower$SCy, SCy_Upper$SCy, dfl$SC[qfi_peak_index])
sum(is.na(SCy_Lower))

#Interpolate SCb lower bound
SCfl$SCy_L = dfl$SC
if (x14 == 1) {
  SCfl$SCy_L = interp1(SCy_Lower$i, SCy_Lower$SCy, dfl$time, "pchip", extrap = NA)
} else if (x14 == 2) {  
  SCfl$SCy_L = interp1(SCy_Lower$i, SCy_Lower$SCy, dfl$time, "linear", extrap = NA)
} else if (x14 == 3) { 
  SCfl$SCy_L = interp1(SCy_Lower$i, SCy_Lower$SCy, dfl$time, "spline", extrap = NA)
}  else { 
  SCfl$SCy_L = interp1(SCy_Lower$i, SCy_Lower$SCy, dfl$time, "cubic", extrap = NA)
}

#Interpolate SCb upper bound
SCfl$SCy_U = dfl$SC
if (x14 == 1) {
  SCfl$SCy_U = interp1(SCy_Upper$i, SCy_Upper$SCy, dfl$time, "pchip", extrap = NA)
} else if (x14 == 2) {  
  SCfl$SCy_U = interp1(SCy_Upper$i, SCy_Upper$SCy, dfl$time, "linear", extrap = NA)
} else if (x14 == 3) { 
  SCfl$SCy_U = interp1(SCy_Upper$i, SCy_Upper$SCy, dfl$time, "spline", extrap = NA)
}  else { 
  SCfl$SCy_U = interp1(SCy_Upper$i, SCy_Upper$SCy, dfl$time, "cubic", extrap = NA)
}

################################################################################
############################## Smooth bounds ###################################
################################################################################
x3000 = 0.001 #temporary smoothing parameter SCy

### Smooth SCy Bounds
valid_idx_SCy_Us <- complete.cases(SCfl$SCy_U, dfl$time)
SCfl$SCy_Us <- NA
SCfl$SCy_Us[valid_idx_SCy_Us] <- loess(SCfl$SCy_U[valid_idx_SCy_Us] ~ dfl$time[valid_idx_SCy_Us], span = x3000)$fitted
SCfl$SCy_Us[is.na(SCfl$SCy_Us)] <- 1
sum(is.na(SCfl$SCy_Us))

valid_idx_SCy_Ls <- complete.cases(SCfl$SCy_L, dfl$time)
SCfl$SCy_Ls[valid_idx_SCy_Ls] <- loess(SCfl$SCy_L[valid_idx_SCy_Ls] ~ dfl$time[valid_idx_SCy_Ls], span = x3000)$fitted
SCfl$SCy_Ls[is.na(SCfl$SCy_Ls)] <- 1
sum(is.na(SCfl$SCy_Ls))

x4000 = 0.001 #temporary smoothing parameter SCb

#Smooth SCb bounds
valid_idx_SCb_Us <- complete.cases(SCfl$SCb_U, dfl$time)
SCfl$SCb_Us <- NA
SCfl$SCb_Us[valid_idx_SCb_Us] <- loess(SCfl$SCb_U[valid_idx_SCb_Us] ~ dfl$time[valid_idx_SCb_Us], span = x4000)$fitted
SCfl$SCb_Us[is.na(SCfl$SCb_Us)] <- 1
sum(is.na(SCfl$SCb_Us))

valid_idx_SCb_Ls <- complete.cases(SCfl$SCb_L, dfl$time)
SCfl$SCb_Ls[valid_idx_SCb_Ls] <- loess(SCfl$SCy_L[valid_idx_SCb_Ls] ~ dfl$time[valid_idx_SCb_Ls], span = x4000)$fitted
SCfl$SCb_Ls[is.na(SCfl$SCb_Ls)] <- 1
sum(is.na(SCfl$SCb_Ls))

################################################################################
################################# Graphics #####################################
################################################################################
###Lower Bounds
ggplot() +
  geom_line(data = SCfl, aes(x = time, y = SCb_Ls, color = "SCb Lower (Smoothed)"), alpha = 1, linewidth = 1)+
  #geom_line(data = SCfl, aes(x = time, y = SCb_Us, color = "SCb Upper (Smoothed)"), alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = SC, color = "Raw SC"), alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCy_Us, color = "SCy Upper (Smoothed)"), alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Ls, color = "SCy Lower (Smoothed)"), alpha = 1, linewidth = 0.5)+
  labs(x = "Time",
       y = "SC (mS/cm)",
       color = "Lower Limts")+
  scale_color_manual(
    values = c(
      "SCb Lower (Smoothed)" = "darkblue",
      "Raw SC" = "grey",
      "SCy Upper (Smoothed)" = "orange"
      #,"SCb Upper (Smoothed)" = "darkblue"
    #  ,"SCy Lower (Smoothed)" = "orange"
    )
  ) +
  #xlim(3000, 7000) +
  #ylim(0,50) + 
  theme_classic() + theme(legend.position = "top",legend.title = element_text(size = 10, face = "bold"),legend.text = element_text(size = 9))

################
###Upper Bounds
ggplot() +
  #geom_line(data = SCfl, aes(x = time, y = SCb_Ls, color = "SCb Lower (Smoothed)"), alpha = 1, linewidth = 1)+
  geom_line(data = SCfl, aes(x = time, y = SCb_Us, color = "SCb Upper (Smoothed)"), alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = SC, color = "Raw SC"), alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Us, color = "SCy Upper (Smoothed)"), alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCy_Ls, color = "SCy Lower (Smoothed)"), alpha = 1, linewidth = 0.5)+
  labs(x = "Time",
       y = "SC (mS/cm)",
       color = "Upper Limts")+
  scale_color_manual(
    values = c(
     # "SCb Lower (Smoothed)" = "darkblue",
      "Raw SC" = "grey"
     # ,"SCy Upper (Smoothed)" = "red"
      ,"SCb Upper (Smoothed)" = "darkblue"
        ,"SCy Lower (Smoothed)" = "orange"
    )
  ) +
  #xlim(3000, 7000) +
  #ylim(0,50) + 
  theme_classic() + theme(legend.position = "top",legend.title = element_text(size = 10, face = "bold"),legend.text = element_text(size = 9))



################################################################################
############################ Estimate SCb, SCy #################################
################################################################################

optimize_SC_estimates <- function(w) {
  SCb_est <- SCfl$SCb_Ls + w * (SCfl$SCb_Us - SCfl$SCb_Ls)
  SCy_est <- SCfl$SCy_Ls + w * (SCfl$SCy_Us - SCfl$SCy_Ls)
  SC_recon <- (dfl$bf * SCb_est + dfl$qf * SCy_est) / (dfl$bf + dfl$qf)
  sqrt(mean((dfl$SC - SC_recon)^2, na.rm = TRUE))
}

opt_result <- optimize(optimize_SC_estimates, interval = c(0, 1))
w_best <- opt_result$minimum

SCfl$SCb_est <- SCfl$SCb_Ls + w_best * (SCfl$SCb_Us - SCfl$SCb_Ls)
SCfl$SCy_est <- SCfl$SCy_Ls + w_best * (SCfl$SCy_Us - SCfl$SCy_Ls)

SCfl$SC_recon <- (dfl$bf * SCfl$SCb_est + dfl$qf * SCfl$SCy_est) / (dfl$bf + dfl$qf)

################################################################################
########################### Graphics Estimates #################################
################################################################################
#SCy bounds and SCy_est
ggplot() +
  #geom_line(data = SCfl, aes(x = time, y = SCb_L), color = "darkblue", alpha = 1, linewidth = 1)+
  #geom_line(data = SCfl, aes(x = time, y = SCb_U), color = "darkblue", alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = SC), color = "lightgrey", alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCy_Us), color = "orange", alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCy_Ls), color = "orange", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCb_est), color = "lightblue", alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCy_est), color = "darkred", alpha = 1, linewidth = 0.5)+
  labs(x = "Time",
       y = "SC (mS/cm)")+
  #xlim(3000, 7000) +
  #ylim(0,50) +
  theme_classic()

#SCb bounds and SCb_est
ggplot() +
  geom_line(data = SCfl, aes(x = time, y = SCb_Ls), color = "green", alpha = 1, linewidth = 1)+
  geom_line(data = SCfl, aes(x = time, y = SCb_Us), color = "lightblue", alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = SC), color = "grey", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Us), color = "orange", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Ls), color = "orange", alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCb_est), color = "darkblue", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_est), color = "darkred", alpha = 1, linewidth = 0.5)+
  labs(x = "Time",
       y = "SC (mS/cm)")+
  #xlim(3000, 7000) +
  #ylim(0,50) +
  theme_classic()


#SCb est and SCy est, SC total
ggplot() +
  #geom_line(data = SCfl, aes(x = time, y = SCb_Ls), color = "lightblue", alpha = 1, linewidth = 1)+
  #geom_line(data = SCfl, aes(x = time, y = SCb_Us), color = "lightblue", alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = SC), color = "grey", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Us), color = "orange", alpha = 1, linewidth = 0.5)+
  #geom_line(data = SCfl, aes(x = time, y = SCy_Ls), color = "orange", alpha = 1, linewidth = 0.5)+
  geom_line(data = SCfl, aes(x = time, y = SCb_est), color = "darkblue", alpha = 1, linewidth = 0.2)+
  geom_line(data = SCfl, aes(x = time, y = SCy_est), color = "darkred", alpha = 1, linewidth = 0.2)+
  labs(x = "Time",
       y = "SC (mS/cm)")+
  #xlim(3000, 7000) +
  #ylim(0,50) +
  theme_classic()

ggplot() +
  geom_line(data = dfl, aes(x = time, y = strflow), color = "grey", alpha = 1, linewidth = 1)+
  geom_line(data = dfl, aes(x = time, y = qf), color = "blue", alpha = 1, linewidth = 0.5)+
  geom_line(data = dfl, aes(x = time, y = bf), color = "red", alpha = 1, linewidth = 0.5)+
  #geom_vline(xintercept = SCqfp$i, color = "black", linetype = "dashed", linewidth = 0.3)+
  labs(x = "Time",
       y = "SC (mS/cm)")+
  #xlim(3000, 7000) +
  ylim(0,1000) +
  theme_classic()



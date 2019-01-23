smallflatlist <-function(obj){

  CIDSdata<- obj$raw_data %>%
    mutate(r45o44 = v45.mV / (v44.mV/2))  %>%
    mutate(r46o44 = v46.mV / (v44.mV/2))  %>%
    mutate(r47o44 = v47.mV / (v44.mV/2))  %>%
    mutate(r48o44 = v48.mV / (v44.mV/2))  %>%
    mutate(r49o44 = v49.mV / (v44.mV/2))  %>%

    do({
      ref_pre <- filter(., type == "standard") %>% select(-type) %>% { setNames(., str_c("pre_", names(.))) }
      ref_post <- filter(., type == "standard") %>% select(-type) %>% { setNames(., str_c("post_", names(.))) }
      filter(., type == "sample") %>%
        mutate(pre_ref = cycle - 1, post_ref = cycle) %>%
        left_join(ref_pre, by = c("pre_ref" = "pre_cycle")) %>%
        left_join(ref_post, by = c("post_ref" = "post_cycle"))
    })   %>%
    mutate(r45O44avgstdr45O44times1000 = (r45o44/((pre_r45o44+post_r45o44)/2)-1)*1000) %>%
    mutate(r46O44avgstdr46O44times1000 = (r46o44/((pre_r46o44+post_r46o44)/2)-1)*1000) %>%
    mutate(r47O44avgstdr47O44times1000 = (r47o44/((pre_r47o44+post_r47o44)/2)-1)*1000) %>%
    mutate(r48O44avgstdr48O44times1000 = (r48o44/((pre_r48o44+post_r48o44)/2)-1)*1000) %>%
    mutate(r49O44avgstdr49O44times1000 = (r49o44/((pre_r49o44+post_r49o44)/2)-1)*1000) %>%
    mutate(sample44ostd44times100 = (100*(v44.mV/(v44.mV+v44.mV)/2))) %>%
    mutate(mismatch = v44.mV/(v44.mV+v44.mV)/2) %>%
    mutate(d18O = obj$vendor_data_table$`d 18O/16O`) %>%
    mutate(d13C = obj$vendor_data_table$`d 13C/12C`) %>%
    mutate(d17O = obj$vendor_data_table$`d 17O/16O`)


  meand45 = mean(CIDSdata$r45O44avgstdr45O44times1000)
  meand46 = mean(CIDSdata$r46O44avgstdr46O44times1000)
  meand47 = mean(CIDSdata$r47O44avgstdr47O44times1000)
  meand48 = mean(CIDSdata$r48O44avgstdr48O44times1000)
  meand49 = mean(CIDSdata$r49O44avgstdr49O44times1000)

  sdd45 <- sd(CIDSdata$r45O44avgstdr45O44times1000)
  sdd46 <- sd(CIDSdata$r46O44avgstdr46O44times1000)
  sdd47 <- sd(CIDSdata$r47O44avgstdr47O44times1000)
  sdd48 <- sd(CIDSdata$r48O44avgstdr48O44times1000)
  sdd49 <- sd(CIDSdata$r49O44avgstdr49O44times1000)

  semd45<-sd(CIDSdata$r45O44avgstdr45O44times1000)/sqrt(length(CIDSdata$r45O44avgstdr45O44times1000))
  semd46<-sd(CIDSdata$r46O44avgstdr46O44times1000)/sqrt(length(CIDSdata$r46O44avgstdr46O44times1000))
  semd47<-sd(CIDSdata$r47O44avgstdr47O44times1000)/sqrt(length(CIDSdata$r47O44avgstdr47O44times1000))
  semd48<-sd(CIDSdata$r48O44avgstdr48O44times1000)/sqrt(length(CIDSdata$r48O44avgstdr48O44times1000))
  semd49<-sd(CIDSdata$r49O44avgstdr49O44times1000)/sqrt(length(CIDSdata$r49O44avgstdr49O44times1000))

  R13VPDB <-	0.011237
  R18VSMOW <-	0.002005
  R17VSMOW <-	0.000380
  R47ZeroCO2 <-	4.65908E-05

  #refR13 <- ((-11.103/1000)+1)*R13VPDB #help seb pull out number from file
  refR13 <- ((obj$method_info$standards$delta_value[1]/1000)+1)*R13VPDB
  refR18 <- ((obj$method_info$standards$delta_value[2]/1000)+1)*R18VSMOW
  #refR18 <- ((35.775/1000)+1)*R18VSMOW
  refR17 <- ((refR18/R18VSMOW)^0.5164)*R17VSMOW

  ref12C <- 1/(1+refR13)
  ref13C <- 1-ref12C
  ref16O <- 1/(1+refR18+refR17)
  ref18O <- ref16O*refR18
  ref17O <- ref16O*refR17
  samplemeand13C <- mean(obj$vendor_data_table$`d 13C/12C`)
  samplemeand18O <- mean(obj$vendor_data_table$`d 18O/16O`)
  sampleR13 <- ((samplemeand13C/1000)+1)*R13VPDB
  sampleR18 <- ((samplemeand18O/1000)+1)*R18VSMOW
  sampleR17 <- ((sampleR18/R18VSMOW)^0.5164)*R17VSMOW

  sample12C <- 1/(1+sampleR13)
  sample13C <- 1-sample12C
  sample16O <- 1/(1+sampleR18+sampleR17)
  sample18O <- sample16O*sampleR18
  sample17O <- sample16O*sampleR17


  #blue box AV-AX
  refmass12.16.16 <- ref12C*ref16O*ref16O
  refmass12.16.17 <- ref12C*ref16O*ref17O*2
  refmass13.16.16 <- ref13C*ref16O*ref16O
  refmass12.16.18 <- ref12C*ref16O*ref18O*2
  refmass12.17.17 <- ref12C*ref17O*ref17O
  refmass13.17.16 <- ref13C*ref17O*ref16O*2
  refmass12.17.18 <- ref12C*ref17O*ref18O*2
  refmass13.16.18 <- ref13C*ref16O*ref18O*2
  refmass13.17.17 <- ref13C*ref17O*ref17O
  refmass12.18.18 <- ref12C*ref18O*ref18O
  refmass13.17.18 <- ref13C*ref17O*ref18O*2
  refmass13.18.18 <- ref13C*ref18O*ref18O

  samplemass12.16.16 <- sample12C*sample16O*sample16O
  samplemass12.16.17 <- sample12C*sample16O*sample17O*2
  samplemass13.16.16 <- sample13C*sample16O*sample16O
  samplemass12.16.18 <- sample12C*sample16O*sample18O*2
  samplemass12.17.17 <- sample12C*sample17O*sample17O
  samplemass13.17.16 <- sample13C*sample17O*sample16O*2
  samplemass12.17.18 <- sample12C*sample17O*sample18O*2
  samplemass13.16.18 <- sample13C*sample16O*sample18O*2
  samplemass13.17.17 <- sample13C*sample17O*sample17O
  samplemass12.18.18 <- sample12C*sample18O*sample18O
  samplemass13.17.18 <- sample13C*sample17O*sample18O*2
  samplemass13.18.18 <- sample13C*sample18O*sample18O

  #blue box AZ-BB
  ref44 <- refmass12.16.16
  ref45 <- refmass12.16.17 + refmass13.16.16
  ref46 <- refmass12.16.18 + refmass12.17.17 + refmass13.17.16
  ref47 <- refmass12.17.18 + refmass13.16.18 + refmass13.17.17
  ref48 <- refmass12.18.18 + refmass13.17.18
  ref49 <- refmass13.18.18


  sample44 <- samplemass12.16.16
  sample45 <- samplemass12.16.17 + samplemass13.16.16
  sample46 <- samplemass12.16.18 + samplemass12.17.17 + samplemass13.17.16
  sample47 <- samplemass12.17.18 + samplemass13.16.18 + samplemass13.17.17
  sample48 <- samplemass12.18.18 + samplemass13.17.18
  sample49 <- samplemass13.18.18

  refR45 <- ref45/ref44
  refR46 <- ref46/ref44
  refR47 <- ref47/ref44
  refR48 <- ref48/ref44
  refR49 <- ref49/ref44

  sampleR45 <- sample45/sample44
  sampleR46 <- sample46/sample44
  sampleR47 <- sample47/sample44
  sampleR48 <- sample48/sample44
  sampleR49 <- sample49/sample44

  # gray box
  #top part already in data
  R45 <- ((meand45/1000)+1)*refR45
  R46 <- ((meand46/1000)+1)*refR46
  R47 <- ((meand47/1000)+1)*refR47
  R48 <- ((meand48/1000)+1)*refR48
  R49 <- ((meand49/1000)+1)*refR49

  #Yellow box
  D45 <-((R45/sampleR45)-1)*1000
  D46 <-((R46/sampleR46)-1)*1000
  D47 <-((R47/sampleR47)-1)*1000
  D48 <-((R48/sampleR48)-1)*1000
  D49 <-((R49/sampleR49)-1)*1000

  D47full <- D47-D46-D45
  D48full <- D48-D46-D46
  D49full <- D49-D46-D46-D45

  #back to geeen box

  d47zeroCO2 <- ((R47/R47ZeroCO2)-1)*1000

  PB <- CIDSdata$v44.mV[1] - CIDSdata$pre_v44.mV[1]

  Identifier1 <-obj$file_info$`Identifier 1`
  Analysis <- obj$file_info$Analysis
  Method <- obj$file_info$Method
  Identifier2 <- 1#obj$file_info$`Identifier 2`
  Preparation <- obj$file_info$Preparation
  file_datetime <-obj$file_info$file_datetime
  LeftPressure <- as.integer((sub("\\]","", sub(".*Pressure \\[", "", obj$file_info$measurement_info))[4]))
  RightPressure <- as.integer((sub("\\]","", sub(".*Pressure \\[", "", obj$file_info$measurement_info))[5]))

  lineofsmallflatlist <- data_frame(Analysis,Method,file_datetime,Identifier1,Identifier2,Preparation, meand45,semd45,sdd45,meand46,semd46,sdd46,meand47,semd47,sdd47,meand48,semd48,sdd48,meand49,semd49,sdd49,mean(D47full),mean(D48full),samplemeand13C,samplemeand18O, PB, LeftPressure, RightPressure)

  AcquisitionCorrectedData <- list()
  AcquisitionCorrectedData$lineofsmallflatlist <- lineofsmallflatlist
  AcquisitionCorrectedData$CIDSdata <- CIDSdata
  AcquisitionCorrectedData$CIDSdata <- obj$file_info

  return(AcquisitionCorrectedData)


}

smallflatlist_fixed_ref_gas <-function(obj){

  CIDSdata<- obj$raw_data %>%

    mutate(r45o44 = v45.mV / (v44.mV))  %>%
    mutate(r46o44 = v46.mV / (v44.mV))  %>%
    mutate(r47o44 = v47.mV / (v44.mV))  %>%
    mutate(r48o44 = v48.mV / (v44.mV))  %>%
    mutate(r49o44 = v49.mV / (v44.mV))  %>%

    do({
      ref_pre <- filter(., type == "standard") %>% select(-type) %>% { setNames(., str_c("pre_", names(.))) }
      ref_post <- filter(., type == "standard") %>% select(-type) %>% { setNames(., str_c("post_", names(.))) }
      filter(., type == "sample") %>%
        mutate(pre_ref = cycle - 1, post_ref = cycle) %>%
        left_join(ref_pre, by = c("pre_ref" = "pre_cycle")) %>%
        left_join(ref_post, by = c("post_ref" = "post_cycle"))
    })   %>%
    mutate(d45 = (r45o44/((pre_r45o44+post_r45o44)/2)-1)*1000) %>%
    mutate(d46 = (r46o44/((pre_r46o44+post_r46o44)/2)-1)*1000) %>%
    mutate(d47 = (r47o44/((pre_r47o44+post_r47o44)/2)-1)*1000) %>%
    mutate(d48 = (r48o44/((pre_r48o44+post_r48o44)/2)-1)*1000) %>%
    mutate(d49 = (r49o44/((pre_r49o44+post_r49o44)/2)-1)*1000) %>%
    mutate(r45O44avgstdr45O44times1000 = (d45)) %>%
    mutate(r46O44avgstdr46O44times1000 = (d46)) %>%
    mutate(r47O44avgstdr47O44times1000 = (d47)) %>%
    mutate(r48O44avgstdr48O44times1000 = (d48)) %>%
    mutate(r49O44avgstdr49O44times1000 = (d49)) %>%
    mutate(sample44ostd44times100 = (100*(v44.mV/(v44.mV+v44.mV)/2))) %>%
    mutate(mismatch = v44.mV/(v44.mV+v44.mV)/2)

  CIDSdatanofirstCyc <- filter(CIDSdata, cycle != 7)
  dO17corrected <- correct_CO2_for_17O(CIDSdatanofirstCyc, d45, d46)
  d13C <- c(mean(dO17corrected$d13.raw-11.103),dO17corrected$d13.raw-11.103)
  d18O <- c(mean(dO17corrected$d18.raw+35.775),dO17corrected$d18.raw+35.775)
  CIDSdata <- CIDSdata %>%
    cbind(d18O)%>%
    cbind(d13C)
  CIDSdata <- filter(CIDSdata, cycle != 7)
  meand45 = mean(CIDSdata$r45O44avgstdr45O44times1000)
  meand46 = mean(CIDSdata$r46O44avgstdr46O44times1000)
  meand47 = mean(CIDSdata$r47O44avgstdr47O44times1000)
  meand48 = mean(CIDSdata$r48O44avgstdr48O44times1000)
  meand49 = mean(CIDSdata$r49O44avgstdr49O44times1000)

  sdd45 <- sd(CIDSdata$r45O44avgstdr45O44times1000)
  sdd46 <- sd(CIDSdata$r46O44avgstdr46O44times1000)
  sdd47 <- sd(CIDSdata$r47O44avgstdr47O44times1000)
  sdd48 <- sd(CIDSdata$r48O44avgstdr48O44times1000)
  sdd49 <- sd(CIDSdata$r49O44avgstdr49O44times1000)

  semd45<-sd(CIDSdata$r45O44avgstdr45O44times1000)/sqrt(length(CIDSdata$r45O44avgstdr45O44times1000))
  semd46<-sd(CIDSdata$r46O44avgstdr46O44times1000)/sqrt(length(CIDSdata$r46O44avgstdr46O44times1000))
  semd47<-sd(CIDSdata$r47O44avgstdr47O44times1000)/sqrt(length(CIDSdata$r47O44avgstdr47O44times1000))
  semd48<-sd(CIDSdata$r48O44avgstdr48O44times1000)/sqrt(length(CIDSdata$r48O44avgstdr48O44times1000))
  semd49<-sd(CIDSdata$r49O44avgstdr49O44times1000)/sqrt(length(CIDSdata$r49O44avgstdr49O44times1000))

  R13VPDB <-	0.011237
  R18VSMOW <-	0.002005
  R17VSMOW <-	0.000380
  R47ZeroCO2 <-	4.65908E-05


  #refR13 <- ((-11.103/1000)+1)*R13VPDB #help seb pull out number from file

  refR13 <- ((-11.103/1000)+1)*R13VPDB
  refR18 <- ((35.775/1000)+1)*R18VSMOW
  #refR18 <- ((35.775/1000)+1)*R18VSMOW
  refR17 <- ((refR18/R18VSMOW)^0.5164)*R17VSMOW

  ref12C <- 1/(1+refR13)
  ref13C <- 1-ref12C
  ref16O <- 1/(1+refR18+refR17)
  ref18O <- ref16O*refR18
  ref17O <- ref16O*refR17
  samplemeand13C <- mean(d13C)
  samplemeand18O <- mean(d18O)
  sampleR13 <- ((samplemeand13C/1000)+1)*R13VPDB
  sampleR18 <- ((samplemeand18O/1000)+1)*R18VSMOW
  sampleR17 <- ((sampleR18/R18VSMOW)^0.5164)*R17VSMOW

  sample12C <- 1/(1+sampleR13)
  sample13C <- 1-sample12C
  sample16O <- 1/(1+sampleR18+sampleR17)
  sample18O <- sample16O*sampleR18
  sample17O <- sample16O*sampleR17


  #blue box AV-AX
  refmass12.16.16 <- ref12C*ref16O*ref16O
  refmass12.16.17 <- ref12C*ref16O*ref17O*2
  refmass13.16.16 <- ref13C*ref16O*ref16O
  refmass12.16.18 <- ref12C*ref16O*ref18O*2
  refmass12.17.17 <- ref12C*ref17O*ref17O
  refmass13.17.16 <- ref13C*ref17O*ref16O*2
  refmass12.17.18 <- ref12C*ref17O*ref18O*2
  refmass13.16.18 <- ref13C*ref16O*ref18O*2
  refmass13.17.17 <- ref13C*ref17O*ref17O
  refmass12.18.18 <- ref12C*ref18O*ref18O
  refmass13.17.18 <- ref13C*ref17O*ref18O*2
  refmass13.18.18 <- ref13C*ref18O*ref18O

  samplemass12.16.16 <- sample12C*sample16O*sample16O
  samplemass12.16.17 <- sample12C*sample16O*sample17O*2
  samplemass13.16.16 <- sample13C*sample16O*sample16O
  samplemass12.16.18 <- sample12C*sample16O*sample18O*2
  samplemass12.17.17 <- sample12C*sample17O*sample17O
  samplemass13.17.16 <- sample13C*sample17O*sample16O*2
  samplemass12.17.18 <- sample12C*sample17O*sample18O*2
  samplemass13.16.18 <- sample13C*sample16O*sample18O*2
  samplemass13.17.17 <- sample13C*sample17O*sample17O
  samplemass12.18.18 <- sample12C*sample18O*sample18O
  samplemass13.17.18 <- sample13C*sample17O*sample18O*2
  samplemass13.18.18 <- sample13C*sample18O*sample18O

  #blue box AZ-BB
  ref44 <- refmass12.16.16
  ref45 <- refmass12.16.17 + refmass13.16.16
  ref46 <- refmass12.16.18 + refmass12.17.17 + refmass13.17.16
  ref47 <- refmass12.17.18 + refmass13.16.18 + refmass13.17.17
  ref48 <- refmass12.18.18 + refmass13.17.18
  ref49 <- refmass13.18.18


  sample44 <- samplemass12.16.16
  sample45 <- samplemass12.16.17 + samplemass13.16.16
  sample46 <- samplemass12.16.18 + samplemass12.17.17 + samplemass13.17.16
  sample47 <- samplemass12.17.18 + samplemass13.16.18 + samplemass13.17.17
  sample48 <- samplemass12.18.18 + samplemass13.17.18
  sample49 <- samplemass13.18.18

  refR45 <- ref45/ref44
  refR46 <- ref46/ref44
  refR47 <- ref47/ref44
  refR48 <- ref48/ref44
  refR49 <- ref49/ref44

  sampleR45 <- sample45/sample44
  sampleR46 <- sample46/sample44
  sampleR47 <- sample47/sample44
  sampleR48 <- sample48/sample44
  sampleR49 <- sample49/sample44

  # gray box
  #top part already in data
  R45 <- ((meand45/1000)+1)*refR45
  R46 <- ((meand46/1000)+1)*refR46
  R47 <- ((meand47/1000)+1)*refR47
  R48 <- ((meand48/1000)+1)*refR48
  R49 <- ((meand49/1000)+1)*refR49

  #Yellow box
  D45 <-((R45/sampleR45)-1)*1000
  D46 <-((R46/sampleR46)-1)*1000
  D47 <-((R47/sampleR47)-1)*1000
  D48 <-((R48/sampleR48)-1)*1000
  D49 <-((R49/sampleR49)-1)*1000

  D47full <- D47-D46-D45
  D48full <- D48-D46-D46
  D49full <- D49-D46-D46-D45

  #back to geeen box

  d47zeroCO2 <- ((R47/R47ZeroCO2)-1)*1000

  PB <- CIDSdata$v44.mV[1] - CIDSdata$pre_v44.mV[1]

  Identifier1 <-obj$file_info$`Identifier 1`
  Analysis <- obj$file_info$Analysis
  Method <- obj$file_info$Method
  Identifier2 <- obj$file_info$`Identifier 2`
  Preparation <- obj$file_info$Preparation
  file_datetime <-obj$file_info$file_datetime
  measurement_info_part <- obj$file_info$measurement_info
  Yield <- try(as.numeric(sub("\",.*","", sub(".*left side ","", measurement_info_part))))
  LeftPressure <- as.numeric(sub("  l_p.*","", sub(".*mBar l ", "", measurement_info_part)))
  RightPressure <- as.numeric(sub("  r_p.*","", sub(".*mBar r ", "", measurement_info_part)))
  LeftBellows <- as.numeric(sub("%.*","", sub(".*l_p", "", measurement_info_part)))
  RightBellows <- as.numeric(sub("%.*","", sub(".*r_p ", "", measurement_info_part)))
  v44.mV <- first(obj$raw_data$v44.mV)

  lineofsmallflatlist <- data_frame(Analysis,Method,file_datetime,Identifier1,Identifier2,Preparation, meand45,semd45,sdd45,meand46,semd46,sdd46,meand47,semd47,sdd47,meand48,semd48,sdd48,meand49,semd49,sdd49,mean(D47full),mean(D48full),samplemeand13C,samplemeand18O, PB,Yield, LeftPressure, RightPressure,LeftBellows,RightBellows, v44.mV)

  AcquisitionCorrectedData <- list()
  AcquisitionCorrectedData$lineofsmallflatlist <- lineofsmallflatlist
  AcquisitionCorrectedData$CIDSdata <- CIDSdata
  AcquisitionCorrectedData$CIDSdata <- obj$file_info

  return(AcquisitionCorrectedData)
}


#' Correct CO2 data for O17
#'
#' This function takes d45 and d46 and corrects the values for O17 to derive raw d15 and d18 values.
#' The equations used are based on those derived by Jan Kaiser and Thomas Röckmann in
#' "Correction of mass spectrometric isotope ratio measurements for isobaric isotopologues of O2, CO, CO2, N2O and SO2" (Rapid Communications in Mass Spectrometry, 2008, 3997--4008).
#' It is important to note that this function does not currently take site preference into consideration
#'
#' @param data the data frame
#' @param d45 column (in permil)
#' @param d46 column (in permil)
#' @param ref_17R the 17O/16O reference ratio, VSMOW by default (value from #from Brand et al, 2010)
#' @param ref_18R the 18O/16O reference ratio, VSMOW by default (value from #from Brand et al, 2010)
#'
#' @param lambda the mass dependent scaling coefficient for the oxygen isotopes, default 0.528 (value from from Brand et al, 2010)
#' @param d_max the maximum +/- delta value to consider in the root finding [in permil], should not need to change this unless samples are heavily enriched
#' @param quiet whether the function should output information messages or be quiet (default is to output)
#' @export
#' @return the data frame with corrected 17O
correct_CO2_for_17O <- function (data, d45, d46, ref_17R = 0.000393, ref_13R=0.011180, ref_18R = 0.00208835, lambda = 0.528, d_max = 1000) {

  if (missing(d45)) stop("please specify the column that holds the d45 data")
  if (missing(d46)) stop("please specify the columm that holds the d46 data")
  d45_quo <- enquo(d45)
  d46_quo <- enquo(d46)

  # fitting parameters
  b <- ref_17R^2 / 2*ref_18R
  c <- ref_17R/ref_13R       #17Rr/13Rr
  d <- b/c

  #' expects raw delta values (i.e. NOT in permil)
  d18_root_eq <- function(d18, d45, d46) {
    d17 <- md_scale_delta(d18, lambda, unit = 1)
    return( # note: comment out terms here to see the effect
      d18 - (d46 + d*((2+c)*d46-(2-2*c)*d17-(2+4*c)*d45*(1+d17)+3*c*d17^2)))
  }

  #' expects raw delta values (i.e. NOT in permil)
  calc_d18 <- function(d45, d46) {
      na_idx <- is.na(d45) | is.na(d46)
      out <- rep(NA_real_, length(d45))
      out[!na_idx] <- 
        mapply(function(.d45, .d46) {
          uniroot(function(x) d18_root_eq(x, .d45, .d46), lower = -d_max/1000, upper = d_max/1000, tol = 1e-9)$root
        }, d45[!na_idx], d46[!na_idx])
      return(out)
    }

  calc_d13 <- function(d18, d45) {
    d17 <- md_scale_delta(d18, lambda, unit = 1)
    return(d45 + 2*c * (d45 - d17))
  }

  df <- data %>%
    mutate(
      .d45 = !!d45_quo,
      .d46 = !!d46_quo,
      .d18 = 1000 * calc_d18 (.d45/1000, .d46/1000),
      p.17Ocor = paste0("scaling=", lambda, "; ref 17R=", ref_17R, ", 18R=", ref_18R, ", 13C="),
      d13.raw = 1000 * calc_d13 (.d18/1000, .d45/1000),
      d18.raw = .d18
    )


  df %>% select(-.d45, -.d46, -.d18) %>% return()
}

#---------- isotope convenience functions (not exported) -----------

# multiplication factor for permil
PERMIL = 1000;

# Calculate the 45R from N2O
calculate_45R <- function (`15R`, `17R`) {
  return (2 * `15R` + `17R`)
}

# Calculate the 46R from N2O
calculate_46R <- function (`15R`, `17R`, `18R`) {
  return (`15R`^2 + `18R` + 2 * `15R` * `17R`)
}

# Convert delta to ratio
delta_to_ratio <- function (delta, ref_ratio, unit = PERMIL) {
  return ((delta/PERMIL + 1) * ref_ratio)
}

# Convert ratio to delta
ratio_to_delta <- function (ratio, ref_ratio, unit = PERMIL) {
  return ( (ratio/ref_ratio - 1) * unit )
}

# Mass scale delta
md_scale_delta <- function(delta, lambda, unit = PERMIL) {
  ((delta/unit + 1) ^ lambda - 1)*unit
}


pick_mi <- function(mi, pattern) {
  x <- str_subset(mi, pattern)
  if (length(x) == 0) NA_character_ else x[1]
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#' @param rawdata the data frame
#' @param quiet whether the function should output information messages or be quiet (default is to output)
#' @export
#' @examples
#' @return the data frame with pre and post cyc as well as clumped isotope information
# parse measurement info
clumpedbyCyc <- function (rawdata, lambda = 0.528){
  raw_data_w_measurement_info <-
    rawdata %>%
    # nest to run the operations on the measurement info only once for each file
    nest(-file_id) %>%
    # find the relevant measurement info
    mutate(
      # pull out the measurement info from the nested data (the same for each raw data row so only need to look at 1)
      mi = map(data, ~.x$measurement_info[[1]]),
      # pick the relevant entries in the measurement info
      mi_select = map(mi,
                      ~data_frame(
                        Yield = pick_mi(.x, "left side"),
                        lp = pick_mi(.x, "l_p"),
                        rp = pick_mi(.x, "r_p"),
                        pc = pick_mi(.x, "PC"),
                        bgrd = pick_mi(.x, "Background")
                      )
      )
    ) %>%
    # unnest the relevant measurement info
    unnest(mi_select) %>%
    # extract combined values into multiple columns
    extract(lp, into = c("LeftPressure", "LeftBellows"),
            regex = "mBar l ([0-9.]+)   l_p ([0-9.]+)") %>%
    extract(rp, into = c("RightPressure", "RightBellows"),
            regex = "mBar r ([0-9.]+)   r_p ([0-9.]+)") %>%
    extract(pc, into = c("PC"),
            regex = "PC \\[([0-9.-]+)") %>%
    extract(bgrd, into = str_c("v", c("44", "45", "46", "47", "47.5", "48", "49"), ".background"),
            regex = str_c(rep("([0-9.-]+) mV,?", 7), collapse = "")) %>%
    # turn into numbers
    mutate_at(
      vars(Yield, LeftPressure, LeftBellows, RightPressure, RightBellows, PC, ends_with("background")),
      funs(parse_number)
    ) %>%
    # discard unnecessary columns
    select(-mi) %>%
    # unnest all data again
    unnest(data)

  isostandards <- did_files %>% iso_get_standards_info() %>% select(file_id, delta_name, delta_value) %>% mutate(delta_name = str_c("ref ", delta_name)) %>% spread(delta_name, delta_value)

  ref_pre <- filter(raw_data_w_measurement_info, type == "standard") %>%
    select(-type, -Analysis) %>%
    # prefix column names with pre
    { setNames(., str_c("pre_", names(.))) }

  ref_post <- filter(raw_data_w_measurement_info, type == "standard") %>%
    select(-type, -Analysis) %>%
    # prefix column names with post
    { setNames(., str_c("post_", names(.))) }


  # now combine the samples with their respective pre and post cycle standards (i.e. to bracket them)
  combined_data <-
    raw_data_w_measurement_info %>%
    filter(., type == "sample") %>%
    mutate(pre_ref = cycle - 1, post_ref = cycle) %>%
    left_join(ref_pre, by = c("file_id" = "pre_file_id", "pre_ref" = "pre_cycle")) %>%
    left_join(ref_post, by = c("file_id" = "post_file_id", "post_ref" = "post_cycle")) %>%
    mutate(
      Analysis = parse_number(str_c(Analysis, ".", cycle)),
      `var V 54.mV` = v54.mV,
      `var V 44.mV` = v44.mV,
      `var R sample 45CO2/44CO2` = (v45.mV/v44.mV),
      `var R standard 45CO2/44CO2` = (pre_v45.mV/pre_v44.mV + post_v45.mV/post_v44.mV)/2,
      `var d 45CO2/44CO2` = (`var R sample 45CO2/44CO2` / `var R standard 45CO2/44CO2` - 1) * 1000,
      `r45o44` = v45.mV / v44.mV,
      `r46o44` = v46.mV / v44.mV,
      `r47o44` = v47.mV / v44.mV,
      `r48o44` = v48.mV / v44.mV,
      `r49o44` = v49.mV / v44.mV,
      `pre_r45o44` = pre_v45.mV / pre_v44.mV,
      `pre_r46o44` = pre_v46.mV / pre_v44.mV,
      `pre_r47o44` = pre_v47.mV / pre_v44.mV,
      `pre_r48o44` = pre_v48.mV / pre_v44.mV,
      `pre_r49o44` = pre_v49.mV / pre_v44.mV,
      `post_r45o44` = post_v45.mV / post_v44.mV,
      `post_r46o44` = post_v46.mV / post_v44.mV,
      `post_r47o44` = post_v47.mV / post_v44.mV,
      `post_r48o44` = post_v48.mV / post_v44.mV,
      `post_r49o44` = post_v49.mV / post_v44.mV,
      `d45` = ((r45o44 /((pre_r45o44+post_r45o44)/2)-1)*1000),
      `d46` = ((r46o44 /((pre_r46o44+post_r46o44)/2)-1)*1000),
      `d47` = ((r47o44 /((pre_r47o44+post_r47o44)/2)-1)*1000),
      `d48` = ((r48o44/((pre_r48o44+post_r48o44)/2)-1)*1000),
      `d49` = ((r49o44/((pre_r49o44+post_r49o44)/2)-1)*1000),
      PB = v44.mV - (pre_v44.mV+post_v44.mV)/2
    )
  combined_data <- left_join(combined_data,isostandards, "file_id")
  combined_data <- correct_CO2_for_17O(combined_data,d45, d46)
  combined_data <- mutate(combined_data, d13C = d13.raw+`ref d 13C/12C`)#11.103
  combined_data <- mutate(combined_data, d18O = d18.raw+`ref d 18O/16O`)#35.775

  R13VPDB <-	0.011237
  R18VSMOW <-	0.002005
  R17VSMOW <-	0.000380
  R47ZeroCO2 <-	4.65908E-05

  combined_data <- combined_data %>% mutate(
    refR13 = ((`ref d 13C/12C`/1000)+1)*R13VPDB,
    refR18 = ((`ref d 18O/16O`/1000)+1)*R18VSMOW,
    refR17 = ((refR18/R18VSMOW)^lambda)*R17VSMOW,
    ref12C = 1/(1+refR13),
    ref13C = 1-ref12C,
    ref16O = 1/(1+refR18+refR17),
    ref18O = ref16O*refR18,
    ref17O = ref16O*refR17,
    sampleR13 = ((d13C/1000)+1)*R13VPDB,
    sampleR18 = ((d18O/1000)+1)*R18VSMOW,
    sampleR17 = ((sampleR18/R18VSMOW)^lambda)*R17VSMOW,
    sample12C = 1/(1+sampleR13),
    sample13C = 1-sample12C,
    sample16O = 1/(1+sampleR18+sampleR17),
    sample18O = sample16O*sampleR18,
    sample17O = sample16O*sampleR17,
    #blue box AV-AX
    refmass12.16.16 = ref12C*ref16O*ref16O,
    refmass12.16.17 = ref12C*ref16O*ref17O*2,
    refmass13.16.16 = ref13C*ref16O*ref16O,
    refmass12.16.18 = ref12C*ref16O*ref18O*2,
    refmass12.17.17 = ref12C*ref17O*ref17O,
    refmass13.17.16 = ref13C*ref17O*ref16O*2,
    refmass12.17.18 = ref12C*ref17O*ref18O*2,
    refmass13.16.18 = ref13C*ref16O*ref18O*2,
    refmass13.17.17 = ref13C*ref17O*ref17O,
    refmass12.18.18 = ref12C*ref18O*ref18O,
    refmass13.17.18 = ref13C*ref17O*ref18O*2,
    refmass13.18.18 = ref13C*ref18O*ref18O,

    samplemass12.16.16 = sample12C*sample16O*sample16O,
    samplemass12.16.17 = sample12C*sample16O*sample17O*2,
    samplemass13.16.16 = sample13C*sample16O*sample16O,
    samplemass12.16.18 = sample12C*sample16O*sample18O*2,
    samplemass12.17.17 = sample12C*sample17O*sample17O,
    samplemass13.17.16 = sample13C*sample17O*sample16O*2,
    samplemass12.17.18 = sample12C*sample17O*sample18O*2,
    samplemass13.16.18 = sample13C*sample16O*sample18O*2,
    samplemass13.17.17 = sample13C*sample17O*sample17O,
    samplemass12.18.18 = sample12C*sample18O*sample18O,
    samplemass13.17.18 = sample13C*sample17O*sample18O*2,
    samplemass13.18.18 = sample13C*sample18O*sample18O,
    #blue box AZ-BB
    ref44 = refmass12.16.16,
    ref45 = refmass12.16.17 + refmass13.16.16,
    ref46 = refmass12.16.18 + refmass12.17.17 + refmass13.17.16,
    ref47 = refmass12.17.18 + refmass13.16.18 + refmass13.17.17,
    ref48 = refmass12.18.18 + refmass13.17.18,
    ref49 = refmass13.18.18,

    sample44 = samplemass12.16.16,
    sample45 = samplemass12.16.17 + samplemass13.16.16,
    sample46 = samplemass12.16.18 + samplemass12.17.17 + samplemass13.17.16,
    sample47 = samplemass12.17.18 + samplemass13.16.18 + samplemass13.17.17,
    sample48 = samplemass12.18.18 + samplemass13.17.18,
    sample49 = samplemass13.18.18,

    refR45 = ref45/ref44,
    refR46 = ref46/ref44,
    refR47 = ref47/ref44,
    refR48 = ref48/ref44,
    refR49 = ref49/ref44,

    sampleR45 = sample45/sample44,
    sampleR46 = sample46/sample44,
    sampleR47 = sample47/sample44,
    sampleR48 = sample48/sample44,
    sampleR49 = sample49/sample44,

    # # gray box
    R45 = ((d45/1000)+1)*refR45,
    R46 = ((d46/1000)+1)*refR46,
    R47 = ((d47/1000)+1)*refR47,
    R48 = ((d48/1000)+1)*refR48,
    R49 = ((d49/1000)+1)*refR49,

    # #Yellow box
    D45 = ((R45/sampleR45)-1)*1000,
    D46 = ((R46/sampleR46)-1)*1000,
    D47 = ((R47/sampleR47)-1)*1000,
    D48 = ((R48/sampleR48)-1)*1000,
    D49 = ((R49/sampleR49)-1)*1000,

    D47full = D47-D46-D45,
    D48full = D48-D46-D46,
    D49full = D49-D46-D46-D45,
    runinfo ="",
    Donotuse = FALSE

  )
}

#--------------------------------------------------------------------------------------------------
#' @param flatlist.Cyc the data frame normaly from clumpedbyCyc
#' @param quiet whether the function should output information messages or be quiet (default is to output)
#' @examples
#' @export
#' @return the data frame with acquisition level clumped isotope information
# parse measurement info

clumpedCyctoAcquisition<-  function(flatlist.Cyc) {
  flatlist.Cyc %>%
    group_by(file_id) %>%
    summarise(
      Analysis = floor(Analysis[1]),
      Identifier1 =`Identifier 1`[1],
      Identifier2 = `Identifier 2`[1],
      file_datetime =file_datetime[1],
      mass=as.numeric(`Identifier 2`[1]),
      Preparation = Preparation[1],
      runinfo = runinfo [1],
      Donotuse = Donotuse[1],
      Yield = Yield[1],
      Method = Method[1],
      d45.stdev= sd(d45),
      d45 = mean(d45),
      d46.stdev = sd(d46),
      d46 = mean(d46),
      d47.stdev= sd(d47),
      d47 = mean(d47),
      d48.stdev = sd(d48),
      d48 = mean(d48),
      d49.stdev = sd(d49),
      d49 = mean(d49),
      D47.stdev = sd(D47full),
      D47full= mean(D47full),
      D48.stdev = sd(D48full),
      D48full= mean(D48full),
      D49.stdev = sd(D49full),
      D49full= mean(D49full),
      d13C.stdev= sd(d13C),
      d13C =  mean(d13C),
      d18O.stdev= sd(d18O),
      d18O.VPDB.min = ((((mean(d18O)-30.86)/1.03086)+1000)/1.00821)-1000,
       #splitting VSMOW to VPDB convertion then acid fractionation correction 
      d18O = mean(d18O),
      d18O.ref = `ref d 18O/16O`[1],
      d13C.ref = `ref d 13C/12C`[1],
      PB= mean(PB),
      LeftPressure= mean(LeftPressure),
      RightPressure = mean(LeftPressure)
    )
}

#--------------------------------------------------------------------------------------------------
#' @param combined_data_acquisition the data frame normaly from clumpedCyctoAcquisition
#' @param quiet whether the function should output information messages or be quiet (default is to output)
#' @examples
#' @export
#' @return the data frame with grouped and averaged acquisition level clumped isotope information
# parse measurement info

clumpedAcquisitiontoRun<-  function(combined_data_acquisition) {
  combined_data_acquisition %>% mutate(
    new_sample =  Preparation != c("", head(Preparation, -1))  | `Identifier1` != c("", head(`Identifier1`,-1)), #
    batch = cumsum(new_sample)
  ) %>%
    group_by(batch) %>% mutate(id = row_number()) %>%  add_tally()
}

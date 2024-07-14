###   error operations ########################################

# function to propagate standard errors while dividing or multiplying values
# arguments are the two means and their SEs
errordiv = function(x1,x2,se1,se2)
{
  r = x1/x2
  t = data.frame(se1/x1,se2/x2)
  ser = r*sqrt(t[,1]^2 + t[,2]^2)
  a = data.frame(freq = numeric(length(r)))
  a$freq = r
  a$se = ser
  return(a)
}

erroradd = function(vec)
{
  err = sqrt(sum(vec^2))
  return(err)
}

# function to take in two means and SEs, simulate 1000 values from each and then get 1000 ratios
simerrordiv = function(x1, x2, se1, se2)
{
  # takes untransformed (link) mean and SE values, and generates normal dist. from 1000 sims,
  # then transformed 
  # after the function, lower and upper quantiles are selected as limits of 95% CI
  tp = data.frame(num = clogloglink(rnorm(1000, x1, se1), inverse = T), 
                  den = clogloglink(rnorm(1000, x2, se2), inverse = T)) %>%
    reframe(rat = num/den, 
            val = num)
  
  return(tp)
}


###   create a set of locations ########################################

# function to select a random GROUP.ID from each location in the dataset 
# (spatial subsampling) so that each set of IDs can form the basis for a dataset to use for analysis
# this is then repeated 1000 times in the script "create_random_groupids.R" 
# to create 1000 datasets, each without pseudoreplication

# locs is a data frame with location, group id info

createrandomlocs = function(locs)
{
  require(tidyverse)
  
  locs1 = locs %>% 
    group_by(LOCALITY.ID, month) %>% sample_n(1)
  
  return(locs1$group.id)
}






### expandbyspecies ########################################

# this function adds absences (in addition to presences) to subsets of the data that
# need to be analyzed - for the analysis of a certain species, every complete list 
# without the species is expanded to have 0s added against that species name

# ensure that the working directory has list of India's birds with scientific names 
# (just a safety mechanism for the function to work for small subsets, needs to be enabled if required)
# only need to input data, the species of interest and the complete list of India's bird species
# also groupspecs if required (a dataframe with all relevant list level info), it is defaulted to data

expandbyspecies = function(data, species)
{
  require(tidyverse)
  
  data <- data %>% 
    mutate(across(contains("gridg"), ~ as.factor(.))) %>% 
    mutate(timegroups = as.factor(timegroups))

  # considers only complete lists
  
  checklistinfo = data %>%
    distinct(gridg1, gridg2, gridg3, gridg4, 
             ALL.SPECIES.REPORTED, OBSERVER.ID, 
             #city,
             #DURATION.MINUTES,EFFORT.DISTANCE.KM,
             group.id, month, year, no.sp, timegroups) %>%
    filter(ALL.SPECIES.REPORTED == 1) %>%
    distinct(group.id, .keep_all = TRUE)
  
  # expand data frame to include the bird species in every list
  expanded = checklistinfo %>% 
    mutate(COMMON.NAME = species) %>% 
    left_join(data) %>%
    dplyr::select(-c("COMMON.NAME","gridg2","gridg4","OBSERVER.ID",
                     "ALL.SPECIES.REPORTED","group.id","year","gridg0")) %>% 
  # deal with NAs (column is character)
  mutate(OBSERVATION.COUNT = case_when(is.na(OBSERVATION.COUNT) ~ 0,
                                       OBSERVATION.COUNT != "0" ~ 1, 
                                       TRUE ~ as.numeric(OBSERVATION.COUNT)))

  return(expanded)
}

# faster version of this function using data.table and dtplyr
# optimising runtime
# previous expandbyspecies() to be retired entirely in next annual update
expand_dt = function(data, species) {

  require(tidyverse)
  require(dtplyr)
  require(data.table)
  

  setDT(data)
  data <- data %>% 
    lazy_dt(immutable = FALSE) |> 
    mutate(across(contains("gridg"), ~ as.factor(.))) |> 
    as.data.table()


  # Get distinct rows and filter based on a condition
  # (using base data.table because lazy_dt with immutable == FALSE would
  # modify data even though we are assigning to checklistinfo.
  # and immutable == TRUE copies the data and this is a huge bottleneck)
  # considers only complete lists
  checklistinfo <- unique(data[, 
      .(gridg1, gridg2, gridg3, gridg4, ALL.SPECIES.REPORTED, OBSERVER.ID, 
        group.id, month, year, no.sp)
      ])[
        # filter
      ALL.SPECIES.REPORTED == 1
    ]
  
  checklistinfo <- checklistinfo[
    , 
    .SD[1], # subset of data
    by = group.id
]
  
    
  # expand data frame to include the bird species in every list
  data2 = checklistinfo %>% 
    lazy_dt(immutable = FALSE) |> 
    mutate(COMMON.NAME = species) %>% 
    left_join(data |> lazy_dt(immutable = FALSE),
              by = c("group.id", "gridg1", "gridg2", "gridg3", "gridg4",
                      "ALL.SPECIES.REPORTED", "OBSERVER.ID", "month", "year", 
                      "no.sp","COMMON.NAME")) %>%
    dplyr::select(-c("COMMON.NAME","gridg2","gridg4","OBSERVER.ID",
                     "ALL.SPECIES.REPORTED","group.id","year","gridg0")) %>% 
    # deal with NAs (column is character)
    mutate(OBSERVATION.COUNT = case_when(is.na(OBSERVATION.COUNT) ~ 0,
                                       OBSERVATION.COUNT != "0" ~ 1, 
                                       TRUE ~ as.numeric(OBSERVATION.COUNT))) |> 
    as_tibble()
  
  return(data2)

}





### run models ########################################

# trends
singlespeciesrun = function(data, species, specieslist, restrictedspecieslist)
{
  require(tidyverse)
  require(merTools)
  
  data1 = data
  
  # get information for the species of interest 
  specieslist2 = specieslist %>% filter(COMMON.NAME == species)
  
  # three different flags for three different model types that will be run.
  # 0 is normal model, with full random effects. depending on restricted species,
  # model changes slightly.
  flag = 0
  if (species %in% restrictedspecieslist$COMMON.NAME)
  {
    flag = 1
    restrictedlist1 = restrictedspecieslist %>% filter(COMMON.NAME == species)
    specieslist2$ht = restrictedlist1$ht
    specieslist2$rt = restrictedlist1$rt
    
    if (restrictedlist1$mixed == 0) {
      flag = 2
    }
  }
  
  data1 = data1 %>%
    filter(COMMON.NAME == species) %>%
    distinct(gridg3, month) %>% 
    left_join(data1)
  
  datay = data1 %>%
    distinct(gridg3, gridg1, group.id, .keep_all = TRUE) %>% 
    group_by(gridg3, gridg1) %>% 
    reframe(medianlla = median(no.sp)) %>%
    group_by(gridg3) %>% 
    reframe(medianlla = mean(medianlla)) %>%
    reframe(medianlla = round(mean(medianlla)))
  
  medianlla = datay$medianlla
  
  
  # expand dataframe to include absences as well
  ed = expand_dt(data1, species) %>% 
    # converting months to seasons
    mutate(month = as.numeric(month)) %>% 
    mutate(month = case_when(month %in% c(12,1,2) ~ "Win",
                             month %in% c(3,4,5) ~ "Sum",
                             month %in% c(6,7,8) ~ "Mon",
                             month %in% c(9,10,11) ~ "Aut")) %>% 
    mutate(month = as.factor(month))


  # the model ---------------------------------------------------------------
  
  if (flag == 0)
  {
    m1 = glmer(OBSERVATION.COUNT ~ month + month:log(no.sp) + (1|gridg3/gridg1), 
               data = ed, family = binomial(link = 'cloglog'), 
               nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
  }
  
  if (flag == 1)
  {
    m1 = glmer(OBSERVATION.COUNT ~ month + month:log(no.sp) + (1|gridg1), 
               data = ed, family = binomial(link = 'cloglog'), 
               nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
  }
  
  if (flag == 2)
  {
    m1 = glm(OBSERVATION.COUNT ~ month + month:log(no.sp), 
             data = ed, family = binomial(link = 'cloglog'))
  }
  

  # predicting from model ---------------------------------------------------

  # prepare a new data file to predict
  ltemp <- ed %>% 
    group_by(month) %>% 
    mutate(no.sp = medianlla,
           # taking the first value but any random value will do because we do not
           # intend to prect random variation across grids
           gridg1 = data1$gridg1[1], 
           gridg3 = data1$gridg3[1])

  f2 <- ltemp %>% 
    mutate(freq = 0, se = 0) %>%
    dplyr::select(freq, se)
  
  
  if (flag != 2)
  {
    #pred = predict(m1, newdata = ltemp, type = "response", re.form = NA, allow.new.levels=TRUE)
    pred = predictInterval(m1, newdata = ltemp, which = "fixed",
                           level = 0.48, type = "linear.prediction")
    f2$freqt = pred$fit
    f2$set = pred$fit-pred$lwr
  }
  
  if (flag == 2)
  {
    pred = predict(m1, newdata = ltemp, type = "link", se.fit = T)
    f2$freqt = pred$fit
    f2$set = pred$se.fit
  }
  
  f1 = f2 %>%
    filter(!is.na(freqt) & !is.na(se)) %>%
    # average across month
    reframe(freq = mean(freqt), se = mean(set))
  
  
  tocomb = c(species, f1$freq, f1$se)
  return(tocomb)
  # each species's tocomb becomes one column in final trends0 output object
  
}

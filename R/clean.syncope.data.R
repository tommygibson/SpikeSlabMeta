#### Clean up syncope data

library(here)
library(tidyverse)

univOR <- read.csv(here("R", "syncope_data.csv"))

# get rid of obs with no info
# get rid of obervation(s) from Oh's paper -- weird estimate for ECG
# get rid of Colivicci, del rosso, martin as well
# Derose was a subset of Gabayan, and those were the only two papers to contribute to the meta analysis for myocardial infarction
# Also, only one study has incontinence, so remove that observation
# We're taking ECG data from quinn's 2011 paper, which uses the same data as his 2004 paper
# so delete ECG info from the 2004 paper
univOR <- filter(univOR, 
                 event_exposure > 0 | Orhat > 0,
                 !(Paper %in% c("Oh", "Colivicchi", "Del Rosso", "Martin")),
                 !(Variable %in% c("Myocardial Infarction", "Incontinence")),
                 !(Variable == "ECG" & Paper == "Quinn"),
                 !(Variable == "Hemoglobin" & Paper == "Thiruganasambandamoorthy"))

#fill in contingency table counts
univOR <- univOR %>% mutate(
  event_noexposure = event_total - event_exposure,
  noevent_noexposure = nonevent_total - noevent_exposure,
  ORhat1 = (event_exposure * noevent_noexposure) / (event_noexposure * noevent_exposure),
  ORhat2 = ifelse(is.na(Orhat), ORhat1, Orhat),
  lnORhat = log(ORhat2),
  ln.lower = log(OR.lower),
  ln.upper = log(OR.upper),
  SE.chi = sqrt(lnORhat ^ 2 / Chisquare),
  SE.counts = sqrt(1 / event_noexposure + 1 / event_exposure + 1 / noevent_exposure + 1 / noevent_noexposure),
  SE.extrap = ((lnORhat - ln.lower) / 1.96 + (ln.upper - lnORhat) / 1.96) / 2,
  SE.lnOR = case_when(!is.na(SE.extrap) ~ SE.extrap,
                      !is.na(SE.counts) ~ SE.counts,
                      !is.na(SE.chi) ~ SE.chi),
  Variable = as.character(Variable),
  Variable = case_when(Variable == "Congestive" ~ "CHF",
                       Variable == "Hemoglobin" ~ "Hematocrit",
                       Variable == "Nonwhite Race" ~ "White Race",
                       TRUE ~ Variable),
  Varnum = as.numeric(factor(Variable, levels = unique(Variable))),
  Paper = as.factor(as.character(Paper)),
  Paper_num = as.numeric(Paper),
  counts = is.na(event_exposure) + 1,
  # reverse things for nonwhite race so that log(OR) is positive
  # i.e. it's white race as the risk factor, rather than nonwhite race
  n_i0 = case_when(Variable == "White Race" ~ event_exposure + noevent_exposure,
                   Variable != "White Race" ~ event_noexposure + noevent_noexposure),
  n_i1 = case_when(Variable == "White Race" ~ event_noexposure + noevent_noexposure,
                   Variable != "White Race" ~ event_exposure + noevent_exposure),
  y_i0 = case_when(Variable == "White Race" ~ event_exposure,
                   Variable != "White Race" ~ event_noexposure),
  y_i1 = case_when(Variable == "White Race" ~ event_noexposure,
                   Variable != "White Race" ~ event_exposure),
  N_i = n_i0 + n_i1)

keeps <- c("Variable", "lnORhat", "SE.lnOR", "n_i0", "n_i1", "y_i0", "y_i1", "N_i", "Paper_num", "Paper", "counts", "Varnum")                 
syncope <- univOR %>% 
  group_by(Variable) %>%
  arrange(Varnum, desc(lnORhat)) %>%
  select(all_of(keeps))

write.csv(syncope, here("R", "syncope.cleaned.csv"), row.names = FALSE)

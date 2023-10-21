library(tidyverse)

# Memory and time usage
dt <- "202310160"
ff <- list.files(paste0("Z:/paper1/outputs/", dt, "/lyra_out"), full.names = T)
ll <- list()
for(i in 1:length(ff)){
  ll[[i]] <-suppressMessages(read_csv(ff[i]) %>% 
                               setNames(c("this")) %>% 
                               filter(str_detect(this, "Mem usage|Wall time"))%>% 
                               separate(this, into = c("Metric", "Value"), sep = " : ") %>% 
                               pivot_wider(names_from = Metric, values_from = Value) %>% 
                               setNames(c("wall_time", "mem_usage")))
}

bind_rows(ll) %>% view()
# max 70 GB for 10 reps
# max 70 GB for 20 reps
# max 19.5 hours for 20 reps
# max 83GB for 30 reps
# max 30 hours for 30 reps

# average round run time 
dt <- "202310150"
ff <- list.files(paste0("Z:/paper1/outputs/", dt, "/lyra_errors"), full.names = T)
ll <- list()
for(i in 1:length(ff)){
  ll[[i]] <-suppressMessages(read_csv(ff[i]) %>% 
                               setNames(c("this")) %>% 
                               filter(str_detect(this, "took ")))
}

bind_rows(ll) %>% view()

read_csv("Z:/paper1/outputs/202310122/lyra_errors/scen2Iinst2_") %>% 
  setNames(c("this")) %>% 
  filter(str_detect(this, "took "))

# takes on average 45 mins - should allow for an hour per rep

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(tictoc)

tic()

proj <- "RNAi"

options(digits.secs = 6)

dates_folds <- list.files(".", pattern="20")

dates_folds <- dates_folds[which(!str_detect(dates_folds, ".Rds"))]

ps <- dates_folds

for (p in ps) {
  message(p)
  
  # Get list of files that have tracking results 
  flist <- list.files(
    path = p,
    pattern = "DLC_resnet50_fly_trackerNov28shuffle1_1030000.csv",
    recursive = TRUE, full.names = TRUE)
  
  flist <- flist[str_detect(flist, proj)] 
  
  read_tracks <- function(fname) {
    # message(fname)
    
    # Base id
    id <- str_remove(fname,
                     "DLC_resnet50_fly_trackerNov28shuffle1_1030000.csv")
    id <- str_split(id, "/", simplify = TRUE)[3]
    
    # Read digitized points
    M <- read_csv(fname,
                  col_types = cols(.default = col_character()),
                  skip = 2)[] %>% 
      mutate(across(.cols = everything(), .fns = as.numeric)) %>% 
      rename(frame = coords)
    
    # Convert to distance
    dists <- data.frame(id = rep(id, times = nrow(M) - 1),
                        frame = numeric(nrow(M) - 1),
                        x = round(M$x[2:nrow(M)], 4),
                        y = round(M$y[2:nrow(M)], 4),
                        d = numeric(nrow(M) - 1),
                        likelihood = numeric(nrow(M) - 1))
    for (ii in 1:nrow(dists)) {
      dists$frame[ii] <- M$frame[ii + 1]
      dists$d[ii] <- dist(M[ii:(ii+1), c(2, 3)]) %>% as.numeric()
      dists$likelihood[ii] <- mean(c(M$likelihood[ii], M$likelihood[ii + 1]))
    }
    
    dists <- dists %>% 
      mutate(d = round(d, 3),
             likelihood = round(likelihood, 3))
    
    # Find the right data file for this trial
    tmp_split <- str_split(id, "_", simplify = TRUE)
    
    run_id <- tmp_split[length(tmp_split)] %>% as.numeric()
    if (run_id < 10) {
      id_data <- str_sub(id, end = -3)
    } else {
      id_data <- str_sub(id, end = -4)
    }
    
    # Load data file and convert to POSIXct time
    data_path <- paste0(p, "/", id_data, ".csv")
    D <- read_csv(data_path,
                  col_types = cols(.default = col_double())) %>%
      mutate(
        Microsecond = sprintf("%06d", Microsecond),
        timestamp = as.POSIXct(paste0(Year, "-", Month, "-", Day, " ",
                                      Hour, ":", Minute,
                                      ":", Second, ".", Microsecond))) %>%
      relocate(timestamp)
    
    # Convert to seconds after the start of the recording
    Therm <- tibble(Second = (D$timestamp - D$timestamp[1])[-1] %>% 
                      as.numeric() %>% round(., 3),
                    Thermistor_Temp = round(D$Thermistor_Temp[-1], 3),
                    Analog = D$Analog[-1])
    
    if(nrow(dists)==nrow(Therm)){}else{
    Therm <- Therm[-1,]
    }
    
    return(bind_cols(dists, Therm))
  }
  
  pp <- map(.x = flist,
            .f = read_tracks)
  
  # pp[[95]] %>% slice(1:10)
  # 
  # pp[[95]] %>% 
  #   ggplot(aes(Second, d, color = likelihood)) +
  #   geom_line()
  # 
  # pp[[95]] %>% 
  #   filter(likelihood > 0.97) %>% 
  #   select(Second, Thermistor_Temp, d) %>% 
  #   pivot_longer(cols = -Second) %>% 
  #   ggplot(aes(Second, value, color = name)) +
  #   geom_line() +
  #   facet_grid(name ~ ., scales = "free_y")
  
  write_rds(pp, file = paste0(p, "_Fly_tracks.Rds"))
  length(pp)
  
}

##### Combine

tracks <- list.files(path = ".", pattern = "_Fly_tracks")

combine_tracks <- function(fname) {
  tr <- read_rds(fname)
  tr_comb <- bind_rows(tr)
  return(tr_comb)
}

all_tracks <- map_dfr(.x = tracks,
                      combine_tracks)

write_csv(all_tracks, paste0("Combined_tracks_",proj,".csv"))
write_rds(all_tracks, paste0("Combined_tracks_",proj,".Rds"))

toc()

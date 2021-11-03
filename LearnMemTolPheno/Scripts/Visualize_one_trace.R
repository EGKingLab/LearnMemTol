

TTdat <- readRDS(file="../ProcessedData/Combined_tracks_RNAi_2.Rds")
str(TTdat)
TTdat[1:5,]
ff <- "2021-07-01_RNAi-50656.25750-3_group1"

TTdat <- read.csv(file="../ProcessedData/RNAi1_tracked_data.csv")


tt1 <- TTdat[TTdat$id==ff,]
tt1 <- tt1[tt1$likelihood >= 0.6,]



plot(tt1$Second, tt1$x)

ggplot(tt1, aes(Second, x, color=likelihood)) +
  geom_point()
  



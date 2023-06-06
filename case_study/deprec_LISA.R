## Using LISA

z <- mb_est$draws_mu$TSLN - HT_mu
z <- z[,-c(982,1374)]

mapping_data <- b_est$summ_mu %>%
  filter(model == "TSLN") %>% 
  left_join(.,map_sa2, by = "ps_area") %>%
  st_as_sf() %>%
  st_transform(4326) %>% 
  filter(!st_is_empty(mapping_data))

nb2 <- poly2nb(mapping_data, queen = T) #object new.shp2 used the sf package
W2 <- nb2mat(nb2, style="W", zero.policy = TRUE)

spatial_lag <- z %*% t(W2)

foo <- function(x)mean(x > 0)
spatial_lag_est <- apply(spatial_lag, 2, foo)
z_est <- apply(z, 2, foo)

test <- data.frame(z_est, spatial_lag_est) %>% 
  mutate(highhigh = ifelse(z_est > 0.9 & spatial_lag_est > 0.9, 1, 0),
         lowlow = ifelse(z_est <0.1 & spatial_lag_est <0.1, 1, 0)) %>% 
  bind_cols(mapping_data) %>% 
  st_as_sf() %>%
  st_transform(4326)

library(tmap)
tmap_mode("view")
tm_shape(test)+
  tm_polygons(col = "median")

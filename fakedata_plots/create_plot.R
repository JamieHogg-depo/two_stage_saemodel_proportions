## -----------------------------------------------------------------------------
## create_plot ## --------------------------------------------------------------
## -----------------------------------------------------------------------------

# Functions ### 
base_folder <- "C:/r_proj/two_stage_saemodel_proportions/fakedata_plots"
export <- F #FALSE

jsave <- function(filename, square = T, square_size = 5000, ratio = c(6,9)){
  if(square){
    ggsave(filename = filename,
           path = paste0(base_folder, "/plots"),
           dpi = 1000,
           width = square_size,
           height = square_size,
           scale = 1,
           units = "px")
  }else{
    total = square_size^2
    a <- sqrt((total*ratio[1])/ratio[2])
    b <- (ratio[2]*a)/ratio[1]
    ggsave(filename = filename,
           path = paste0(base_folder, "/plots"),
           dpi = 1000,
           width = round(b),
           height = round(a),
           scale = 1,
           units = "px")
  }
}

# Colors 
hue_pal()(5) # gets the first few default colors from ggplot
jcol <- data.frame(model = c("TSLN", "ELN", "LOG", "HT"),
           color = c("#C77CFF", "#7CAE00", "#00BFC4", "#ae3200"))

# Create subset data
prev_median_wide <- dplyr::select(summ_mu, median, model, ps_area) %>% 
						pivot_wider(names_from = model, values_from = median)
prev_ci_wide <- dplyr::select(summ_mu, cisize, model, ps_area) %>% 
						pivot_wider(names_from = model, values_from = cisize)

# SA4 level - HT vs modeled ####
summ_mu_sa4 %>% 
  ggplot(aes(y = median, ymin = lower, ymax = upper,
             x = HT, xmin = HT_lower, xmax = HT_upper))+
  theme_bw()+geom_abline(col="red")+
  geom_errorbar(col = "grey")+geom_errorbarh(col = "grey")+
  geom_point()+
  facet_grid(.~model)+
  labs(y = "Modelled prevalence estimate",
       x = "Direct prevalence estimate")
if(export) jsave("sa4_directvsmodeled.png", square = F)

# SA4 level - ARB, RRMSE #### 
(summ_mu_sa4 %>% 
  dplyr::select(ps_sa4, n, RRMSE, model) %>% 
  pivot_wider(names_from = model, values_from = RRMSE) %>% 
  arrange(TSLN) %>% mutate(x = 1:nrow(.)) %>% dplyr::select(-ps_sa4) %>% 
  pivot_longer(-c(x, n)) %>% 
  ggplot(aes(y = 10*value, x = x, color = name))+
  theme_bw()+
  geom_path()+geom_point(aes(size = n))+
  scale_color_manual(values = jcol$color,
                    breaks = jcol$model)+
  labs(y = "RRMSE (x10)", x = "", col = "",
       size = "Sample size")+
  theme(legend.position = "bottom",
        legend.box = "vertical")+
   scale_size_binned(n.breaks = 5))+
(summ_mu_sa4 %>% 
  dplyr::select(ps_sa4, n, ARB, model) %>% 
  pivot_wider(names_from = model, values_from = ARB) %>% 
  arrange(TSLN) %>% mutate(x = 1:nrow(.)) %>% dplyr::select(-ps_sa4) %>% 
  pivot_longer(-c(x, n)) %>% 
  ggplot(aes(y = 10*value, x = x, color = name))+
  theme_bw()+
  geom_path()+geom_point(aes(size = n))+
  scale_color_manual(values = jcol$color,
                     breaks = jcol$model)+
  labs(y = "ARB (x10)", x = "", col = "",
       size = "Sample size")+
  theme(legend.position = "bottom",
        legend.box = "vertical")+
   scale_size_binned(n.breaks = 5))
if(export) jsave("sa4_rrmse_arb.png", square = F)

# What about relative measures??

# Violin plots ####
summ_mu %>% 
  ggplot(aes(y = model, fill = model, x = median))+
  theme_bw()+
  geom_violin()+
  labs(y = "",
       x = "Modeled prevalence estimate")+
  theme(legend.position = "none")+
  scale_fill_manual(values = jcol$color,
                    breaks = jcol$model)
if(export) jsave("violin.png", square = F)

# Boxplot of posterior medians - IRSD ####
prev_median_wide %>%
  left_join(.,dplyr::select(aux2, ps_area, ABS_irsd_decile_nation, ra_sa2_3c, N_persons),
				by = "ps_area") %>%
	pivot_longer(-c(ps_area, ABS_irsd_decile_nation, ra_sa2_3c, N_persons)) %>%
  # rename IRSD categories
  mutate(ABS_irsd_decile_nation = factor(ABS_irsd_decile_nation, levels = 1:10,
                                         labels = c("Least advantaged", 
                                                    as.character(2:9), 
                                                    "Most advantaged"))) %>% 
	ggplot(aes(x = value, fill = name, y = ABS_irsd_decile_nation))+
	theme_bw()+
	geom_boxplot()+
  labs(x = "Modeled prevalence estimate",
       y = "IRSD deciles",
       fill = "")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = jcol$color,
                    breaks = jcol$model)
if(export) jsave("boxplot_byseifa.png", square = F)

# Boxplot of posterior medians - Remoteness ####
prev_median_wide %>%
  left_join(.,dplyr::select(aux2, ps_area, ABS_irsd_decile_nation, ra_sa2_3c, N_persons),
            by = "ps_area") %>%
  pivot_longer(-c(ps_area, ABS_irsd_decile_nation, ra_sa2_3c, N_persons)) %>%
  ggplot(aes(x = value, fill = name, y = ra_sa2_3c))+
  theme_bw()+
  geom_boxplot()+
  #facet_grid(.~ra_sa2_3c, labeller = label_wrap_gen(multi_line = TRUE))+
  labs(x = "Modeled prevalence estimate",
       y = "",
       fill = "")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = jcol$color,
                    breaks = jcol$model)
if(export) jsave("boxplot_byremoteness.png", square = F)

# Compare estimates from models ####
dplyr::select(summ_mu, median, lower, upper, model, ps_area) %>% 
  pivot_wider(names_from = model, values_from = c(median, lower, upper)) %>% 
  ggplot(aes(y = median_TSLN, ymin = lower_TSLN, ymax = upper_TSLN,
             x = median_ELN, xmin = lower_ELN, xmax = upper_ELN))+
  theme_bw()+
  geom_abline(col = "red")+
  geom_errorbar(col="grey")+
  geom_errorbarh(col="grey")+
  geom_point()+xlim(0,0.5)+ylim(0,0.5)+
  labs(y = "TSLN model",
       x = "ELN model")
if(export) jsave("scatter_TSLNvsELN.png", square = T)

dplyr::select(summ_mu, median, lower, upper, model, ps_area) %>% 
  pivot_wider(names_from = model, values_from = c(median, lower, upper)) %>% 
  ggplot(aes(y = median_TSLN, ymin = lower_TSLN, ymax = upper_TSLN,
             x = median_LOG, xmin = lower_LOG, xmax = upper_LOG))+
  theme_bw()+
  geom_abline(col = "red")+
  geom_errorbar(col="grey")+
  geom_errorbarh(col="grey")+
  geom_point()+xlim(0,0.5)+ylim(0,0.5)+
  labs(y = "TSLN model",
       x = "LOG model")
if(export) jsave("scatter_TSLNvsLOG.png", square = T)

dplyr::select(summ_mu, median, lower, upper, model, ps_area) %>% 
  pivot_wider(names_from = model, values_from = c(median, lower, upper)) %>% 
  ggplot(aes(y = median_ELN, ymin = lower_ELN, ymax = upper_ELN,
             x = median_LOG, xmin = lower_LOG, xmax = upper_LOG))+
  theme_bw()+
  geom_abline(col = "red")+
  geom_errorbar(col="grey")+
  geom_errorbarh(col="grey")+
  geom_point()+xlim(0,0.5)+ylim(0,0.5)+
  labs(y = "ELN model",
       x = "LOG model")
if(export) jsave("scatter_ELNvsLOG.png", square = T)

# Caterpillar plots - by model - with CI  #### ---------------------------------

# Prevalence
summ_mu %>% 
  group_by(model) %>% 
  mutate(rank = rank(median, ties.method = "first")) %>%
  ggplot(aes(y = median, ymin = lower, ymax = upper,
             x = rank, 
             col = median))+theme_bw()+
  geom_errorbar(col = "grey")+
  geom_point()+
  facet_wrap(.~model)+
  scale_color_viridis_c(begin = 0.3, end = 1, 
                        direction = 1,
                        option = "B")+
  geom_hline(yintercept = HT_mu)+
  labs(y = "Modeled prevalence estimates",
       x = "")+
  theme(legend.position = "none")
if(export) jsave("cat_prev_medianci.png", square = F)

# Caterpillar plots - by population #### ---------------------------------------

# Prevalence
prev_median_wide %>%
	left_join(.,dplyr::select(aux2, N_persons, ps_area),
	by = "ps_area") %>%
	left_join(.,dplyr::select(sample_agg, ps_area, HT),
	by = "ps_area") %>%
	arrange(TSLN) %>%
	mutate(x = 1:nrow(.)) %>%
	pivot_longer(-c(ps_area, N_persons, x)) %>%
	ggplot(aes(y = value, x = x, col = name))+theme_bw()+
	geom_point()+
	facet_wrap(.~cut_number(N_persons, 4))+
  scale_color_manual(values = jcol$color,
                    breaks = jcol$model)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
if(export) jsave("catBN_prev_median.png", square = F)
	
# Prevalence (CI)
prev_ci_wide %>%
	left_join(.,dplyr::select(aux2, N_persons, ps_area),
	by = "ps_area") %>%
  left_join(.,dplyr::select(direct_est, ps_area, cisize),
            by = "ps_area") %>%
	arrange(TSLN) %>%
	mutate(x = 1:nrow(.)) %>%
	pivot_longer(-c(ps_area, N_persons, x)) %>%
	ggplot(aes(y = value, x = x, col = name))+
	geom_point()+
	facet_wrap(.~cut_number(N_persons, 4))+
  scale_color_manual(values = jcol$color,
                     breaks = jcol$model)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
if(export) jsave("catBN_prev_cisize.png", square = F)

# ORs
dplyr::select(summ_or, median, model, ps_area) %>%
	pivot_wider(names_from = model, values_from = median) %>%
	left_join(.,dplyr::select(aux2, N_persons, ps_area),
	by = "ps_area") %>%
	arrange(TSLN) %>%
	mutate(x = 1:nrow(.)) %>%
	pivot_longer(-c(ps_area, N_persons, x)) %>%
	ggplot(aes(y = value, x = x, col = name))+
	geom_point()+
	facet_wrap(.~cut_number(N_persons, 4))+
  scale_color_manual(values = jcol$color,
                     breaks = jcol$model)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
if(export) jsave("catBN_or_median.png", square = F)
	
# ORs (CI)
dplyr::select(summ_or, cisize, model, ps_area) %>%
	pivot_wider(names_from = model, values_from = cisize) %>%
	left_join(.,dplyr::select(aux2, N_persons, ps_area),
	by = "ps_area") %>%
	arrange(TSLN) %>%
	mutate(x = 1:nrow(.)) %>%
	pivot_longer(-c(ps_area, N_persons, x)) %>%
	ggplot(aes(y = value, x = x, col = name))+
	geom_point()+
	facet_wrap(.~cut_number(N_persons, 4))+
  scale_color_manual(values = jcol$color,
                     breaks = jcol$model)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
if(export) jsave("catBN_or_cisize.png", square = F)

# City Insets #### -------------------------------------------------------------

lims <- data.frame(
  xmin = c(152.6, 150.35, 144.5, 115.45, 138.1, 146.8, 148.6, 130.3),
  xmax = c(153.6, 151.35, 145.5, 116.45, 139.1, 147.8, 149.6, 131.3),
  ymin = -c(28, 34.4, 38.4, 32.5, 35.4, 43.4, 35.8, 13),
  ymax = -c(27, 33.4, 37.4, 31.5, 34.4, 42.4, 34.8, 12),
  city = c("Brisbane", "Sydney", "Melbourne", "Perth", "Adelaide", "Hobart", "Canberra", "Darwin")
)

# Maps: Prevalence #### --------------------------------------------------------
    # posterior medians (left)
    # size of CI (right)
    # grid for three models (along y axis)

direct_est <- sample_agg %>% 
  dplyr::select(ps_area, HT, cisize) %>% 
  rename(median = HT) %>% 
  mutate(model = "Direct")
mapping_data <- summ_mu %>%
  bind_rows(direct_est) %>% 
  left_join(.,map_sa2, by = "ps_area") %>%
  st_as_sf() %>%
  st_transform(4326)

# Create base map for prevalence
(bm_mu <- mapping_data %>% 
   ggplot(aes(fill = median))+
   theme_void()+
   geom_sf(col = NA)+
   facet_grid(.~model)+
   scale_fill_viridis_c(begin = 0.3, end = 1, 
                        direction = 1,
                        option = "B")+
   labs(fill = "Prevalence")+
   theme(legend.position = "right", legend.key.height = unit(0.5, "cm")))

# Create base map for CI prevalence
(bm_muci <- mapping_data %>% 
    ggplot(aes(fill = cisize))+
    theme_void()+
    geom_sf(col = NA)+
    facet_grid(.~model)+
    scale_fill_viridis_c(begin = 0, end = 0.8, 
                         direction = -1,
                         option = "D")+
    labs(fill = "HDI size")+
    theme(legend.position = "right", legend.key.height = unit(0.5, "cm"),
          strip.background = element_blank(),
          strip.text.x = element_blank()))

# Export full maps
bm_mu/bm_muci
if(export) jsave("map_mu.png")

# Subset to capital cities
cities <- lims[c(1,2,3,7),]
for(i in 1:nrow(cities)){
  mu <- bm_mu +
    xlim(cities$xmin[i], cities$xmax[i]) +
    ylim(cities$ymin[i], cities$ymax[i]) +
    ggtitle(label = cities$city[i])
  muci <- bm_muci +
    xlim(cities$xmin[i], cities$xmax[i]) +
    ylim(cities$ymin[i], cities$ymax[i])
  mu/muci
  jsave(paste0("map_insets/map_mu_", cities$city[i], ".png"), square = F)
}

# Maps: ORs #### --------------------------------------------------------------

# SETUP
cut_offs <- c(1/1.5, 1.5)
direct_est <- sample_agg %>% 
  dplyr::select(ps_area, OR, OR_lower, OR_upper) %>% 
  rename(median = OR) %>% 
  mutate(model = "Direct",
         cisize = OR_upper - OR_lower) %>% 
  dplyr::select(-c(OR_lower, OR_upper))
mapping_data <- summ_or %>%
  bind_rows(direct_est) %>% 
	left_join(.,map_sa2, by = "ps_area") %>%
	st_as_sf() %>%
	st_transform(4326) %>%
	mutate(median = ifelse(median > cut_offs[2], cut_offs[2], median),
	       median = ifelse(median < cut_offs[1], cut_offs[1], median))
		   
# define fill colours
Fill.colours <- c("#2C7BB6", "#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C", "#D7191C")
End <- log(1.6)
Breaks.fill <- c(1/1.5, 1/1.25, 1, 1.25, 1.5)
Fill.values <- c(-End, log(Breaks.fill), End)

# Create base map for ORs
bm_or <- mapping_data %>%
	ggplot(aes(fill = log(median)))+
	theme_void()+
	geom_sf(col = NA)+
	facet_grid(.~model)+
	scale_fill_gradientn(colors = Fill.colours,
							values = rescale(Fill.values),
							labels = as.character(round(Breaks.fill, 3)),
							breaks = log(Breaks.fill),
							limits = range(Fill.values))+
	labs(fill = "OR")+
	theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),)
	
# Create base map for cisize of ORs
bm_orci <- mapping_data %>%
	ggplot(aes(fill = cisize))+
	theme_void()+
	geom_sf(col = NA)+
	facet_grid(.~model)+
  scale_fill_viridis_c(begin = 0, end = 0.8, 
                       direction = -1,
                       option = "D")+
  labs(fill = "HDI size")+
  theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Maps: EPs for ORs #### ----------------------------------------------

# SETUP
mapping_data <- summ_or %>%
  bind_rows(data.frame(model = "Direct", ps_area = 1:1695)) %>% 
	left_join(.,map_sa2, by = "ps_area") %>%
	st_as_sf() %>%
	st_transform(4326)
	
# Create base map for EPs
(bm_ep <- mapping_data %>%
	ggplot(aes(fill = EP))+
	theme_void()+
	geom_sf(col = NA)+
	facet_grid(.~model)+
  scale_fill_viridis_c(begin = 0, end = 0.8,
                       direction = -1,
                       option = "C", 
                       limits = c(0.01,0.99),
                       #oob = squish,
                       breaks = c(0,0.2,0.5,0.8,1),
                       labels = as.character(c(0,0.2,0.5,0.8,1)),
                       trans = "logit")+
	labs(fill = "EP")+
	theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
	      strip.background = element_blank(),
	      strip.text.x = element_blank()))
	
# Export full maps
bm_or/bm_ep/bm_orci
if(export) jsave("map_or.png")

# Subset to capital cities
cities <- lims[c(1,2,3,7),]
for(i in 1:nrow(cities)){
  or <- bm_or +
    xlim(cities$xmin[i], cities$xmax[i]) +
    ylim(cities$ymin[i], cities$ymax[i]) +
    ggtitle(label = cities$city[i])
  orci <- bm_orci +
    xlim(cities$xmin[i], cities$xmax[i]) +
    ylim(cities$ymin[i], cities$ymax[i])
  orep <- bm_ep +
    xlim(cities$xmin[i], cities$xmax[i]) +
    ylim(cities$ymin[i], cities$ymax[i])
  (or/orep/orci)
  jsave(paste0("map_insets/map_or_", cities$city[i], ".png"), square = F)
}

## ----------------------------------------------------------------------------
## Alterntaive #### -----------------------------------------------------------
## ----------------------------------------------------------------------------

## Quantils #### ---------------------------------------------------------------

# SETUP
# cut_offs <- unname(quantile(summ_mu$median, probs = c(0.1,0.9)))
# mapping_data <- summ_mu %>% 
#   bind_rows(direct_est) %>% 
#   left_join(.,map_sa2, by = "ps_area") %>% 
#   st_as_sf() %>% 
#   st_transform(4326) %>% 
#   mutate(median = ifelse(median > cut_offs[2], cut_offs[2], median),
#          median = ifelse(median < cut_offs[1], cut_offs[1], median))
# 
# # define fill colours
# fill.colors <- c('#FEEDDE','#FDD0A2','#FDAE6B','#FD8D3C','#F16913','#D94801','#8C2D04')
# End <- max(summ_mu$median)
# Start <- min(summ_mu$median)
# Fill.values <- unname(quantile(summ_mu$median, probs = seq(0.1, 0.9, length.out  = 7)))
# Breaks.fill <- Fill.values[c(-1,-7)]
# 
# # Create base map for prevalence
# bm_mu <- mapping_data %>% 
#   ggplot(aes(fill = median))+
#   theme_void()+
#   geom_sf(col = NA)+
#   facet_grid(.~model)+
#   scale_fill_gradientn(colors = fill.colors,
#                        values = rescale(Fill.values),
#                        labels = as.character(round(Breaks.fill, 2)),
#                        breaks = Breaks.fill,
#                        limits = range(Fill.values))+
#   labs(fill = "Prevalence")+
#   theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))
# 
# # Create base map for cisize of prevalence
# bm_muci <- mapping_data %>%
#   ggplot(aes(fill = cisize))+
#   theme_void()+
#   geom_sf(col = NA)+
#   facet_grid(.~model)+
#   scale_fill_gradient2(low = "light blue",
#                        high = "blue",
#                        limits = c(0, 0.3), oob = squish)+
#   labs(fill = "HDI size")+
#   theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
#         strip.background = element_blank(),
#         strip.text.x = element_blank())
# 
# # Export full maps
# bm_mu/bm_muci
# if(export) jsave("map_mu.png")
# 
# # Subset to capital cities
# cities <- lims[c(1,2,3,7),]
# for(i in 1:nrow(cities)){
#   mu <- bm_mu +
#     #guides(fill = "none")+ # removes legend
#     xlim(cities$xmin[i], cities$xmax[i]) +
#     ylim(cities$ymin[i], cities$ymax[i]) +
#     ggtitle(label = cities$city[i])
#   muci <- bm_muci +
#     #guides(fill = "none")+ # removes legend
#     xlim(cities$xmin[i], cities$xmax[i]) +
#     ylim(cities$ymin[i], cities$ymax[i])
#   mu/muci
#   jsave(paste0("map_prevalence_", cities$city[i], ".png"), square = F)
# }
# 
# # Add black box to map
# bm_mu +
#   annotate("rect",
#            xmin = lims$xmin[2] - 0.1, xmax = lims$xmax[2] + 0.1,
#            ymin = lims$ymin[2] - 0.1, ymax = lims$ymax[2] + 0.1,
#            colour = "black", fill = NA)
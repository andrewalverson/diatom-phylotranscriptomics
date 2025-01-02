library(here)
library(readxl)
library(tidyverse)
library(viridis)
library(cowplot)
library(wesanderson)

voucher_metadata <- read_excel(here('cultures-metadata', 'voucher-list.xlsx'))
voucher_metadata <- voucher_metadata %>% filter(!is.na(label))

# load eukaryotic and stramenopile busco summaries
busco.stram <- read.csv(here('busco', 'busco_stramenopiles_summary.csv'))
busco.euk   <- read.csv(here('busco', 'busco_eukaryote_summary.csv'))
omark <- read.csv(here('omark', 'omark_results.csv'))

# join the voucher and busco tables
dat <- left_join(voucher_metadata, busco.stram, by = c("label" = "tip")) 
dat <- left_join(dat, busco.euk, by = c("label" = "tip"))
dat <- left_join(dat, omark, by = "label") 

# add the fragmented + complete BUSCOs
dat <- dat %>% mutate(euk.busco.total = euk.complete+euk.fragmented)
dat <- dat %>% mutate(stram.busco.total = stram.complete+stram.fragmented)

# calculate total percent HOGs found
dat <- dat %>% mutate(hogs.total = 100-hogs.missing)

# set the color palette
pal <- wes_palette("Darjeeling1", n = 5)
pal <- pal[-c(1,3)]

# calculate means to report in the paper
dat %>% group_by(data.source) %>% 
  summarise(Mean = mean(euk.busco.total, na.rm = TRUE))

dat %>% group_by(data.source) %>% 
  summarise(Mean = mean(stram.busco.total, na.rm = TRUE))

dat %>% group_by(data.source) %>% 
  summarise(Mean = mean(hogs.total, na.rm = TRUE))

# plot busco for stramenopile db
busco.stram.plot <- dat %>% 
  group_by(data.source) %>%
  ggplot(aes(x = data.source, y = stram.busco.total, color=data.source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 1.25, stroke = 0.005) +
  scale_color_manual(values=pal) +
  theme_bw() +
  labs(color="Data source") +
  scale_x_discrete(labels = c("Genome", "MMETSP", "This\nstudy")) +
  ylab("Percent Stramenopile BUSCOs") +
  ylim(0, 102) +
  theme(axis.text.x = element_text(size = 6, color = "grey30"),
        axis.text.y = element_text(size = 6, color = "grey30"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "grey30", size = 6.5),
        axis.ticks = element_line(color = "grey30", linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        legend.position= "none",
        legend.title = element_text(size = 6, colour = "grey30"),
        legend.text = element_text(size = 6, color = "grey30"),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1)
  )

# ggsave(busco.plot, file = here("figs", "busco-stram.pdf"), height = 4.35, width = 4.35, units = 'cm')

# plot busco for eukaryotes db
busco.euk.plot <- dat %>% 
  group_by(data.source) %>%
  ggplot(aes(x = data.source, y = euk.busco.total, color=data.source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 1.25, stroke = 0.005) +
  scale_color_manual(values=pal) +
  theme_bw() +
  labs(color="Data source") +
  scale_x_discrete(labels = c("Genome", "MMETSP", "This\nstudy")) +
  ylab("Percent Eukaryotic BUSCOs") +
  ylim(0, 102) +
  theme(axis.text.x = element_text(size = 6, color = "grey30"),
        axis.text.y = element_text(size = 6, color = "grey30"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "grey30", size = 6.5),
        axis.ticks = element_line(color = "grey30", linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        legend.position= "none",
        legend.title = element_text(size = 6, colour = "grey30"),
        legend.text = element_text(size = 6, color = "grey30"),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1)
  )

# ggsave(busco.plot, file = here("figs", "busco-stram.pdf"), height = 4.35, width = 4.35, units = 'cm')

# plot omark HOGs
omark.plot <- dat %>% 
  group_by(data.source) %>% 
  ggplot(aes(x = data.source, y = hogs.total, color=data.source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 1.25, stroke = 0.01) +
  scale_color_manual(values=pal) +
  theme_bw() +
  labs(color="Data source") +
  scale_x_discrete(labels = c("Genome", "MMETSP", "This\nstudy")) +
  ylab("Percent Conserved HOGs") +
  ylim(0, 102) +
  theme(axis.text.x = element_text(size = 6, color = "grey30"),
        axis.text.y = element_text(size = 6, color = "grey30"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "grey30", size = 6.5),
        axis.ticks = element_line(color = "grey30", linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        legend.position= "none",
        legend.title = element_text(size = 6, colour = "grey30"),
        legend.text = element_text(size = 6, color = "grey30"),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1)
  )

# ggsave(omark.plot, file = here("figs", "omark-hogs.pub.pdf"), height = 4.35, width = 4.35, units = 'cm')

# 
figure.publication <- plot_grid(busco.euk.plot, omark.plot)
ggsave(figure.publication, file = here("figs", "busco-omark.pub.pdf"), height = 4.35, width = 8.7, units = 'cm')

# save png files for powerpoint
# ggsave(euk_plot, filename = here("figs", "busco-scores-euk.png"),
#        height = 4, width = 4)
# ggsave(stram_plot, filename = here("figs", "busco-scores-stram.png"),
#        height = 4, width = 4)
# ggsave(figure, filename = here("figs", "busco-scores.png"),
#        height = 4, width = 8)
# ggsave(figure, filename = here("figs", "busco-scores.pdf"),
#        height = 4, width = 8)


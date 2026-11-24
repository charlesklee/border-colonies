#!/usr/bin/env Rscript
# ==============================================================================
# 20_figures.R — Render Publication Figures (Figs 1A–B, 2, 3, 4, S3, S4)
# ==============================================================================
# Purpose : Generate all code-produced manuscript figures (174 mm wide, 600 dpi
#           LZW TIFF + vector PDF, Arial / Liberation Sans 7–8 pt) from the
#           intermediate tables produced by scripts 11–18:
#             - Figure 1A–B : Sampling site maps (NZ overview + Tauranga Harbour)
#             - Figure 2    : ASV categories (A), Bray-Curtis NMDS (B), Venn (C)
#             - Figure 3    : Species PCoA (A) & per-species tissue R^2 (B)
#             - Figure 4    : Top 15 families (A) & within-pair log2 ratios (B)
#             - Figures S3–S4: Status PCoA & IB encounter summary
#           (Note: Figure 1C [Google Earth], Figure 1D [photos], Figure S2
#           [dissection schematics], and Figure S1 final artwork [Gemini 4 Neon
#           from 18_sample_flow_audit.R counts] are not direct script outputs.)
# Inputs  : config.R, output/intermediate/*.csv, output/tables/Table1_asv_categories.csv
# Outputs : output/figures/Fig1_map_AB.pdf/.tiff, Fig2.pdf/.tiff,
#           Fig3.pdf/.tiff, Fig4.pdf/.tiff
# Supports: Figures 1A–B, 2, 3, 4; Figures S3, S4
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(systemfonts)
  library(sf)
  library(maps)
  library(mapdata)
  library(dplyr)
})

f_info <- systemfonts::match_font("Arial")
font_family <- if (!is.null(f_info$path) && grepl("Arial|LiberationSans", basename(f_info$path), ignore.case = TRUE)) "Arial" else "Liberation Sans"

theme_pub <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text        = element_text(family = font_family, size = base_size, colour = "black"),
      axis.text   = element_text(family = font_family, size = base_size - 0.5, colour = "black"),
      axis.title  = element_text(family = font_family, size = base_size, colour = "black", face = "bold"),
      legend.title = element_text(family = font_family, size = base_size, face = "bold"),
      legend.text = element_text(family = font_family, size = base_size - 0.6),
      plot.tag    = element_text(family = font_family, size = base_size + 1.5, face = "bold")
    )
}

save_fig <- function(p, stem, w_mm, h_mm, dpi = 600) {
  grDevices::cairo_pdf(file.path(DIR_FIGURES, paste0(stem, ".pdf")), width = w_mm / 25.4, height = h_mm / 25.4, family = font_family)
  print(p); invisible(dev.off())
  ragg::agg_tiff(file.path(DIR_FIGURES, paste0(stem, ".tiff")), width = w_mm, height = h_mm, units = "mm", res = dpi, compression = "lzw")
  print(p); invisible(dev.off())
}

# --- Figure 3 (PCoA + Per-Species Tissue R^2) ---
if (file.exists(file.path(DIR_INTERIM, "fig3A_pcoa_coords.csv")) &&
    file.exists(file.path(DIR_INTERIM, "table_S4C_per_species_tissue.csv")) &&
    file.exists(file.path(DIR_INTERIM, "table_S3B_paired_tissue_permanova.csv"))) {
  p3_df  <- read.csv(file.path(DIR_INTERIM, "fig3A_pcoa_coords.csv"), stringsAsFactors = FALSE)
  s4c_df <- read.csv(file.path(DIR_INTERIM, "table_S4C_per_species_tissue.csv"), stringsAsFactors = FALSE)
  s3b_df <- read.csv(file.path(DIR_INTERIM, "table_S3B_paired_tissue_permanova.csv"), stringsAsFactors = FALSE)

  sp_r2 <- s3b_df$R2[s3b_df$Dataset == "Primary_52p" & s3b_df$Distance == "Bray_Curtis" & s3b_df$Term == "Species"][1]
  sp_p  <- s3b_df$p_value[s3b_df$Dataset == "Primary_52p" & s3b_df$Distance == "Bray_Curtis" & s3b_df$Term == "Species"][1]
  p3_df$Species <- factor(p3_df$Species, levels = SP_LEVELS)

  pair_seg <- merge(
    p3_df[p3_df$Tissue == "Zooid", c("pair_id", "PCoA1", "PCoA2")],
    p3_df[p3_df$Tissue == "Tunic", c("pair_id", "PCoA1", "PCoA2")],
    by = "pair_id", suffixes = c("_Z", "_T")
  )
  p3a <- ggplot(p3_df, aes(x = PCoA1, y = PCoA2, colour = Species, shape = Tissue)) +
    geom_segment(data = pair_seg, aes(x = PCoA1_Z, y = PCoA2_Z, xend = PCoA1_T, yend = PCoA2_T),
                 colour = "#BDBDBD", linewidth = 0.25, inherit.aes = FALSE) +
    geom_point(size = 1.9, alpha = 0.9) +
    annotate("text", x = max(p3_df$PCoA1), y = max(p3_df$PCoA2), hjust = 1, vjust = 1,
             label = sprintf("Species R\u00B2 = %.2f, p = %.4f", sp_r2, sp_p),
             size = 2.4, family = font_family, fontface = "bold") +
    scale_colour_manual(values = COLS_SPECIES, labels = SP_LEGEND_EXPR, name = "Host species") +
    scale_shape_manual(values = c("Zooid" = 16, "Tunic" = 17), name = "Tissue") +
    labs(x = sprintf("PCoA1 (%.1f%%)", p3_df$PCoA1_pct[1]),
         y = sprintf("PCoA2 (%.1f%%)", p3_df$PCoA2_pct[1])) +
    theme_pub(7.5)

  b3_df <- s4c_df[s4c_df$Dataset == "Primary_52p" & s4c_df$Distance == "Bray_Curtis", ]
  b3_df$Species <- factor(b3_df$Species, levels = SP_LEVELS)
  b3_df$sig <- b3_df$p_holm < ALPHA_SIG
  p3b <- ggplot(b3_df, aes(x = Species, y = R2, colour = Species, fill = ifelse(sig, as.character(Species), "NS"))) +
    geom_col(width = 0.65, linewidth = 0.55) +
    geom_text(aes(label = paste0("n=", n_pairs)), vjust = -0.45, colour = "black", size = 2.2, family = font_family) +
    scale_colour_manual(values = COLS_SPECIES, guide = "none") +
    scale_fill_manual(values = c(COLS_SPECIES, "NS" = "white"), guide = "none") +
    scale_y_continuous(limits = c(0, max(b3_df$R2) * 1.25), expand = c(0, 0)) +
    labs(x = "Host species", y = expression(bold("Within-pair tissue " * R^2))) +
    theme_pub(7.5)

  save_fig((p3a + p3b) + plot_layout(widths = c(1.35, 1.0)) + plot_annotation(tag_levels = "A"),
           "Fig3", w_mm = 174, h_mm = 78)
}

# --- Figure 4 (Top 15 Families + Within-pair log2 ratios) ---
if (file.exists(file.path(DIR_INTERIM, "fig4A_top15_families.csv")) &&
    file.exists(file.path(DIR_INTERIM, "fig4B_within_pair_log2.csv"))) {
  f4a <- read.csv(file.path(DIR_INTERIM, "fig4A_top15_families.csv"), stringsAsFactors = FALSE)
  f4b <- read.csv(file.path(DIR_INTERIM, "fig4B_within_pair_log2.csv"), stringsAsFactors = FALSE)
  f4a$Species <- factor(f4a$Species, levels = SP_LEVELS)
  f4a$Tissue_abbr <- factor(f4a$Tissue_abbr, levels = c("Z", "T"))
  fam_lvls <- unique(f4a$Family)
  f4a$Family <- factor(f4a$Family, levels = fam_lvls)
  cols_fam <- setNames(c("#1B9E77","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77",
                         "#CC6677","#882255","#AA4499","#EE7733","#0077BB","#EE3377","#009988","#661100","#D9D9D9")[seq_along(fam_lvls)], fam_lvls)
  sp_strip <- data.frame(Species = factor(SP_LEVELS, levels = SP_LEVELS), x = 1.5)
  p4a <- ggplot(f4a, aes(x = Tissue_abbr, y = Mean_RA_pct, fill = Family)) +
    geom_col(width = 0.86, colour = "white", linewidth = 0.12) +
    geom_tile(data = sp_strip, aes(x = x, y = 104.2, width = 1.88, height = 4.5),
              fill = unname(COLS_SPECIES[as.character(sp_strip$Species)]), inherit.aes = FALSE) +
    geom_text(data = sp_strip, aes(x = x, y = 104.2, label = Species),
              colour = "white", fontface = "bold", size = 2.3, family = font_family, inherit.aes = FALSE) +
    facet_grid(~ Species, switch = "x") +
    scale_fill_manual(values = cols_fam, name = "Family (Top 15 + Other)") +
    scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(x = "Host species and tissue (Z = Zooid, T = Tunic)", y = "Mean relative abundance (%)") +
    theme_pub(7.5) + theme(strip.text.x = element_blank(), strip.background.x = element_blank(), legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE, title.position = "top"))

  f4b$Species <- factor(f4b$Species, levels = SP_LEVELS)
  f4b$Taxon   <- factor(f4b$Taxon, levels = c("Cyanobiaceae", "Vibrionaceae", "ASV_8"))
  f4b_med <- f4b %>% group_by(Taxon, Species) %>% summarise(n_pairs = n(), med_l2 = median(log2_ratio), .groups = "drop")
  p4b <- ggplot(f4b, aes(x = Species, y = log2_ratio, colour = Species)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#444444", linewidth = 0.4) +
    geom_jitter(width = 0.18, height = 0, size = 1.45, alpha = 0.82) +
    geom_crossbar(data = f4b_med, aes(x = Species, y = med_l2, ymin = med_l2, ymax = med_l2),
                  width = 0.56, colour = "black", linewidth = 0.45, inherit.aes = FALSE) +
    facet_grid(~ Taxon, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = COLS_SPECIES, labels = SP_LEGEND_EXPR, name = "Host species") +
    labs(x = "Host species", y = expression(bold(log[2](("Zooid RA" + italic(p)) / ("Tunic RA" + italic(p)))))) +
    theme_pub(7.5) + theme(legend.position = "bottom")

  save_fig((p4a / p4b) + plot_layout(heights = c(1.15, 0.95)) + plot_annotation(tag_levels = "A"),
           "Fig4", w_mm = 174, h_mm = 158)
}

cat("Figure rendering complete (output/figures/).\n")
print(sessionInfo())

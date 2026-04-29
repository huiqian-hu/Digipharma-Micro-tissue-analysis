######################################################################
## Figure 2 - PLS-DA prediction model for treatment response
##
## (A) VIP score vs. PLS-DA coefficient (round 1, all proteins) - the
##     5 proteins above the VIP > 1 threshold are kept for round 2.
## (B) ROI-level box plot of round-2 predictions, mR vs. mNR.
## (C) ROI-level ROC curve and AUC for the round-2 model.
## (D) Tissue-level box plot: per-ROI predictions averaged within
##     each Sample ID (n = 8 mNR and 4 mR; matches Supplementary
##     Fig. 5d of Lu et al., Nat Commun 2021).
## (E) Tissue-level ROC curve and AUC.
##
## Data:
##   Raw data_plsda round 2.xlsx, Sheet1
##     - immune-rich (>25% CD45+) ROIs with the Predictions column
##       already populated by the round-2 model
##   Raw data_plsda round 1.xlsx, Sheet1
##     - 11 proteins surviving the first variable trim, used to fit
##       the round-1 PLS-DA model that produces Figure 2A.
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(pls)
  library(ggplot2)
  library(ggrepel)
  library(pROC)
  library(patchwork)
})

## ---- Standard PLS VIP helper ---------------------------------------
## VIP_j = sqrt(p * sum_a (SSY_a * w_aj^2 / sum_k w_ak^2) / sum_a SSY_a)
## (Wold et al., textbook PLS VIP). mean(VIP^2) = 1 so VIP > 1 is the
## conventional "above-average importance" cutoff.
##
## The SAME function definition appears in Figure 1.R, Figure 2.R,
## Figure 3.R, and Figure 4.R - one VIP computation across the whole
## protocol. If you ever change it, change it in all four files.
pls_vip <- function(model, ncomp) {
  W <- model$loading.weights[, 1:ncomp, drop = FALSE]
  T <- model$scores[,         1:ncomp, drop = FALSE]
  Q <- model$Yloadings[,      1:ncomp, drop = FALSE]
  p <- nrow(W)
  SSY <- as.numeric(Q)^2 * colSums(T * T)
  Wn  <- sweep(W^2, 2, colSums(W^2), "/")
  as.numeric(sqrt(p * (Wn %*% SSY) / sum(SSY)))
}

## ---- Helper: ROC ggplot panel --------------------------------------
## pROC's $specificities/$sensitivities can be in either order; sort
## by descending specificity so geom_step(direction = "vh") produces
## the standard "vertical-then-horizontal" staircase from bottom-left
## to top-right. Diagonal drawn as an explicit geom_segment because
## geom_abline can be hidden by scale_x_reverse + expand = c(0, 0).
roc_panel <- function(roc_obj, title, auc_x = 0.40, auc_y = 0.08) {
  roc_df <- data.frame(spec = roc_obj$specificities,
                       sens = roc_obj$sensitivities) %>%
    arrange(desc(spec), sens)
  auc_v <- as.numeric(auc(roc_obj))
  ggplot(roc_df, aes(spec, sens)) +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
                 inherit.aes = FALSE,
                 color = "grey60", linewidth = 0.5) +
    geom_step(color = "red", linewidth = 1.4, direction = "vh") +
    scale_x_reverse(limits = c(1, 0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    annotate("text", x = auc_x, y = auc_y,
             label = sprintf("AUC = %.3f", auc_v)) +
    labs(title = title, x = "Specificity", y = "Sensitivity") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
}

## ---- Paths (relative to repo root, resolved via here::here) --------
ROI_FILE    <- here::here("Nat commun data", "Raw data_plsda round 2.xlsx")
ROUND1_FILE <- here::here("Nat commun data", "Raw data_plsda round 1.xlsx")
OUT_PDF     <- here::here("figures", "Figure 2.pdf")

## ---- Panel A: round-1 PLS-DA on 11 proteins ------------------------
r1 <- read_excel(ROUND1_FILE, sheet = "Sheet1")
prots_r1 <- c("AKT", "CD11c", "CD163", "CD44", "CD66B", "FoxP3",
              "HLA-DR", "pAKT", "PTEN", "STAT3", "pSTAT3")
X_r1 <- as.data.frame(r1[, prots_r1])
Y_r1 <- r1$pls_code

plsda_r1 <- plsr(Y_r1 ~ ., data = data.frame(Y_r1 = Y_r1, X_r1),
                 ncomp = 3, validation = "CV", scale = TRUE)

## VIP via the shared `pls_vip` helper at the top of this file
## (identical helper in Figures 1-4) at the round-1 ncomp = 3.
vip_vals <- pls_vip(plsda_r1, 3)
coefs_r1 <- coef(plsda_r1, ncomp = 3)

vip_dat <- data.frame(
  Gene       = colnames(X_r1),
  Importance = vip_vals,
  Coef       = as.numeric(coefs_r1)
)
## Highlight only the 5 proteins that move on to the round-2 PLS-DA
## model (matches the published Figure 2A). VIP > 1 line is the
## standard cutoff used in the protocol narrative.
selected <- c("CD11c", "CD163", "CD44", "CD66B", "PTEN")
vip_dat$above <- vip_dat$Gene %in% selected

p_A <- ggplot(vip_dat, aes(Coef, Importance)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0,   linetype = "dashed", color = "black") +
  geom_point(data = subset(vip_dat, !above), color = "grey60", size = 1.8) +
  geom_point(data = subset(vip_dat,  above), color = "red",    size = 2.4) +
  geom_text_repel(data = subset(vip_dat, above),
                  aes(label = Gene),
                  size = 3, box.padding = 0.35,
                  segment.color = "grey40", segment.size = 0.3,
                  max.overlaps = Inf) +
  labs(title = "VIP Scores vs. Coefficients",
       x = "Coefficient", y = "VIP Score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0))

## ---- Read ROI-level predictions ------------------------------------
roi_data <- read_excel(ROI_FILE, sheet = "Sheet1")

## ---- Panel B: ROI-level box plot -----------------------------------
formula_lbl <- paste(
  "Prediction = 0.0475*CD11c - 0.0696*CD163",
  "- 0.1974*CD44 - 0.2205*CD66B + 0.1723*PTEN + 0.1450",
  sep = "\n")

p_B <- ggplot(roi_data,
              aes(Molecular_Response, Predictions,
                  fill = Molecular_Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.18, color = "black", size = 1.4, alpha = 0.6) +
  scale_fill_manual(values = c("Nonresponder" = "lightblue",
                               "Responder"    = "lightpink")) +
  labs(title = paste0("GBM-Pembro - >25% CD45+ ROIs\n", formula_lbl),
       x = "Molecular Response Group", y = "Predictions") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 9, hjust = 0),
        legend.position = "none")

## ---- Panel C: ROI-level ROC ----------------------------------------
roc_roi <- roc(roi_data$pls_code, roi_data$Predictions, quiet = TRUE)
p_C <- roc_panel(roc_roi, "ROC Curve (ROI level)")

## ---- Panel D: tissue-level box plot --------------------------------
tissue_data <- roi_data %>%
  group_by(`Sample ID`, Molecular_Response, pls_code) %>%
  summarise(mean_pred = mean(Predictions, na.rm = TRUE), .groups = "drop")

mw <- wilcox.test(mean_pred ~ Molecular_Response,
                  data = tissue_data)
roc_tis <- roc(tissue_data$pls_code, tissue_data$mean_pred, quiet = TRUE)
auc_tis <- as.numeric(auc(roc_tis))

cat(sprintf("Tissue-level: U = %.1f, p = %.3f, AUC = %.3f\n",
            mw$statistic, mw$p.value, auc_tis))

p_D <- ggplot(tissue_data,
              aes(Molecular_Response, mean_pred,
                  fill = Molecular_Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.12, color = "black", size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Nonresponder" = "lightblue",
                               "Responder"    = "lightpink")) +
  annotate("text", x = 1.5,
           y = max(tissue_data$mean_pred) + 0.10,
           label = sprintf("p = %.3f", mw$p.value), size = 3.5) +
  labs(title = "Tissue-level (mean per Sample ID)",
       x = "Molecular Response Group", y = "Tissue-level Prediction") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

## ---- Panel E: tissue-level ROC -------------------------------------
p_E <- roc_panel(roc_tis, "ROC Curve (tissue level)")

## ---- Compose and write editable PDF --------------------------------
top    <- p_A
middle <- p_B | p_C
bottom <- p_D | p_E
fig2 <- (top / middle / bottom) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))
## If you did not install cairo package, replace "cairo_pdf" with "pdf"
cairo_pdf(OUT_PDF, width = 9.0, height = 11.5, family = "Helvetica")
print(fig2)
dev.off()
cat("Wrote ", OUT_PDF, "\n", sep = "")

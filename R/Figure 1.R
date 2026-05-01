######################################################################
## Figure 1 - PLS-R analysis of CD8A in immune-rich GBM ROIs
##
## (A) Observed vs. predicted CD8A from a 10-component PLS-R model
## (B) MSEP curve vs. number of components (cross-validated)
## (C) VIP score vs. PLS coefficient for each protein
##
## Data: Raw data.xlsx, sheet "GBM" (168 ROIs from Lu et al., 2021).
## Filter: keep ROIs with Ratio of CD45+ cells >= 0.25 (immune rich).
## Y: CD8A. X: all measured proteins except CD8A and CD3.
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(pls)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

## ---- Standard PLS VIP helper ---------------------------------------
## VIP_j = sqrt(p * sum_a (SSY_a * w_aj^2 / sum_k w_ak^2) / sum_a SSY_a)
## (Wold et al., the textbook PLS VIP). With this convention
## mean(VIP^2) = 1, so VIP > 1 is the conventional "above-average
## importance" cutoff.
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

## ---- Paths (relative to repo root, resolved via here::here) --------
DATA    <- here::here("Nat commun data", "Raw data.xlsx")
OUT_PDF <- here::here("figures", "Figure 1.pdf")

## ---- 1. Load and filter to immune-rich (>25% CD45+) ----------------
gbm <- read_excel(DATA, sheet = "GBM")
gbm_ir <- gbm %>% filter(`Ratio of CD45+ cells` >= 0.25)

X <- gbm_ir %>% select(-c(1:5)) %>% select(-CD8A, -CD3)
Y <- gbm_ir$CD8A

## ---- 2. Fit PLS-R with leave-one-out CV ----------------------------
ncomp_max <- 20
plsr.mod  <- plsr(Y ~ ., data = data.frame(Y = Y, X), ncomp = ncomp_max,
                  validation = "CV", scale = TRUE)

ncomp_use <- 10
preds <- as.numeric(predict(plsr.mod, ncomp = ncomp_use))
cat(sprintf("R2 (Y, ncomp=%d) = %.3f\n",
            ncomp_use,
            1 - sum((Y - preds)^2) / sum((Y - mean(Y))^2)))

## ---- 3. Panel A: observed vs. predicted ----------------------------
df_pred <- data.frame(
  Observed  = Y,
  Predicted = preds,
  Response  = factor(gbm_ir$Molecular_Response,
                     levels = c("Responder", "Nonresponder"))
)

p_A <- ggplot(df_pred, aes(Observed, Predicted, shape = Response)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.6) +
  geom_point(size = 2.2, stroke = 0.8) +
  scale_shape_manual(values = c("Responder" = 3, "Nonresponder" = 1)) +
  labs(title = "GBM-Pembro - >25% CD45+ cells ROIs",
       x = "Observed Values", y = "Predicted Values") +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.12, 0.88),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(hjust = 0.5))

## ---- 4. Panel B: MSEP curve ----------------------------------------
msep_obj <- MSEP(plsr.mod, estimate = c("train", "CV"))
msep_df  <- data.frame(
  ncomp = rep(0:ncomp_max, 2),
  MSEP  = c(as.numeric(msep_obj$val["train", , ]),
            as.numeric(msep_obj$val["CV",    , ])),
  type  = rep(c("train", "CV"), each = ncomp_max + 1)
)

p_B <- ggplot(msep_df, aes(ncomp, MSEP, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("train" = "black", "CV" = "red")) +
  scale_linetype_manual(values = c("train" = "solid", "CV" = "dashed")) +
  labs(x = "number of components", y = "MSEP") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

## ---- 5. Panel C: VIP vs. coefficient -------------------------------
## VIP via the shared `pls_vip` helper at the top of this file
## (identical helper in Figures 1-4) evaluated at ncomp_use = 10.
vip_vals <- pls_vip(plsr.mod, ncomp_use)
coefs    <- coef(plsr.mod, ncomp = ncomp_use)

dat_C <- data.frame(
  Gene       = colnames(X),
  Importance = vip_vals,
  Coef       = as.numeric(coefs)
)
dat_C$above <- dat_C$Importance > 1.0

p_C <- ggplot(dat_C, aes(Coef, Importance)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0,   linetype = "dashed", color = "black") +
  geom_point(data = subset(dat_C, !above), color = "grey60", size = 1.8) +
  geom_point(data = subset(dat_C,  above), color = "red",    size = 2.2) +
  geom_text_repel(data = subset(dat_C, above),
                  aes(label = Gene),
                  size = 3, box.padding = 0.35,
                  segment.color = "grey40", segment.size = 0.3,
                  max.overlaps = Inf) +
  labs(title = "VIP Scores vs. Coefficients",
       x = "Coefficient", y = "VIP Score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0))

## ---- 6. Compose and write editable PDF -----------------------------
fig1 <- (p_A / p_B / p_C) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

pdf(OUT_PDF, width = 6.8, height = 11.5, family = "Helvetica")
print(fig1)
dev.off()
cat("Wrote ", OUT_PDF, "\n", sep = "")

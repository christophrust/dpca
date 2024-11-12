data(fredmd)
fredmd <- scale(fredmd)

freqs <- -50:50 / 50 * pi
res <- dpca::dpca(fredmd, freqs = freqs, qsel = TRUE, q = 10)

## eigenvalues
matplot(x = freqs, y = t(Re(res$eig$values)), type = "l")

## q selection
cat(sprintf("Number of selected dynamic components: %s\n", res$HL_select$q))

## sample variability of the criterion S^2_C (Hallin & Liska 2007, equation 10)
pdf("hl_q_graph.pdf", width = 10, height = 7)
par(mar = c(5, 4, 4, 6))
plot(
  x = res$HL_select$penalty_scales, y = res$HL_select$q_path,
  type = "n", col = "blue", ylab = "Number of selected factors",
  xlab = "c", lwd = 2, xlim = c(0, 1)
)
par(new = TRUE)
plot(
  x = res$HL_select$penalty_scales, y = res$HL_select$sample_var, type = "l",
  axes = FALSE, bty = "n", xlab = "", ylab = "",
  col = "red", lwd = 1.5, xlim = c(0, 1)
)
axis(4)
mtext(expression(paste("Sample variance of ", S[c])), side = 4, padj = 3)
par(new = TRUE)
plot(
  x = res$HL_select$penalty_scales, y = res$HL_select$q_path,
  axes = FALSE, bty = "n", xlab = "", ylab = "",
  type = "l", col = "blue", lwd = 3, xlim = c(0, 1)
)
dev.off()


r_selection <- select_r(x = fredmd, crit = "IC2", max_r = 20)

pdf("hl_r_graph.pdf", width = 10, height = 7)
par(mar = c(5, 4, 4, 6))
plot(
  x = r_selection$penalty_scales, y = r_selection$q_path,
  type = "n", ylab = "Number of selected factors",
  xlab = "c", lwd = 2
)
par(new = TRUE)
plot(
  x = r_selection$penalty_scales, y = r_selection$sample_var, type = "l",
  axes = FALSE, bty = "n", xlab = "", ylab = "",
  col = "red", lwd = 1.5
)
axis(4)
mtext(expression(paste("Sample variance of ", S[c])), side = 4, padj = 3)
par(new = TRUE)
plot(
  x = r_selection$penalty_scales, y = r_selection$q_path,
  axes = FALSE, bty = "n", xlab = "", ylab = "",
  type = "l", col = "blue", lwd = 3
)
dev.off()




## common component vs. observed data
pdf("dcc_vs_observed.pdf", width = 10, height = 7)
par(mfrow = c(1, 2))
indx <- which(colnames(fredmd) == "INDPRO")
matplot(
  y = t(rbind(t(fredmd[, indx]), res$dcc[indx, ])), x = as.numeric(time(fredmd)),
  t = "l", col = c("red", "blue"), main = "INDPRO", ylab = "", xlab = "Time",
  ylim = c(-5, 5), lty = 1
)
indx <- which(colnames(fredmd) == "USTRADE")
matplot(
  y = t(rbind(t(fredmd[, indx]), res$dcc[indx, ])), x = as.numeric(time(fredmd)),
  t = "l", col = c("red", "blue"), ylab = "", main = "USTRADE", xlab = "Time",
  ylim = c(-5, 5), lty = 1
)
dev.off()

## common component vs. observed data
pdf("scc_vs_observed.pdf", width = 10, height = 7)
par(mfrow = c(1, 2))
indx <- which(colnames(fredmd) == "INDPRO")
matplot(
  y = t(rbind(t(fredmd[, indx]), res_spca$cc[indx, ])), x = as.numeric(time(fredmd)),
  t = "l", col = c("red", "blue"), main = "INDPRO", ylab = "", xlab = "Time",
  ylim = c(-5, 5), lty = 1
)
indx <- which(colnames(fredmd) == "USTRADE")
matplot(
  y = t(rbind(t(fredmd[, indx]), res_spca$cc[indx, ])), x = as.numeric(time(fredmd)),
  t = "l", col = c("red", "blue"), ylab = "", main = "USTRADE", xlab = "Time",
  ylim = c(-5, 5), lty = 1
)
dev.off()

su_dpca <- summary(lm(as.vector(t(fredmd)) ~ as.vector(res$dcc)))
su_spca <- summary(lm(as.vector(t(fredmd)) ~ as.vector(res_spca$cc)))

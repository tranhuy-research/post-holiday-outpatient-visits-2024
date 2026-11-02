# OUTPATIENT VISITS AND PUBLIC HOLIDAYS, 2024
options(width = 120)
input <- "outpatient_visits_2024.csv"

# 1. Read data and validate the daily series
d <- read.csv(input, colClasses = c(date = "character", visits = "numeric"))
stopifnot(identical(names(d), c("date", "visits")))
raw_date <- d$date
d$date <- as.Date(raw_date, format = "%Y-%m-%d")
stopifnot(!anyNA(d), all(format(d$date, "%Y-%m-%d") == raw_date),
          !anyDuplicated(d$date), all(is.finite(d$visits)),
          all(d$visits >= 0 & d$visits == floor(d$visits)))
d <- d[order(d$date), ]
rownames(d) <- NULL
calendar <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
stopifnot(nrow(d) == 366L, identical(d$date, calendar))
d$day_id <- seq_len(nrow(d))
d$dow <- factor(as.integer(format(d$date, "%u")), levels = c(7, 1:6),
                labels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
d$is_weekend <- as.integer(d$dow %in% c("Sat", "Sun"))

# 2. Define holiday and exposure windows
holidays <- data.frame(
  id = 1:5,
  holiday = c("New Year", "Tet", "Hung Kings", "30 Apr - 1 May", "National Day"),
  start = as.Date(c("2024-01-01", "2024-02-08", "2024-04-18", "2024-04-27", "2024-08-31")),
  end = as.Date(c("2024-01-01", "2024-02-14", "2024-04-18", "2024-05-01", "2024-09-03")))
pre <- paste0("pre_m", 3:1)
for (h in 3:1) d[[paste0("pre_m", h)]] <- as.integer(d$date %in% (holidays$start - h))
d$post_k <- 0L
for (i in 1:5) {
  d[[paste0("H", i)]] <- as.integer(d$date >= holidays$start[i] & d$date <= holidays$end[i])
  k <- as.integer(d$date - holidays$end[i])
  d[[paste0("P", i)]] <- as.integer(k %in% 1:7)
  d$post_k[k %in% 1:7] <- k[k %in% 1:7]
}
d$is_holiday <- rowSums(d[paste0("H", 1:5)])
d$post_on <- rowSums(d[paste0("P", 1:5)])
d$post_slope <- ifelse(d$post_on == 1, d$post_k - 1L, 0L)
is_pre <- rowSums(d[pre]) > 0
d$reference <- !is_pre & d$is_holiday == 0 & d$post_on == 0
qc <- c(colSums(d[pre]), is_holiday = sum(d$is_holiday), post_on = sum(d$post_on),
        pre_post_overlap = sum(is_pre & d$post_on == 1), reference = sum(d$reference))
stopifnot(all(qc == c(4, 4, 4, 18, 35, 2, 303)),
          all(d$is_holiday <= 1), all(d$post_on <= 1),
          identical(d$date[is_pre & d$post_on == 1], as.Date(c("2024-04-24", "2024-04-25"))))
print(qc)

# 3. Summarise daily outpatient visits
describe <- function(x) c(n = length(x), mean = mean(x), SD = sd(x),
                         median = median(x), Q1 = unname(quantile(x, .25)),
                         Q3 = unname(quantile(x, .75)))
groups <- c(list("All days" = rep(TRUE, nrow(d)), Reference = d$reference),
            setNames(lapply(pre, function(v) d[[v]] == 1), paste0("Pre -", 3:1)),
            list(Holiday = d$is_holiday == 1),
            setNames(lapply(1:7, function(k) d$post_k == k), paste0("Post +", 1:7)))
table1 <- data.frame(group = names(groups),
                     do.call(rbind, lapply(groups, function(g) describe(d$visits[g]))),
                     row.names = NULL)
table1_display <- data.frame(
  Window = table1$group, Days = as.integer(table1$n),
  `Mean (SD)` = sprintf("%.0f (%.0f)", table1$mean, table1$SD),
  `Median [Q1-Q3]` = sprintf("%.0f [%.0f-%.0f]", table1$median, table1$Q1, table1$Q3),
  check.names = FALSE)
show_table1 <- function() {
  cat("\nTABLE 1. Daily outpatient visits by holiday window, 2024\n",
      "Unit: visits/day; sample size: calendar days.\n\n", sep = "")
  print(table1_display, row.names = FALSE, right = FALSE)
  cat("\nReference: outside all Pre, Holiday and Post windows.\n",
      "Pre: days before holiday start; Post: days after holiday end.\n",
      "24-25 Apr belong to both Pre and Post; rows are not mutually exclusive.\n",
      "New Year Pre days fall outside 2024. Descriptive values are unadjusted.\n", sep = "")
  invisible(table1_display)
}
show_table1()
# Descriptive statistics by weekday.
table_dow <- data.frame(dow = levels(d$dow),
  do.call(rbind, lapply(levels(d$dow), function(v) describe(d$visits[d$dow == v]))),
  row.names = NULL)
table_weekend <- data.frame(
  group = c("Weekdays", "Weekend"),
  do.call(rbind, lapply(0:1, function(v) describe(d$visits[d$is_weekend == v]))),
  row.names = NULL)
print(table_weekend, digits = 4, row.names = FALSE)

# 4. Define temporal splines, models and contrasts
# Keep full-year spline knots fixed across leave-one-holiday-out analyses.
for (df in c(5, 7, 9)) {
  nm <- paste0("t", df, "_", seq_len(df))
  d[nm] <- as.data.frame(splines::ns(d$day_id, df = df))
}
fit_model <- function(data = d, df = 7, exposure = c("is_holiday", "post_on", "post_slope")) {
  f <- reformulate(c(pre, exposure, "dow", paste0("t", df, "_", seq_len(df))), "visits")
  m <- glm(f, data = data, family = quasipoisson(link = "log"), na.action = na.fail)
  stopifnot(m$converged, !anyNA(coef(m)), m$rank == length(coef(m)))
  m
}

# Calendar-time HAC covariance with Bartlett weights and maximum lag 14.
# Quasi-Poisson dispersion cancels in the sandwich covariance.
# Use zero scores for excluded dates to preserve calendar-day spacing.
hac_lag <- 14L
vcov_calendar <- function(m, day_id, lag = hac_lag) {
  X <- model.matrix(m)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(day_id) == n, !anyDuplicated(day_id))
  U <- matrix(0, nrow = 366, ncol = p)
  U[day_id, ] <- X * as.numeric(m$y - fitted(m))
  meat <- crossprod(U)
  for (j in seq_len(lag)) {
    a <- crossprod(U[(j + 1):366, , drop = FALSE], U[1:(366 - j), , drop = FALSE])
    meat <- meat + (1 - j / (lag + 1)) * (a + t(a))
  }
  bread <- solve(crossprod(X, X * as.numeric(fitted(m))))
  V <- n / (n - p) * bread %*% meat %*% bread
  dimnames(V) <- list(colnames(X), colnames(X))
  (V + t(V)) / 2
}
contrast <- function(m, V, L, label) {
  est <- drop(L %*% coef(m))
  se <- sqrt(pmax(0, rowSums((L %*% V) * L)))
  data.frame(term = label, RR = exp(est), lower = exp(est - qnorm(.975) * se),
             upper = exp(est + qnorm(.975) * se), SE_logRR = se)
}
unit_L <- function(m, terms) {
  L <- matrix(0, length(terms), length(coef(m)), dimnames = list(NULL, names(coef(m))))
  L[cbind(seq_along(terms), match(terms, colnames(L)))] <- 1
  L
}
post_L <- function(m) {
  L <- unit_L(m, rep("post_on", 7))
  L[, "post_slope"] <- 0:6
  L
}
rr_m3 <- function(m, V) {
  L <- rbind(unit_L(m, c(pre, "is_holiday")), post_L(m))
  contrast(m, V, L, c("Pre -3", "Pre -2", "Pre -1", "Holiday", paste0("Post +", 1:7)))
}

# 5. Fit the primary model and compare HAC with ordinary covariance
m3 <- fit_model()
V3 <- vcov_calendar(m3, d$day_id)
table2 <- rr_m3(m3, V3)
table2_ordinary <- rr_m3(m3, vcov(m3))
print(table2, digits = 4)
se_comparison <- data.frame(term = table2$term, SE_quasi = table2_ordinary$SE_logRR,
                            SE_HAC = table2$SE_logRR,
                            SE_ratio = table2$SE_logRR / table2_ordinary$SE_logRR,
                            lower_quasi = table2_ordinary$lower,
                            upper_quasi = table2_ordinary$upper)
print(se_comparison, digits = 4)
Lpre <- unit_L(m3, pre)
bpre <- drop(Lpre %*% coef(m3))
W <- drop(t(bpre) %*% solve(Lpre %*% V3 %*% t(Lpre), bpre))
pre_joint <- data.frame(test = "Joint Pre: HAC Wald", chi_square = W, df = 3,
                        p = pchisq(W, df = 3, lower.tail = FALSE))
print(pre_joint)

# Sum absolute differences across all 35 post-holiday dates.
# Keep other covariates fixed and use HAC covariance for delta-method intervals.
idx <- which(d$post_on == 1)
X <- model.matrix(m3)[idx, , drop = FALSE]
X0 <- X
X0[, c("post_on", "post_slope")] <- 0
mu <- exp(drop(X %*% coef(m3)))
mu0 <- exp(drop(X0 %*% coef(m3)))
excess <- sum(mu - mu0)
gradient <- colSums(X * mu - X0 * mu0)
excess_se <- sqrt(drop(t(gradient) %*% V3 %*% gradient))
absolute_post <- data.frame(days = length(idx), excess = excess, SE_HAC = excess_se,
                            lower = excess - qnorm(.975) * excess_se,
                            upper = excess + qnorm(.975) * excess_se)
print(absolute_post, digits = 5)
cat("Model-based contrasts describe differences in expected visit counts.\n")
# Estimate absolute differences separately for post-holiday days +1 to +7.
# Each estimate sums five dates, one from each holiday period.
post_gradients <- t(vapply(1:7, function(k) {
  j <- which(d$post_k[idx] == k)
  colSums(X[j, , drop = FALSE] * mu[j] - X0[j, , drop = FALSE] * mu0[j])
}, numeric(ncol(X))))
absolute_post_cov <- post_gradients %*% V3 %*% t(post_gradients)
absolute_post_daily <- data.frame(
  day = 1:7,
  days = vapply(1:7, function(k) sum(d$post_k[idx] == k), integer(1)),
  difference = vapply(1:7, function(k) sum((mu - mu0)[d$post_k[idx] == k]), numeric(1)),
  SE_HAC = sqrt(pmax(0, diag(absolute_post_cov))))
absolute_post_daily$lower <- with(absolute_post_daily, difference - qnorm(.975) * SE_HAC)
absolute_post_daily$upper <- with(absolute_post_daily, difference + qnorm(.975) * SE_HAC)
# Use cross-day covariance to calculate the confidence interval for the total.
stopifnot(all(absolute_post_daily$days == 5),
          isTRUE(all.equal(sum(absolute_post_daily$difference), excess)),
          isTRUE(all.equal(unname(colSums(post_gradients)), unname(gradient))),
          isTRUE(all.equal(sqrt(sum(absolute_post_cov)), excess_se)))
print(absolute_post_daily, digits = 4, row.names = FALSE)

# Print the figure legend separately from the plot.
figure1_legend <- paste(
  "Figure 1. Relative and absolute post-holiday differences in recorded outpatient visits, 2024.",
  "Both panels share aligned calendar-day positions +1 to +7 after holiday end.",
  "(a) Adjusted Post visit ratios from M3; vertical bars are pointwise 95% HAC confidence intervals (Bartlett lag 14 days).",
  "RR is plotted on a logarithmic scale. Reference days lie outside all Pre, Holiday and Post windows.",
  "All three Pre indicators and the Holiday indicator remain in M3; only their estimates are omitted from this figure.",
  "The Post line follows the linear response on log(RR) imposed by M3.",
  "(b) Model-estimated absolute differences at each Post day, summed across the five holiday periods (five calendar dates per point).",
  "Contrasts set post_on and post_slope to zero while retaining Pre indicators, weekday and time spline (df=7).",
  "Confidence intervals use the delta method with the same HAC covariance; the total uses all cross-day covariances.",
  "The two overlapping dates (24-25 April) retain both exposure components. These are model contrasts, not confirmed catch-up patients.")

# Draw two aligned panels using estimates from the primary model.
# Open a separate graphics window when the plotting area is too small.
plot_main <- function(new_window = FALSE) {
  if (new_window || (interactive() &&
      (dev.size("in")[1] < 6.5 || dev.size("in")[2] < 6.5 ||
       dev.size("in")[1] / dev.size("in")[2] > 1.3)))
    dev.new(width = 8, height = 7.6, noRStudioGD = TRUE)
  op <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(op) }, add = TRUE)
  ink <- "#29343B"; muted <- "#667780"; post_col <- "#496B73"
  grid_col <- "#E7EBED"; ref_col <- "#8C989E"
  par(family = "sans", col = ink, col.axis = ink, col.lab = ink, fg = ink,
      cex = 1, oma = c(.3, 0, .3, .2), xaxs = "i", yaxs = "i", tcl = -.2)
  layout(matrix(1:2, ncol = 1), heights = c(1, 1))

  # (a) Post-holiday ratios adjusted for pre-holiday and holiday indicators.
  par(mar = c(3.2, 4.5, 3.2, .8), mgp = c(2.7, .7, 0))
  p <- table2[5:11, ]
  xp <- 1:7
  # Use identical horizontal limits and margins to align days across panels.
  lim <- range(c(1, p$lower, p$upper))
  pad <- .10 * diff(log(lim))
  lim <- exp(log(lim) + c(-pad, pad))
  plot(xp, p$RR, type = "n", log = "y", xlim = c(.55, 7.45), ylim = lim,
       axes = FALSE, xlab = "", ylab = "Adjusted RR")
  yt <- sort(unique(c(pretty(lim, n = 4), 1)))
  yt <- yt[yt >= lim[1] & yt <= lim[2]]
  abline(h = yt[yt != 1], col = grid_col, lwd = .6)
  abline(h = 1, col = ref_col, lwd = 1.2, lty = 2)
  axis(2, at = yt, labels = format(yt, trim = TRUE, scientific = FALSE),
       las = 1, lwd = 0, lwd.ticks = .7, cex.axis = .9)
  axis(1, at = xp, labels = paste0("+", 1:7),
       lwd = 0, lwd.ticks = .7, cex.axis = .88)
  mtext("Days after holiday end", 1, line = 2.0, cex = .9)
  mtext("a  Adjusted visit ratios", 3, line = 1.8, adj = 0, font = 2, cex = 1.06)
  lines(xp, p$RR, col = post_col, lwd = 2)
  arrows(xp, p$lower, xp, p$upper, angle = 90, code = 3,
         length = .045, col = post_col, lwd = 1.5)
  points(xp, p$RR, pch = 21, bg = post_col, col = post_col, cex = 1.25, lwd = 1.5)

  # (b) Daily absolute differences with confidence intervals.
  par(mar = c(3.2, 4.5, 3.2, .8), mgp = c(2.7, .7, 0))
  a <- absolute_post_daily
  yr <- range(c(0, a$lower, a$upper))
  yr <- yr + c(-.08, .13) * diff(yr)
  plot(a$day, a$difference, type = "n", xlim = c(.55, 7.45), ylim = yr,
       axes = FALSE, xlab = "", ylab = "Difference in recorded visits")
  yt <- pretty(yr, n = 4); yt <- yt[yt >= yr[1] & yt <= yr[2]]
  abline(h = yt[yt != 0], col = grid_col, lwd = .6)
  abline(h = 0, col = ref_col, lty = 2, lwd = 1.2)
  axis(2, at = yt, labels = format(yt, big.mark = ",", trim = TRUE, scientific = FALSE),
       las = 1, lwd = 0, lwd.ticks = .7, cex.axis = .9)
  axis(1, at = 1:7, labels = paste0("+", 1:7), lwd = 0, lwd.ticks = .7, cex.axis = .9)
  arrows(a$day, a$lower, a$day, a$upper, angle = 90, code = 3,
         length = .045, col = post_col, lwd = 1.5)
  points(a$day, a$difference, pch = 21, col = post_col, bg = post_col, cex = 1.25)
  mtext("Days after holiday end (summed across 5 holiday periods)", 1, line = 2.0, cex = .9)
  mtext("b  Absolute post-holiday differences", 3, line = 1.8, adj = 0, font = 2, cex = 1.06)
  invisible(list(relative = p, absolute_daily = absolute_post_daily,
                 absolute_total = absolute_post))
}

# 6. Examine model diagnostics
r <- residuals(m3, type = "pearson")
diagnostics <- data.frame(converged = m3$converged, n = nobs(m3),
                           residual_df = df.residual(m3),
                           phi_Pearson = sum(r^2) / df.residual(m3), HAC_lag = hac_lag)
print(diagnostics)
ac <- acf(r, lag.max = 30, plot = FALSE)
acf_table <- data.frame(lag = 1:30, ACF = as.numeric(ac$acf)[-1])
print(acf_table, digits = 3)
# Use Ljung-Box tests as approximate diagnostics for GLM residuals.
print(Box.test(r, lag = 14, type = "Ljung-Box"))
print(Box.test(r, lag = 28, type = "Ljung-Box"))
if (any(abs(acf_table$ACF[15:30]) > 1.96 / sqrt(nrow(d))))
  warning("Residual autocorrelation beyond lag 14 exceeds the reference bounds; assess HAC bandwidth adequacy.")
influence_table <- data.frame(date = d$date, visits = d$visits, Pearson = r,
                               leverage = hatvalues(m3), Cook = cooks.distance(m3))
print(head(influence_table[order(influence_table$Cook, decreasing = TRUE), ], 10), digits = 4)
plot_diagnostics <- function() {
  if (interactive() && (dev.size("in")[1] < 5 || dev.size("in")[2] < 4))
    dev.new(width = 7, height = 5, noRStudioGD = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
  plot(d$date, r, type = "l", xlab = "Date", ylab = "Pearson residual"); abline(h = 0, lty = 2)
  plot(fitted(m3), r, xlab = "Fitted visits", ylab = "Pearson residual"); abline(h = 0, lty = 2)
  acf(r, lag.max = 30, main = "Pearson residual ACF")
  plot(d$date, influence_table$Cook, type = "h", xlab = "Date", ylab = "Cook distance")
  invisible(influence_table)
}
if (interactive()) plot_diagnostics()

# 7. Explore holiday-specific associations with a shared post-holiday slope
m4 <- fit_model(exposure = c(paste0("H", 1:5), paste0("P", 1:5), "post_slope"))
V4 <- vcov_calendar(m4, d$day_id)
table_m4 <- cbind(holiday = rep(holidays$holiday, 2),
                  contrast(m4, V4, unit_L(m4, c(paste0("H", 1:5), paste0("P", 1:5))),
                           c(paste0("Holiday H", 1:5), paste0("Post +1 P", 1:5))))
print(table_m4, digits = 4)
cat("M4 estimates refer to the five observed holiday periods in 2024.\n")

# 8. Run sensitivity analyses S1-S3
sensitivity <- function(m, data, label) {
  out <- rr_m3(m, vcov_calendar(m, data$day_id))
  data.frame(analysis = label, n = nrow(data), out)
}
s1 <- do.call(rbind, lapply(c(5, 9), function(df)
  sensitivity(fit_model(df = df), d, paste0("S1: time df=", df))))

# S2: Use three post-holiday spline parameters, including the intercept.
# Construct the spline within days +1 to +7 and set it to zero elsewhere.
Bpost <- splines::ns(1:7, df = 3, intercept = TRUE, Boundary.knots = c(1, 7))
post_names <- paste0("post_ns", 1:3)
d[post_names] <- 0
d[idx, post_names] <- Bpost[d$post_k[idx], ]
m_s2 <- fit_model(exposure = c("is_holiday", post_names))
L2 <- matrix(0, 7, length(coef(m_s2)), dimnames = list(NULL, names(coef(m_s2))))
L2[, post_names] <- Bpost
s2 <- contrast(m_s2, vcov_calendar(m_s2, d$day_id), L2, paste0("Post +", 1:7))
s2_comparison <- data.frame(day = 1:7, RR_linear = table2$RR[5:11],
                            RR_spline = s2$RR, lower_spline = s2$lower, upper_spline = s2$upper)

# S3: Exclude each complete pre-holiday, holiday and post-holiday window.
# April exclusions also remove dates in overlapping exposure windows.
s3_models <- vector("list", 5)
s3 <- do.call(rbind, lapply(1:5, function(i) {
  keep <- !(d$date >= holidays$start[i] - 3 & d$date <= holidays$end[i] + 7)
  di <- d[keep, ]
  m <- fit_model(data = di)
  s3_models[[i]] <<- m
  sensitivity(m, di, paste0("S3: omit ", holidays$holiday[i]))
}))
print(s1, digits = 4, row.names = FALSE)
print(s2_comparison, digits = 4, row.names = FALSE)
print(s3, digits = 4, row.names = FALSE)
cat("The primary model assumes a linear post-holiday trajectory on the log-RR scale.\n",
    "Only five holiday periods were observed at one hospital in one year; intervals may underrepresent between-period uncertainty.\n")

# Display Table 1 and the primary figure after running the script.
# Inspect table_dow for weekday summaries; call plot_main(new_window = TRUE) to redraw.
show_table1()
cat("\nFIGURE 1 LEGEND\n", figure1_legend, "\n", sep = "")
if (interactive()) plot_main()


# Sample size for testing a correlation using Fisher's r-to-z transform
# From: https://homepage.univie.ac.at/robin.ristl/samplesize.php?test=correlation
# Two-sided test + Bonferroni correction + drop-out adjustment

r_true   <- 0.40   # expected true correlation
r_null   <- 0.00   # correlation under H0
alpha    <- 0.05   # nominal alpha
m_tests  <- 2      # Bonferroni divisor (number of tests)
power    <- 0.80
dropout  <- 0.15

# Bonferroni-corrected alpha (two-sided test)
alpha_adj <- alpha / m_tests

# AL v.0.2

# Sample size calculation

# Fisher z transform
z_true <- atanh(r_true)
z_null <- atanh(r_null)
delta  <- abs(z_true - z_null)

# Normal quantiles
z_alpha <- qnorm(1 - alpha_adj/2)  # two-sided
z_beta  <- qnorm(power)

# Required sample size (approx.)
n_no_dropout <- ceiling(3 + ((z_alpha + z_beta) / delta)^2)

# Inflate for expected drop-out
n_with_dropout <- ceiling(n_no_dropout / (1 - dropout))

cat("alpha_adj =", alpha_adj, "\n")
cat("n (no dropout) =", n_no_dropout, "\n")
cat("n (with 15% dropout) =", n_with_dropout, "\n")


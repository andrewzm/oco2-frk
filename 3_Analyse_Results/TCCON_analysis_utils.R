## Compute variance even with one obs in sample (variance of one obs is zero)
safe.var <- function(x)
  if(length(x) == 1) 0 else var(x)

## Convert the first letter in string to upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x
}

## Find the 95% empirical coverage from data (d) and prediction distribution
cov95 <- function(d,mu,sd) {
  mean((d > mu - 1.96*sd) & (d < mu + 1.96*sd))
}



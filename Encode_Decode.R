# Implemnation in R of Laurent Boué's article, "Real numbers, data science
# and chaos: How to fit any dataset with a single parameter". See
# arXiv: 1904.12320v1

library(pracma)
library(stringr)
library(Rmpfr)
library(purrr)
library(ggplot2)

#
# Functions
#

# page 8
dyadicMap <- function(x) {
  (2 * x) %% 1
}

# page 17
decimalToBinary <- function(x, tau = 8) {
  # added 0 in case initial value is 0
  # added 1 to more easily handle min-max standardization
  if (x == 1) {
    d2b <- paste0('0.', paste(rep('1', tau), collapse = ''))
  } else if (x == 0) {
    d2b <- paste0('0.', paste(rep('0', tau - 1), collapse = ''), '1')
  } else {
    d2b <- bits(x, k = tau)
  }
  str_sub(d2b, 3,-1)
}

# page 17
binaryToDecimal <- function(x) {
  #bits <- strsplit(x, '') %>% unlist()
  #b2d <- mpfr(sum(2 ^ -(which(bits == 1))), precBits = str_length(x))
  #return(b2d)
  mpfr(paste0('0.', x), precBits = nchar(x), base = 2)
}

# page 9
dyadicDecoder <- function(k, tau = 8, di) {
  dd <- 2 ^ ((k - 1) * tau) * di
  as.numeric(dd - floor(dd))
}

# page 14
phiInv <- function(z) {
  asin(sqrt(z)) / (2 * pi)
}

# page 14
decimalToBinary_phiInv <- function(z, tau = 8) {
  decimalToBinary(phiInv(z), tau = tau)
}

# page 15
phi <- function(alpha) {
  big_pi <- Const('pi', prec = getPrec(alpha))
  sin(2 * big_pi * alpha) ^ 2
}

# page 15
logisticDecoder <- function(k, tau = 8, di) {
  ld <- sin(2 ^ ((k - 1) * tau) * asin(sqrt(di))) ^ 2
  as.numeric(ld)
}

# A taylor series to calculate sine
sin_t <- function(x) {
  # makes series same data type as x
  series <- x * 0
  last_val <- x * 0
  n <- 0
  repeat {
    # taylor series
    # do these steps in pairs: one positive, one negative
    step <- (2 * n) + 1
    series <- series + ((-1) ^ n) * (x ^ step) / fact(step)
    
    n <- n + 1
    step <- (2 * n) + 1
    series <- series + ((-1) ^ n) * (x ^ step) / fact(step)
    
    # look to see if series is changing, if not we can end expansion
    if (last_val == series)
      break
    #if (n > 15) break
    last_val <- series
    n <- n + 1
  }
  #attr(series, 'n') <- n
  series
}

sin_p <- function(x) {
  # WARNING: This function takes *much* longer to converge than sin_t
  # makes series same data type as x
  series <- x ^ 0
  last_val <- x * 0
  n <- 1
  repeat {
    series <- series * (1 - ((x ^ 2) / (n * pi) ^ 2))
    if ((x * last_val) == (x * series))
      break
    #if (n > 1000) break
    last_val <- series
    n <- n + 1
  }
  #attr(series, 'n') <- n
  x * series
}

# a taylor series to calculate arcsine
asin_t <- function(x) {
  # makes series same data type as x
  series <- x * 0
  last_val <- x * 0
  n <- 0
  repeat {
    # taylor series
    step <- (2 * n) + 1
    num <- fact(2 * n) * (x ^ step)
    denom <- (2 ^ (2 * n)) * (fact(n) ^ 2) * (step)
    series <- series + num / denom
    # look to see if series is changing, if not we can end expansion
    if (last_val == series)
      break
    #if (n > 20) break
    last_val <- series
    n <- n + 1
  }
  #attr(series, 'n') <- n
  series
}

# Nilicantha's formula for pi
# # pi = 3 + 4/(2*3*4) – 4/(4*5*6) + 4/(6*7*8) – 4/(8*9*10) + 4/(10*11*12) – (4/(12*13*14) …
pi_calc <- function(digits = 22) {
  # TODO: Need to create a mfpr data type if I want to go big
  #options(digits = 22) # 22 is max in double-precision, for printing only
  pi_calc <- 3
  x <- 2
  n <- 0
  last_val <- 0
  
  repeat {
    n <- n + 1
    y <- x + 1
    z <- x + 2
    combine_with_pi_calc <- 4 / (x * y * z)
    if (n %% 2 == 0) {
      pi_calc <- pi_calc - combine_with_pi_calc
    } else {
      pi_calc <- pi_calc + combine_with_pi_calc
    }
    #print(paste(pi_calc, last_val))
    if (last_val == pi_calc)
      break
    if (n == 10)
      break
    last_val <- pi_calc
    x <- z
    #print(pi_calc)
  }
  pi_calc
  #attr(pi_calc, 'x') <- x
  #attr(pi_calc, 'n') <- n
  
  return(pi_calc)
  
}


#
# Examples
#

xs <- runif(50, 0, 1)
xs_backup <- xs
xs[2:3] <- 0

xs <-
  c(
    0.54881350,
    0.71518937,
    0.60276338,
    0.54488318,
    0.42365480,
    0.64589411,
    0.43758721,
    0.89177300,
    0.96366276,
    0.38344152,
    0.79172504,
    0.52889492,
    0.56804456,
    0.92559664,
    0.07103606,
    0.08712930,
    0.02021840,
    0.83261985,
    0.77815675,
    0.87001215,
    0.97861834,
    0.79915856,
    0.46147936,
    0.78052918,
    0.11827443,
    0.63992102,
    0.14335329,
    0.94466892,
    0.52184832,
    0.41466194,
    0.26455561,
    0.77423369,
    0.45615033,
    0.56843395,
    0.01878980,
    0.61763550,
    0.61209572,
    0.61693400,
    0.94374808,
    0.68182030,
    0.35950790,
    0.43703195,
    0.69763120,
    0.06022547,
    0.66676672,
    0.67063787,
    0.21038256,
    0.12892630,
    0.31542835,
    0.36371077
  )

# t is tau
t <- 8
n <- length(xs)

# # #
# Double-precision values
# All R platforms are required to work with a precision of 53 bits.
precision_binary <- n * t
(precision_decimal <-
    as.integer(log(2) / log(10) * precision_binary))
# # #

# Dyadic Decoder
# page 8
if (TRUE) {
  binaryInitial <-
    map(xs, decimalToBinary, tau = t) %>% paste(collapse = '')
  print(binaryInitial)
  
  # page 9
  decimalInitial <- binaryToDecimal(binaryInitial)
  
  decodedValues <-
    dyadicDecoder(seq(1:n), tau = t, di = decimalInitial)
  
  maxError <- 2 ^ -t
  
  errors <- abs(decodedValues - xs)
  normalizedErrors <- abs(decodedValues - xs) / maxError
  
  df <- data.frame(
    x = seq(1:n),
    initial_Value = xs,
    decoded_Value = decodedValues,
    error = errors,
    normalized_Error = normalizedErrors
  )
  
  df
  
  ggplot(df, aes(x = x)) +
    geom_line(aes(y = initial_Value), color = 'orange') +
    geom_point(aes(y = decoded_Value), shape = 'star') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, len = 6)) +
    scale_x_continuous(limits = c(1, n), breaks = seq(1, n, len = n))
  
  ggplot(df, aes(x = x)) +
    geom_point(
      aes(y = normalizedErrors),
      color = 'darkgreen',
      shape = 'diamond',
      size = 2
    ) +
    geom_hline(yintercept = 1,
               color = 'red',
               linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, len = 6)) +
    scale_x_continuous(limits = c(1, n), breaks = as.integer(seq(0, n, len = 6)))
  
}

# Conjugate Decoder
# page 14
if (TRUE) {
  conjugateInitial_binary <-
    map(xs, decimalToBinary_phiInv, tau = t) %>%
    paste(collapse = '')
  print(conjugateInitial_binary)
  
  conjugateInitial <- binaryToDecimal(conjugateInitial_binary)
  conjugateInitial
  
  # Good until here
  
  # page 15
  decimalInitial <- phi(conjugateInitial)
  decimalInitial
  
  decodedValues <-
    logisticDecoder(seq(1:n), tau = t, di = decimalInitial)
  
  maxError <- pi / (2 ^ (t - 1))
  
  errors <- abs(decodedValues - xs)
  normalizedErrors <- abs(decodedValues - xs) / maxError
  
  df <- data.frame(
    x = seq(1:n),
    initial_Value = xs,
    decoded_Value = decodedValues,
    error = errors,
    normalized_Error = normalizedErrors
  )
  
  df
  
  ggplot(df, aes(x = x)) +
    geom_line(aes(y = initial_Value), color = 'orange') +
    geom_point(aes(y = decoded_Value), shape = 'star') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, len = 6)) +
    scale_x_continuous(limits = c(1, n), breaks = seq(1, n, len = n))
  
  ggplot(df, aes(x = x)) +
    geom_point(
      aes(y = normalizedErrors),
      color = 'darkgreen',
      shape = 'diamond',
      size = 2
    ) +
    geom_hline(yintercept = 1,
               color = 'red',
               linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, len = 6)) +
    scale_x_continuous(limits = c(1, n), breaks = as.integer(seq(0, n, len = 6)))
}

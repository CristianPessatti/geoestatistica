require(tidyverse)

dados <- data.frame(
    x1 = c(1,4,4,6,8,9),
    x2 = c(8,1,6,5,3,10),
    y = c(21,18,25,28,30,37)
)

novos_dados <- data.frame(
    x1 = c(5,1,4,7),
    x2 = c(5,7,10,7)
)

require(geoR)

dg <- as.geodata(dados, coords.col=c(1,2), data.col=3)

points(dg)

tal2 = 4
sigma2 = 25
phi = 2

calculate_sigma <- function(x, tal2, sigma2, phi) {
    D <- as.matrix(dist(x))
    
    SIGMA <- (tal2 * diag(nrow(x))) + sigma2 * exp(-D / phi)
    
    result <- list(S = SIGMA, S.t = t(SIGMA), S.inv = solve(SIGMA))
    
    return(result)
}

estimate_mu <- function(SIGMA, y) {
    I <- rep(1, length(y))
    
    # mu = (1'S⁻¹1)⁻¹ (1'S⁻¹y)
    
    IS.D <- crossprod(I, SIGMA$S.inv)
    
    mu_est <- solve(IS.D %*% I) %*% (IS.D %*% y)
    return(as.numeric(mu_est))
}

require(mvtnorm)

predict_geo_data <- function(x, y, newdata, tal2, sigma2, phi) {
    Sy <- calculate_sigma(x, tal2, sigma2, phi)
    
    mu_hat <- estimate_mu(Sy, y)
    
    Y.Yp <- rbind(x, newdata)
    
    S <- calculate_sigma(Y.Yp, tal2, sigma2, phi)
    
    rows_y <- nrow(x)
    rows_yp <- nrow(newdata)
    yp_start <- (nrow(x) + 1)
    tot_rows <- rows_y + rows_yp
    
    Syp <- S$S[yp_start:tot_rows, yp_start:tot_rows]
    Sy.yp <- S$S[1:rows_y, yp_start:tot_rows]
    Syp.y <- t(Sy.yp)
    
    mus_yp <- mu_hat * rep(1, nrow(newdata))
    mus_y <- mu_hat * rep(1, nrow(x))
    
    mu <- mus_yp - Syp.y %*% Sy$S.inv %*% (y - mus_y)
    sig <- Syp - (Syp.y %*% Sy$S.inv %*% Sy.yp)
    
    result <- list(
        MU = mu,
        SIGMA = sig
    )
    
    return(result)
}

predict_geo_data(
    x = dg$coords, 
    y = dg$data,
    newdata = novos_dados,
    tal2 = tal2,
    sigma2 = sigma2,
    phi = phi
)

# PENDENTE: USAR A NORMAL MULTIVARIADA PARA ACHAR OS Y PREDITOS

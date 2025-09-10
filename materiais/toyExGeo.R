X <- rbind(
    c(4,1),
    c(8,3), 
    c(6,5), 
    c(4,6), 
    c(1,8), 
    c(9,10))
y <- c(18, 30, 28, 25, 21, 37)
n <- length(y)

Xp <- rbind(
    c(5,5), 
    c(1,7), 
    c(4,10), 
    c(8,7))

plot(X, pch=19, xlim=c(0,11), ylim=c(0,11))
text(X, labels = as.character(y), pos=3)
text(X, labels = as.character(1:n), pos=1)
abline(h=1:10, v=1:10, col="grey80", lty=3)
points(Xp, pch=19, col = "blue")

t2 <- 4
s2 <- 25
phi <- 3

Oy <- rep(1,6)
dist(X)
args(dist)
(D <- as.matrix(dist(X, diag = T, up = T)))
(S <- t2*diag(6) + s2*exp(-D/phi))

##
## Calculando a mÃ©dia geral de trÃªs formas diferentes
##

solve(t(Oy)%*%solve(S)%*%Oy)%*%(t(Oy)%*%solve(S)%*% y)

iS.Oy <- solve(S, Oy)
(mu <- drop(solve(crossprod(iS.Oy, Oy), crossprod(iS.Oy, y))))
mean(y)

cS <- chol(S)
icS1 <- forwardsolve(t(cS), Oy)
icSy <- forwardsolve(t(cS), y)
sum(icS1*icSy)/sum(icS1^2)
crossprod(icS1,icSy)/crossprod(icS1)
solve(crossprod(icS1),crossprod(icS1,icSy))

##
## PrediÃ§Ã£o em novos pontos
## (esperanÃ§a e variÃ¢ncia da condicional)

Dy.y <- as.matrix(dist(X, diag = TRUE, up = TRUE))
Dyp.yp <- as.matrix(dist(Xp, diag = TRUE, up = TRUE))
Dy.yp <- geoR::loccoords(X, Xp)

Sy.y <- t2*diag(6) + s2*exp(-Dy.y/phi)
Syp.yp <- t2*diag(4) + s2*exp(-Dyp.yp/phi)
Sy.yp <- s2*exp(-Dy.yp/phi)

Oyp <- rep(1, 4)

(ypred <- drop(mu*Oyp + crossprod(Sy.yp, solve(Sy.y, (y-mu*Oy)))))
(Spred <- Syp.yp - crossprod(Sy.yp, solve(Sy.y, Sy.yp)))

text(Xp, labels = sprintf("%.1f", ypred), pos=3, col = "blue")

##
## FunÃ§Ã£o de (log) verossimilhaÃ§a
##
-(n/2)*log(2*pi) - (1/2)*determinant(S, logarithm = TRUE)$modulus - (1/2)*t(y-mu*Oy)%*%solve(S, y-mu*Oy)

-(n/2)*log(2*pi) - (1/2)*determinant(S, logarithm = TRUE)$modulus - (1/2)*mahalanobis(y, mu*Oy, S)

cS <- chol(S)
-(n/2)*log(2*pi) - sum(log(diag(cS))) - (1/2)*t(y-mu*Oy)%*%solve(S, y-mu*Oy)

cS <- chol(S)
z <- forwardsolve(t(cS), y-mu*Oy)
-(n/2)*log(2*pi) - sum(log(diag(cS))) - (1/2)*crossprod(z)

mvtnorm::dmvnorm(y, mean = mu*Oy, sigma = S, log = TRUE)


## comparando tempos
X <- matrix(rnorm(2*N), ncol=2)
y <- rnorm(N)
n <- length(y)
Oy <- rep(1,N)
S <- t2 * diag(N) + s2*exp(-as.matrix(dist(X))/phi)
iS.Oy <- solve(S, Oy)
mu <- drop(solve(crossprod(iS.Oy, Oy), crossprod(iS.Oy, y)))

system.time({
    -(N/2)*log(2*pi) - (1/2)*determinant(S, logarithm = TRUE)$modulus - (1/2)*t(y-mu*Oy)%*%solve(S, y-mu*Oy)
})

system.time(
    -(N/2)*log(2*pi) - (1/2)*determinant(S, logarithm = TRUE)$modulus - (1/2)*mahalanobis(y, mu*Oy, S)
)

system.time({
    cS <- chol(S)
    -(N/2)*log(2*pi) - sum(log(diag(cS))) - (1/2)*t(y-mu*Oy)%*%solve(S, y-mu*Oy)
})

system.time({
    cS <- chol(S)
    z <- forwardsolve(t(cS), y-mu*Oy)
    -(N/2)*log(2*pi) - sum(log(diag(cS))) - (1/2)*crossprod(z)
})

system.time(
    mvtnorm::dmvnorm(y, mean = mu*Oy, sigma = S, log = TRUE)
)


##
dados <- read.table("http://www.leg.ufpr.br/~paulojus/CE097/dados/dadosIntro.txt", head=T)
coords <- dados[,1:2]
y <- dados[,3, drop = TRUE]

## funÃ§Ã£o de (Log)verossimilhanÃ§a (escrita e forma ingenua)
## aqui apenas para modelo exponencial de covariÃ¢ncia
ltheta <- function(pars, coords, y, printpars=FALSE){    
    ## pars = mu, sigma, phi, tau
    n <- nrow(coords)
    D <- as.matrix(dist(coords, upper=T, diag=T))
    SIGMA <- pars[4]^2 * diag(n) + pars[2]^2 * exp(-D/pars[3])
    cholS <- chol(SIGMA)
    res <- y - pars[1]
    z <- forwardsolve(t(cholS), res)
    logVero <- drop(-0.5*(n*log(2*pi) + 2*sum(log(diag(cholS))) + crossprod(z)))
    if(printpars) print(c(pars, logVero))
    return(logVero)
}

ltheta1 <- function(pars, coords, y, printpars=FALSE){    
    ## pars = mu, sigma, phi, tau
    n <- nrow(coords)
    D <- as.matrix(dist(coords, upper=T, diag=T))
    SIGMA <- pars[4]^2 * diag(n) + pars[2]^2 * exp(-D/pars[3])
    logVero <- drop(mvtnorm::dmvnorm(y, mean=rep(pars[1], n), sigma=SIGMA, log=TRUE))
    if(printpars) print(c(pars, logVero))
    return(logVero)
}

lthetaConc <- function(pars, coords, y, printpars=FALSE){    
    ## pars : phi e nu = tau/sigma
    n <- nrow(coords)
    D <- as.matrix(dist(coords, upper=T, diag=T))
    R <- pars[2]^2 * diag(n) + exp(-D/pars[1])
    cholR <- chol(R)
    icholR1 <- forwardsolve(t(cholR), rep(1,n))
    icholRy <- forwardsolve(t(cholR), y)
    mu <- sum(icholR1*icholRy)/sum(icholR1^2)
    res <- y - mu
    icholRres <- forwardsolve(t(cholR), res)
    Q <- crossprod(icholRres)
    sigma2 <- Q/n
    logdetR <-  2*sum(log(diag(cholR)))
    logVero <- drop(-0.5*(n*log(2*pi) + logdetR + n*log(sigma2) + n))
    if(printpars) print(c(pars, logVero))
    attr(logVero, "pars") <- c(mu = mu, sigma = sqrt(sigma2), phi = pars[1], tau = pars[2]*sqrt(sigma2))
    return(logVero)
}

ltheta(c(58.53, 32.1, 20, 5), coords = coords, y = y)
ltheta1(c(58.53, 32.1, 20, 5), coords = coords, y = y)
lthetaConc(c(20, 0), coords = coords, y = y)

lthetaConc(c(10, 0.5), coords = coords, y = y)
ltheta(c(59.63, 19.08, 10, 9.54), coords = coords, y = y)
ltheta1(c(59.63, 19.08, 10, 9.54), coords = coords, y = y)


args(optim)
parsML <- optim(c(50, 15, 20, 3), ltheta, coords = coords, y = y, printpars=FALSE, 
                method="L-BFGS-B", lower=c(-Inf, 0, 0, 0), 
                control=list(fnscale=-1))         
parsML$value
parsML$par
parsML$par^c(1,2,1,2)

parsMLConc <- optim(c(20, 5), lthetaConc, coords = coords, y = y, printpars=FALSE, 
                    method="L-BFGS-B", lower=c(0, 0), 
                    control=list(fnscale=-1))         
parsMLConc$par
parsMLConc$value
lthetaConc(parsMLConc$par,coords = coords, y = y)
attr(lthetaConc(parsMLConc$par, coords = coords, y = y), "pars")^c(1,2,1,2)

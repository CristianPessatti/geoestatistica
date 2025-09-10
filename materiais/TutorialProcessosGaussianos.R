## Processo Gaussiano (e kernels)
## 1. Understanding Gaussian processes
## 2. Fitting a Gaussian process kernel
## 3. Gaussian process kernels

## Acesso por:
## https://peterroelants.github.io/posts/gaussian-process-tutorial/

## RevisÃ£o de normal multivariada
## https://peterroelants.github.io/posts/multivariate-normal-primer/

##
## Parte 1: Understanding Gaussian processes
##
np <- 41                                 # nÃºmero de pontos
x <- as.matrix(seq(-4, 4, length=np))    # pontos no eixo x

mu <- rep(0, np)                     # mÃ©dia zero
D <- as.matrix(dist(x, diag = T, up = T))# matriz de distÃ¢ncias (euclidiana)
S <- exp(-((D^2)/(2*1^2)))                   # matriz de covariÃ¢ncias (kernel RBF/Gaussiano)

## simulando da NMV (1 realizaÃ§Ã£o de np pontos)
cholS <- chol(S)       # decomposiÃ§Ã£o de Cholesky - problema com RBF/Gaussiano
y <- drop(mu + t(cholS)%*%rnorm(np))  # simulaÃ§Ã£o
length(y)
y
plot(x, y, type = "o")
## gerando 5 simulaÃ§Ãµes (cada coluna Ã© uma simulaÃ§Ã£o)
y <- replicate(5, drop(mu + t(cholS)%*%rnorm(np)))
y
dim(y)
matplot(x, y, type = "b")

## Usando pacotes para simular da NMV
##
## MASS::mvrnorm
y <- MASS::mvrnorm(1, mu, S)
plot(x, y, type = "o")
y <- MASS::mvrnorm(5, mu, S)
dim(y)
matplot(x, t(y), type = "o")

## mvtnorm::rmvnorm
y <- mvtnorm::rmvnorm(1, mu, S, method = "chol")
plot(x, y, type = "o")
y <- mvtnorm::rmvnorm(5, mu, S)
dim(y)
matplot(x, t(y), type = "o")

image(S, asp = 1, hcl.colors = 50)
args(image)
?image

tf <- function(m) t(m)[, nrow(m):1]
imageM <- function(m, grid = max(dim(m)) <= 25, asp = (nrow(m)-1)/(ncol(m)-1), ...) {
    image(tf(m), asp=asp, axes = FALSE, ...)
    mAxis <- function(side, at, ...) # using 'j'
        axis(side, at=at, labels=as.character(j+1L), col="gray", col.axis=1, ...)
    n <- ncol(m); n1 <- n-1L; j <- 0L:n1; mAxis(1, at= j/n1)
    if(grid) abline(v = (0:n - .5)/n1, col="gray77", lty="dotted")
    n <- nrow(m); n1 <- n-1L; j <- 0L:n1; mAxis(2, at=1-j/n1, las=1)
    if(grid) abline(h = (0:n - .5)/n1, col="gray77", lty="dotted")
}
imageM(S, asp=1, hcl.colors(50))

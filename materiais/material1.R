##
## Fundamentos de geoestatÃ­stica com o pacote geoR
##

## importando dados
dados <- read.table("http://www.leg.ufpr.br/~paulojus/CE097/dados/dadosIntro.txt", head=T)
dados
head(dados)
str(dados)

## instalando o pacote geoR (se necessÃ¡rio) e suas dependencias
#install.packages("geoR", dep=TRUE)
## carregando o pacote geoR
require(geoR)
dG <- as.geodata(dados, coords.col=c(1,2), data.col=3)
class(dG)

## nomes de outros conjuntos de dados disponÃ­veis no pacote
data(package="geoR")

## visualizando os dados
points(dG)
args(points.geodata)
points.geodata(dG,pt.div="quint")
points.geodata(dG,pt.div="quint", cex.min=1, cex.max=1)

## um resumo dos dados
summary(dG)
summary(dG, lambda=0)  
dG	
## uma outra visualizaÃ§Ã£o dos dados
x11()
plot(dG)
args(plot.geodata)
plot(dG, lambda=0)
plot(dG, lowess=TRUE)

## outras visualizaÃ§Ãµes com points.geodata()
points(dG)
points(dG, lambda=0)
points(dG, pt.div="quint", cex.min=1, cex.max=1)
?points.geodata
points(dG, pt.div="quint", cex.min=1, cex.max=1,permute=TRUE)

## obtendo um variograma (empÃ­rico)
v <- variog(dG)	
plot(v)		

## outras opÃ§Ãµes de variograma
vC <- variog(dG, option="cloud")	
vC <- variog(dG, max.dist=100, option="cloud")	
plot(vC)		

x11()
par(mfrow=c(2,2))
v <- variog(dG, max.dist=100)	
plot(v)		
v <- variog(dG, max.dist=80)	
plot(v)		
v <- variog(dG, uvec=seq(0, 80, by=10))	
plot(v)		
v <- variog(dG, uvec=seq(0, 80, by=5))	
plot(v)		
par(mfrow=c(1,1))

## escolhendo um deles
v <- variog(dG, uvec=seq(0, 80, by=10))	
plot(v)		

## envelope do variograma sob permutaÃ§Ã£o dos dados
v.mc <- variog.mc.env(dG, obj.variog=v)
plot(v, env=v.mc)

## "escolhendo" funcao de variograma (teÃ³rico)
v <- variog(dG, max.dist=80)	
plot(v)		
ef <- eyefit(v)
ef

## definindo um "grid" de prediÃ§Ã£o
gr <- expand.grid(seq(0,100, len=51), seq(0,100, len=51))
## visualizando a malha na qual serÃ¡ feita a prediÃ§Ã£o
points(dG)
points(gr, pch=19, cex=0.25, col=2)	

kc.ef <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=ef))

## visualizando os valores preditos na forma de um mapa (com diferentes opÃ§Ãµes de cores)
x11()
par(mfrow=c(2,2), mar=c(2,2,1,1))
image(kc.ef)
image(kc.ef, col=terrain.colors(21))
image(kc.ef, col=gray(seq(1,0, l=21)))
image(kc.ef, col=gray(seq(1,0, l=4)))
par(mfrow=c(1,1))

summary(kc.ef$predict)
image(kc.ef, col=terrain.colors(4), breaks=c(0, 30, 80, 100, 150))

image(kc.ef, val = sqrt(kc.ef$krige.var), coords.data=dG$coords)


##
## Fim da anÃ¡lise inicial/bÃ¡sica
##

##
## Podemos revisitar passos e usar outros mÃ©todos 
##
plot(v)		
ef <- eyefit(v)  ## ajustar e salvar 3 escolhas diferentes de modelos
ef

vf1 <- variofit(v, ini=ef[[1]])
vf1
vf2 <- variofit(v, ini=ef[[2]])
vf2
vf3 <- variofit(v, ini=ef[[3]])
vf3

par(mfrow=c(1,3))
plot(v)
lines(ef[[1]], col=2)
lines(vf1, col=4)

plot(v)
lines(ef[[2]], col=2)
lines(vf2, col=4)

plot(v)
lines(ef[[3]], col=2)
lines(vf3, col=4)
par(mfrow=c(1,1))


## estimando a media
ONES <- rep(1, nrow(dados))
D <- as.matrix(dist(dG$coords))
names(vf1)
vf1$cov.pars
SIGMA <- with(vf1, nugget * diag(nrow(dados)) + cov.pars[1] * exp(-D/cov.pars[2]))
dim(SIGMA)
iS.ONES <- solve(SIGMA, ONES)
(mu.gls <- solve(crossprod(iS.ONES, ONES), crossprod(iS.ONES, dG$data)))

## funÃ§Ã£o de (Log)verossimilhanÃ§a (escrita e forma ingenua)
ltheta <- function(pars, gd, printpars=FALSE){    
    ## pars = mu, sigma, phi, tau
    n <- nrow(gd$coords)
    D <- as.matrix(dist(gd$coords, upper=T, diag=T))
    SIGMA <- pars[4]^2 * diag(n) + pars[2]^2 * exp(-D/pars[3])
    logdetS <-  determinant(SIGMA, log=TRUE)$modulus
    res <- gd$data - pars[1]
    Q <- crossprod(res, solve(SIGMA, res))
    logVero <- -0.5*(n*log(2*pi) + logdetS + Q)
    if(printpars) print(c(pars, logVero))
    return(logVero)
}

lthetaConc <- function(pars, gd, printpars=FALSE){    
    ## pars : phi e nu = tau/sigma
    n <- nrow(gd$coords)
    D <- as.matrix(dist(gd$coords, upper=T, diag=T))
    SIGMA <- pars[1]^1 * diag(n) + exp(-D/pars[1])
    logdetS <-  determinant(SIGMA, log=TRUE)$modulus
    #mu <- 
        res <- gd$data - pars[1]
    Q <- crossprod(res, solve(SIGMA, res))
    logVero <- -0.5*(n*log(2*pi) + logdetS + Q)
    if(printpars) print(c(pars, logVero))
    return(logVero)
}

names(vf1)
(pars.vf1 <- with(vf1, 
                  c(mu.gls, sqrt(cov.pars[1]), cov.pars[2], sqrt(nugget))))

ltheta(pars.vf1, gd = dG)
ltheta(c(50, 15, 20, 0), gd = dG)

args(optim)

parsML <- optim(c(50, 15, 20, 0), ltheta, gd = dG, printpars=T, 
                method="L-BFGS-B", lower=c(-Inf, 0, 0, 0), 
                control=list(fnscale=-1))         
parsML$par^c(1,2,1,2)
mu.gls
vf1

## usando a funcao da geoR
ml <- likfit(dG, ini=c(300, 20))
ml

## definindo um "grid" de prediÃ§Ã£o
gr <- expand.grid(seq(0,100, len=51), seq(0,100, len=51))
points(dG)
points(gr, pch=19, cex=0.25, col=2)	


## fazendo a prediÃ§Ã£o espacial (krigagem)

## com opÃ§Ãµes de ajusta "a olho"(eyefit)
kc.ef1 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=ef[[1]]))
kc.ef2 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=ef[[2]]))
kc.ef3 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=ef[[3]]))

## com um variograma ajustado (variofit)
kc.vf1 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=vf1))
kc.vf2 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=vf2))
kc.vf3 <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=vf3))

## com parÃ¢metros ajustados por verossimilhanÃ§a (likfit)
kc.ml <- krige.conv(dG, loc=gr, krige=krige.control(obj.m=ml))
names(kc.ml)

## visualizando os valores preditos na forma de um mapa
par(mfrow=c(2,3), mar=c(2,2,1,1))
image(kc.ef1, col=terrain.colors(21))
image(kc.ef2, col=terrain.colors(21))
image(kc.ef3, col=terrain.colors(21))
image(kc.vf1, col=terrain.colors(21))
image(kc.vf2, col=terrain.colors(21))
image(kc.ml, col=terrain.colors(21))
par(mfrow=c(1,1), mar=c(2,2,1,1))

## escala comum
ZL <- range(
    range(kc.ef1$pred),
    range(kc.ef2$pred),
    range(kc.ef3$pred),
    range(kc.vf1$pred),
    range(kc.vf2$pred),
    range(kc.ml$pred)
)
par(mfrow=c(2,3), mar=c(2,2,1,1))
image(kc.ef1, col=terrain.colors(21), zlim=ZL)
image(kc.ef2, col=terrain.colors(21), zlim=ZL)
image(kc.ef3, col=terrain.colors(21), zlim=ZL)
image(kc.vf1, col=terrain.colors(21), zlim=ZL)
image(kc.vf2, col=terrain.colors(21), zlim=ZL)
image(kc.ml, col=terrain.colors(21), zlim=ZL)
par(mfrow=c(1,1))

## visualizando os erros padrÃ£o de prediÃ§Ã£o (mapa de incerteza de prediÃ§Ã£o)
ZL <- range(
    range(sqrt(kc.ef1$krige.var)),
    range(sqrt(kc.ef2$krige.var)),
    range(sqrt(kc.ef3$krige.var)),
    range(sqrt(kc.vf1$krige.var)),
    range(sqrt(kc.vf2$krige.var)),
    range(sqrt(kc.ml$krige.var))
)
par(mfrow=c(2,3), mar=c(2,2,1,1))
image(kc.ef1, val=sqrt(kc.ef1$krige.var), coords.data=dG$coords, zlim=ZL)
image(kc.ef2, val=sqrt(kc.ef2$krige.var), coords.data=dG$coords, zlim=ZL)
image(kc.ef3, val=sqrt(kc.ef3$krige.var), coords.data=dG$coords, zlim=ZL)
image(kc.vf1, val=sqrt(kc.vf1$krige.var), coords.data=dG$coords, zlim=ZL)
image(kc.vf2, val=sqrt(kc.vf2$krige.var), coords.data=dG$coords, zlim=ZL)
image(kc.ml, val=sqrt(kc.ml$krige.var), coords.data=dG$coords, zlim=ZL)
par(mfrow=c(1,1))

## voltando a variogramas
plot(v)
lines(ml, col="blue")
lines(vf1, col="red")
lines(vf2, col="green")
legend("topleft", legend=c("ml","vf1","vf2"),
       col=c("blue","red","green"), lty=1)

par(mfrow=c(1,3))
plot(kc.ml$pred, kc.vf1$pred, asp = 1); abline(0,1)
plot(kc.ml$pred, kc.vf2$pred, asp = 1); abline(0,1)
plot(kc.vf1$pred, kc.vf2$pred, asp = 1); abline(0,1)
par(mfrow=c(1,1))


## outras opÃ§Ãµes de visualizaÃ§Ã£o
par(mfrow=c(2,2))
contour(kc.ml)
image(kc.ml, col=terrain.colors(21))
points(dG, add=T)
contour(kc.ml, add=T, nlev=21)
persp(kc.ml)
persp(kc.ml, theta=20, phi=35)
par(mfrow=c(1,1))

##
## exportando prediÃ§Ãµes
##
## exportando somente um vetor dos atributos (formato texto)
write(kc.ml$pred, file="kc1.txt", ncol=1)
file.show("kc.ml.txt")    ## tecle c para terminar a exibiÃ§Ã£o do arquivo

## exportando somente os atributos como matriz (formato texto)
write(kc.ml$pred, ncol=51, file="kc2.txt")

## exportando atributos e locaÃ§Ãµes de prediÃ§Ã£o
names(kc.ml)
write.table(cbind(gr, data.frame(kc.ml[c(1,2)])),file="kc3.txt",
            row.names = FALSE)
## ou...
write.csv(cbind(gr, data.frame(kc.ml[c(1,2)])),file="kc4.csv",
          row.names = FALSE)
## ou...
write.csv2(cbind(gr, data.frame(kc.ml[c(1,2)])),file="kc5.csv",
           row.names = FALSE)

##
## simulacoes condicionais
##
args(krige.conv)
args(output.control)

kc.ml <- krige.conv(dG, loc=gr, 
                    krige=krige.control(obj.m=ml),
                    output = output.control(n.pred=1000))
names(kc.ml)
dim(kc.ml$simulations)

## a krigagem Ã© aproximada por simulaÃ§Ã£o condicional
predSim <- apply(kc.ml$simulations, 1, mean)
plot(kc.ml$pred, predSim); abline(0,1)
summary(kc.ml$pred - predSim)


## mapa de prob de estar acima de 65
limiar <- apply(kc.ml$simulations, 1, function(x) mean(x>65))
image(kc.ml, val=limiar)
contour(kc.ml, val=limiar)

## suponha agora que definimos como regiao de alto risco
## as quye possuiram prob > 0.7 de estar acima de 65
limiar.alto <- ifelse(limiar < 0.7, 0, 1)
image(kc.ml, val=limiar.alto, col=c(1,0))

## agora definindo 5 niveis de risco em fc da prob
image(kc.ml, val=limiar, breaks=seq(0, 1, by=0.2),
      col=c("blue","green","yellow","red","black"))

##
## definindo um polÃ­gono de interesse
points(dG)
bor <- locator(type="l")
polygon(bor, col=4)

## completar aqui a krigagem dentro do polygono


plot(dG)
v <- variog(dG, max.dist=80)
plot(v)

args(variog)
v0 <- variog(dG, max.dist=80, direction=0)
v45 <- variog(dG, max.dist=80, direction=pi/4)
v90 <- variog(dG, max.dist=80, direction=pi/2)
v135 <- variog(dG, max.dist=80, direction=3*pi/4)

lines(v0)
lines(v45, col=2)
lines(v90, col=3)
lines(v135,col=4)

v4 <- variog4(dG, max.dist=100)
plot(v4)
plot(v4, omn=T)


ml.iso <- likfit(dG, ini=c(500, 20))
ml.aniso <- likfit(dG, ini=c(500, 20), fix.psiR=F, fix.psiA=F)
ml.anisoR <- likfit(dG, ini=c(500, 20), fix.psiR=F, fix.psiA=T, psiA=0)

logLik(ml.iso)
logLik(ml.anisoR)
logLik(ml.aniso)


## Inf. Bayesiana
summary(dG)

args(krige.bayes)
args(model.control)
MC <- model.control()
args(prior.control)
PC <- prior.control(phi.discrete=seq(0, 40, len=41))
args(output.control)
OC <- output.control(n.post=1000)

kb <- krige.bayes(dG, model=MC, prior=PC, output=OC)
names(kb)
plot(kb)

## outra priori
PC <- prior.control(phi.discrete=seq(0, 40, len=41), 
                    phi.prior="rec")
kb <- krige.bayes(dG, model=MC, prior=PC, output=OC)
plot(kb)


## priori tb tem tau^2_Rel
PC <- prior.control(phi.discrete=seq(0, 50, len=26), 
                    phi.prior="rec",
                    tausq.rel.discrete=seq(0, 1.5, length=16),
                    tausq.rel.prior="rec")
kb <- krige.bayes(dG, model=MC, prior=PC, output=OC)
par(mfrow=c(1,2))
plot(kb)

##kb <- krige.bayes(dG, loc=gr, model=MC, prior=PC, output=OC)
##save(kb, file = "kb.RData")
download.file("http://www.leg.ufpr.br/~paulojus/dadosgeo/kb.RData", destfile="kb.RData")
load("kb.RData")

names(kb)
names(kb$post)
names(kb$pred)

image(kb)

p65bayes <- apply(kb$pred$simul, 1, function(x) mean(x>65))
image(kb, values = p65bayes)

## comparando as predicoes
image(kc.ml)
image(kb)

plot(kc.ml, kb$predictive$pred)




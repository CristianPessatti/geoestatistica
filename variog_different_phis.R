require(geoR)
require(purrr)

phis <- c(0,1,5,10)

dados_gen <- map(phis,
                 function(x){
                     grf(nx=20,ny=20, grid='reg',
                          mean=0,
                          nugget=0,
                          cov.pars=c(1,x))
                })
    
    
par(mfrow=c(2,2))
image(dados_gen[[1]], axes = FALSE, main = paste0("Phi = ", phis[1]))
image(dados_gen[[2]], axes = FALSE, main = paste0("Phi = ", phis[2]))
image(dados_gen[[3]], axes = FALSE, main = paste0("Phi = ", phis[3]))
image(dados_gen[[4]], axes = FALSE, main = paste0("Phi = ", phis[4]))

par(mfrow=c(2,2))
v1 <- variog(dados_gen[[1]])
plot(v1)

v2 <- variog(dados_gen[[2]])
plot(v2)

v3 <- variog(dados_gen[[3]])
plot(v3)

v4 <- variog(dados_gen[[4]])
plot(v4)

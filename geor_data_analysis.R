require(tidyverse)
require(geoR)

# PARANA ====================

data("parana")
parana
plot(parana)
?parana

points(parana)

v_parana <- variog(parana, max.dist = 500)

plot(v_parana)

vf_parana <- variofit(v_parana)

plot(v_parana)
lines(vf_parana, col=2)

gr_parana <- expand.grid(seq(0,800, len=100), seq(0,600, len=100))

kc_parana <- krige.conv(parana, loc=gr_parana, krige = krige.control(obj.model = vf_parana))

image(kc_parana)

# ca20 =========================
data("ca20")
ca20
plot(ca20)
?ca20

points(ca20)

v_ca20 <- variog(ca20)

vf_ca20 <- variofit(v_ca20)

plot(v_ca20)
lines(vf_ca20, col=2)

gr_ca20 <- expand.grid(seq(4800,6100, len=100), seq(4800,5800, len=100))

kc_ca20 <- krige.conv(ca20, loc=gr_ca20, krige = krige.control(obj.model = vf_ca20))

image(kc_ca20)

# KSAT ==============

data("Ksat")
Ksat
plot(Ksat)
?Ksat

points(Ksat)

v_Ksat <- variog(Ksat)

plot(v_Ksat)

v_Ksat <- variog(Ksat, lambda = 0)

plot(v_Ksat)

vf_Ksat <- variofit(v_Ksat)

plot(v_Ksat)
lines(vf_Ksat, col=2)

gr_Ksat <- expand.grid(seq(0,15, len=100), seq(0,22, len=100))

kc_Ksat <- krige.conv(Ksat, loc=gr_Ksat, krige = krige.control(obj.model = vf_Ksat))

image(kc_Ksat)

# WOLFCAMP =====================

data('wolfcamp')
wolfcamp
plot(wolfcamp, trend='1st')
?wolfcamp

points(wolfcamp)


par(mfrow = c(1,1))
v_wolfcamp <- variog(wolfcamp)

plot(v_wolfcamp)

v_wolfcamp1 <- variog(wolfcamp, trend = '1st')

plot(v_wolfcamp1)

v_wolfcamp2 <- variog(wolfcamp, trend = '2nd')

plot(v_wolfcamp2)

vf_wolfcamp <- variofit(v_wolfcamp)

plot(v_wolfcamp)
lines(vf_wolfcamp, col=2)

gr_wolfcamp <- expand.grid(seq(-200,200, len=100), seq(-200,200, len=100))

kc_wolfcamp <- krige.conv(wolfcamp, loc=gr_wolfcamp, krige = krige.control(obj.model = vf_wolfcamp))

image(kc_wolfcamp)


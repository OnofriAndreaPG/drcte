rm(list = ls())
data(chickweed)
head(chickweed)
data(chickweed0)
head(chickweed0)
moda <- drm(count ~ start + end, data = chickweed,
           fct = LL.3(), type = "event")
modb <- drm(count ~ time, data = chickweed0,
           fct = LL.3())
moda <- drmte(count ~ start + end, data = chickweed,
            fct = LL.3())
moda$fit$method
mod2 <- lm(count ~ 1, data = chickweed)
devtools::load_all()
newMod <- as.drc(mod2)
plotData(newMod, xlim = c(0, 300))

scaledH <- (newMod$"fit"$"hessian") / (2 * rse(object, TRUE))

deviance(moda)
stats:::deviance.default
weighted.residuals(modb)
stats:::deviance.lm
deviance(newMod)
deviance(mod)
vcov(mod2)
vcov(newMod)
summary(newMod)
df.residual(newMod)
newMod$df.residual

newMod$fit$value
newMod$type

newMod$fit$hessian
mod$fit$hessian
vcov(newMod)
newMod$type
mod$type
mod2$deviance
drc:::vcDisc(newMod)
vcov(mod)
summary(newMod)
(newMod$"fit"$"hessian") / (2 * drc:::rse(newMod, TRUE))

drc:::dev


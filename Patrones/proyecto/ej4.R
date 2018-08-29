azt <- read.table('azt.txt', header=TRUE)

azt.m1 <- glm(cbind(Y, N) ~ Race + AZT, family=binomial, data=azt)
azt.m2 <- glm(cbind(Y, N) ~ Race * AZT, family=binomial, data=azt)
azt.m3 <- glm(cbind(Y, N) ~ Race, family=binomial, data=azt)
azt.m4 <- glm(cbind(Y, N) ~ AZT, family=binomial, data=azt)

summary(azt.m1)
summary(azt.m2)
summary(azt.m3)
summary(azt.m4)
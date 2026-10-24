sqrt(mean((mydat$mhat-m_true)^2))
sd(mydat$mhat)/sd(m_true)
cor(mydat$mhat,m_true)
sd(m_true)

sqrt(mean((mydat$ehat-e0)^2))
sd(mydat$ehat)/sd(e0)
cor(mydat$ehat,e0)
sd(e0)

hist(e0)
hist(mydat$ehat)
hist((mydat$ehat-e0)^2)
hist((mydat$ehat-e0))
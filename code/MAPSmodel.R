library(ggplot2)
library(dplyr)

data <- read.table("./SynonymousPropSingleton_byTrimer_NEW.txt", header=TRUE)

withProp <- data %>%
        mutate(ProportionSingleton=Count_singleton/Count_variants)

# Plot proportion singleton again mu_snp

png("ProportionSingleton_by_mu.png")
ggplot(withProp, aes(x=mu_snp, y=ProportionSingleton)) +
        geom_point()
dev.off()

# Create linear model

linearMod <- lm(ProportionSingleton ~ mu_snp, data=withProp, weights = withProp$Count_variants)

summary(linearMod)$coefficients

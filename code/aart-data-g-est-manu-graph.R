################################################################################
# this file reproduces Figure~1
# it should be run after 'aart-data-analysis-36mo.R'
# and 'aart-data-analysis-60mo.R'
################################################################################

require(ggplot2)
require(extrafont)

#import data

#36 months
aart.est.36mo = read.csv("./results/aart-est-36mo.csv")
aart.est.36mo = aart.est.36mo[,-1]
aart.est.36mo = aart.est.36mo[1:33,]
aart.est.36mo$age = c(4:36)
aart.est.36mo$ci.high = aart.est.36mo$point.est + qnorm(0.975) * aart.est.36mo$std.error
aart.est.36mo$ci.low = aart.est.36mo$point.est - qnorm(0.975) * aart.est.36mo$std.error

#60 months
aart.est.60mo = read.csv("./results/aart-est-60mo.csv")
aart.est.60mo = aart.est.60mo[,-1]
aart.est.60mo = aart.est.60mo[1:59,]
aart.est.60mo$age = c(3:61)
aart.est.60mo$ci.high = aart.est.60mo$point.est + qnorm(0.975) * aart.est.60mo$std.error
aart.est.60mo$ci.low = aart.est.60mo$point.est - qnorm(0.975) * aart.est.60mo$std.error

#prepare plot df
plot_df.36mo = data.frame("loan.term" = rep("36-months", nrow(aart.est.36mo)),
                          "point.est" = aart.est.36mo$point.est,
                          "age" = aart.est.36mo$age,
                          "ci.low" = aart.est.36mo$ci.low,
                          "ci.high" = aart.est.36mo$ci.high)

plot_df.60mo = data.frame("loan.term" = rep("60-months", nrow(aart.est.60mo)),
                          "point.est" = aart.est.60mo$point.est,
                          "age" = aart.est.60mo$age,
                          "ci.low" = aart.est.60mo$ci.low,
                          "ci.high" = aart.est.60mo$ci.high)


plot_df = rbind(plot_df.36mo, plot_df.60mo)

ggplot() +
  geom_line(data=plot_df, aes(x=age, y=point.est), color="blue") +
  geom_point(data=plot_df, aes(x=age, y=point.est), color="blue") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=ci.low, ymax=ci.high),
              fill="lightblue", alpha=0.5) +
  xlab("Loan Age") + ylab("Estimated Probability Mass") +
  theme_bw() +
  theme(axis.title.x=element_text(size=10, family="Times New Roman"),
        axis.title.y=element_text(size=10,family="Times New Roman"),
        strip.text=element_text(size=10,family="Times New Roman"),
        axis.text=element_text(size=10,family="Times New Roman")) +
  facet_grid(cols = vars(loan.term))

ggsave("./results/aart-g-est.pdf",height=4,width=6,device = cairo_pdf)

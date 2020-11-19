# VICMATLAB MERRA-2 downscaling computation time estimation

72 in 17 hours

17*60/72

# average of 14.2 minutes per model day. 

# the time per day increases with the number of days, exponential
# it could takea  very long time to finish the processing

dat <- data.frame(x = c(3, 26, 72), y = c(7, 90, 17*60))
lf <- lm(log(y)~log(x), dat)

plot(y~x, dat)
yy <- exp(predict(lf))
lines(dat$x, yy, col = "blue")

yyy <- exp(predict(lf, newdata = data.frame(x = 365)))/60 # time to complete (hr)

# At this rate, it will take approximately 140 hours to complete the downscaling
# This is nearly 6 days.
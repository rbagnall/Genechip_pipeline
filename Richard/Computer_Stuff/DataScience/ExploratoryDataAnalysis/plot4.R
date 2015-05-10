#!/usr/bin/Rscript
# Rscript to draw four plots of data from household_power_consumption.txt file
# Full data set available for download from https://d396qusza40orc.cloudfront.net/exdata%2Fdata%2Fhousehold_power_consumption.zip

# required library is sqldf
require(sqldf)

# check to see if raw data exists in the current working directory and if not, exit with a comment
if(!file.exists("household_power_consumption.txt")){
    stop("The household_power_consumption.txt file must be present in your working directory", call. =FALSE)
}

#read in data from 1st to 2nd February, 2007 only
df <- read.csv.sql("household_power_consumption.txt", sql = "select * from file where (Date = '1/2/2007' or Date = '2/2/2007') ", header = TRUE, sep=";")

# convert date column to class:date
df$Date <- as.Date(df$Date, "%d/%m/%Y")

# convert time data to POSIX
df$Time <- strptime(paste(df$Date, df$Time), "%Y-%m-%d %H:%M:%S")

# start png device
png(filename= "plot4.png", bg="white", width=480, height=480, units="px")

# set up for four plots per image
par(mfrow = c(2, 2))


# top left plot
plot(df$Time, df$Global_active_power, type="l", ylab="Global Active Power", xlab="")

# top right plot
plot(df$Time, df$Voltage, type="l", ylab="Voltage", xlab="datetime", yaxt="n")

# axis labels
axis(2, at=c(234, 236, 238, 240, 242, 244, 246), labels=c("234", "", "238", "", "242", "", "246"), las=0)

#bottom left plot
plot(df$Time, df$Sub_metering_1, type="n", ylab="Energy sub metering", xlab="")

# plot 3 sets of data, using different colours for each
lines(df$Time, df$Sub_metering_1, type="l")
lines(df$Time, df$Sub_metering_2, type="l", col="red")
lines(df$Time, df$Sub_metering_3, type="l", col="blue")

# add legend, without a border (bty="n")
legend("topright", c("Sub_metering_1", "Sub_metering_2", "Sub_metering_3"), col=c("black", "red", "blue"), lty=1, bty="n")

#bottom right plot
plot(df$Time, df$Global_reactive_power, type="l", xlab="datetime", ylab="Global_reactive_power")

dev.off()


#!/usr/bin/Rscript
# Rscript to draw a plot of Global_active_power and Time data from household_power_consumption.txt file
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

# convert time column to POSIX
df$Time <- strptime(paste(df$Date, df$Time), "%Y-%m-%d %H:%M:%S")

# start png device
png(filename= "plot2.png", bg="white", width=480, height=480, units="px")

# plot Global_active_power and Time data
plot(df$Time, df$Global_active_power, type="l", ylab="Global Active Power (kilowatts)", xlab="")

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Author: Zhaodong Hao
# Create data: 2018-05-30
# Description: Plotting chromosome idiogram
# Attention: you can use Inkscape to open the output file (.svg) and save as a PDF file
# Contact: haozd1992@163.com
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# set the working dictionary as the path of 'source'd file
whereFrom=sys.calls()[[1]]
# This should be an expression that looks something like
# source("pathname/myfilename.R")
whereFrom=as.character(whereFrom[2]) # get the pathname/filename
whereFrom=paste(getwd(),whereFrom,sep="/") # prefix it with the current working directory
pathnameIndex=gregexpr(".*/",whereFrom) # we want the string up to the final '/'
pathnameLength=attr(pathnameIndex[[1]],"match.length")
whereFrom=substr(whereFrom,1,pathnameLength-1)
setwd(whereFrom)

# read the karyotype files (核型文件五列，第一列染色体，第二、三列染色体长度，第四、五列着丝粒位置)
karyotype <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F, colClasses = c("character", "integer", "integer", "integer", "integer"))

# A4 paper (210mm * 297mm)
# left, right, top, bottom ---margins--- 20mm, 20mm, 25mm, 25mm (plot region: 170mm * 247mm)
# 1mm = 3.543307px


### rect model: <rect x y rx ry width height style="fill:#fdcdac;stroke:black;stroke-width:2;opacity:0.5"/>
# x,y 矩形左上角点的坐标； rx,ry 圆角矩形；width,height 矩形的宽和高（从左上角的点向右向下延伸）； fill 填充颜色；stroke 线的颜色；stroke-width 线的宽度；opacity 透明度

# karyotype rect

# 170mm, n个染色体，每个染色体+空隙宽度为(170/n)mm，染色体与空隙等宽，染色体宽度[(170/n)/2]mm，换算成px

chr_width <- (floor(170 / nrow(karyotype))/2) * 3.543307  # 170可以根据染色体条数适当调整，这里24条染色体设定的画布宽度为170mm，可以多试试几次，选择最合适的画布宽度

#

for (i in 1:nrow(karyotype)) {
  karyotype[i, 6] <- 20 * 3.543307 + (i-1) * (chr_width * 2)
}
names(karyotype)[6] <- "x1"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 7] <- ((150 + 25) - (150 * karyotype[i,3]) / max(karyotype$End)) * 3.543307 + chr_width/2
}
names(karyotype)[7] <- "y1"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 8] <- karyotype[i, 6] + chr_width
}
names(karyotype)[8] <- "x2"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 9] <- karyotype[i, 7]
}
names(karyotype)[9] <- "y2"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 10] <- karyotype[i, 8]
}
names(karyotype)[10] <- "x3"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 11] <- karyotype[i, 9] - chr_width/2 + ((150 * karyotype[i,4]) / max(karyotype$End)) * 3.543307
}
names(karyotype)[11] <- "y3"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 12] <- karyotype[i, 6]
}
names(karyotype)[12] <- "x4"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 13] <- karyotype[i, 9] - chr_width/2 + ((150 * karyotype[i,5]) / max(karyotype$End)) * 3.543307
}
names(karyotype)[13] <- "y4"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 14] <- karyotype[i, 6]
}
names(karyotype)[14] <- "x5"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 15] <- karyotype[i, 9] - chr_width/2 + ((150 * karyotype[i,3]) / max(karyotype$End)) * 3.543307 - chr_width/2
}
names(karyotype)[15] <- "y5"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 16] <- karyotype[i, 8]
}
names(karyotype)[16] <- "x6"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 17] <- karyotype[i, 15]
}
names(karyotype)[17] <- "y6"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 18] <- karyotype[i, 8]
}
names(karyotype)[18] <- "x7"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 19] <- karyotype[i, 13]
}
names(karyotype)[19] <- "y7"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 20] <- karyotype[i, 6]
}
names(karyotype)[20] <- "x8"

for (i in 1:nrow(karyotype)) {
  karyotype[i, 21] <- karyotype[i, 11]
}
names(karyotype)[21] <- "y8"


karyotype$path = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1, " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2, " L", karyotype$x3, ",", karyotype$y3, " L", karyotype$x4, ",", karyotype$y4, " L", karyotype$x5, ",", karyotype$y5, " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x6, ",", karyotype$y6, " L", karyotype$x7, ",", karyotype$y7, " L", karyotype$x8, ",", karyotype$y8, " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")
karyotype$text = paste("<text x=\"", (karyotype$x1 + karyotype$x6)/2 - nchar(karyotype$Chr) * 2.1, "\" y=\"", (150 + 25) * 3.543307 + 15, "\" font-size=\"9\" fill=\"black\" >", karyotype$Chr, "</text>", sep = "")


# read the data_1 files (染色体上的barplot文件，四列，第一列染色体，第二、三列染色体上的位置，第四列对应的数值)
data <- read.table("data_1.txt", sep = "\t", header = T, stringsAsFactors = F, colClasses = c("character", "integer", "integer", "integer"))
require(graphics)

# 自带渐变色
# color <- as.data.frame(topo.colors(max(data$Value), alpha = 1))
# color <- as.data.frame(heat.colors(max(data$Value), alpha = 1))
# color <- as.data.frame(cm.colors(max(data$Value), alpha = 1))
# color <- as.data.frame(terrain.colors(max(data$Value), alpha = 1))
# color <- as.data.frame(rainbow(max(data$Value), s = 1, v = 1, start = 0, end = max(1, max(data$Value) - 1)/max(data$Value), alpha = 1))

# 手动设置渐变色（其中"#8dd3c7", "#ffffb3", "#fb8072"可以根据需要进行修改，注意同时修改下面的data1_legend）
color <- as.data.frame(colorRampPalette(c("#8dd3c7", "#ffffb3", "#fb8072"))(max(data$Value)))

for (i in 1:nrow(data)){
  data[i,5] <- substring(color[data[i,4], 1], 1, 7)
}
names(data)[5] <- "color"



for (i in 1:nrow(data)){
  if (data[i,3] <= max(karyotype$End) * chr_width/2/3.543307 / 150){
    data[i, 6] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] < max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] > max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] < karyotype[data[i,1],4]) {
    data[i, 6] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (max(karyotype$End) * chr_width/2/3.543307 / 150 <= data[i,2] & data[i,3] <= karyotype[data[i,1],4]){
    data[i, 6] <- karyotype[data[i,1],6]
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] > karyotype[data[i,1],4] & data[i,3] <= karyotype[data[i,1],5]) {
    data[i, 6] <- karyotype[data[i,1],6]
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],4] <= data[i,2] & data[i,3] <= karyotype[data[i,1],5]){
    data[i, 6] <- karyotype[data[i,1],6] + ((data[i,2] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],4] & data[i,2] < karyotype[data[i,1],5] & data[i,3] > karyotype[data[i,1],5]) {
    data[i, 6] <- karyotype[data[i,1],6] + ((data[i,2] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],5] <= data[i,2] & data[i,3] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 6] <- karyotype[data[i,1],6]
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) & data[i,3] > karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) ){
    data[i, 6] <- karyotype[data[i,1],6]
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 6] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 7] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
}

names(data)[5] <- "x1"
names(data)[6] <- "y1"

for (i in 1:nrow(data)){
  if (data[i,3] <= max(karyotype$End) * chr_width/2/3.543307 / 150){
    data[i, 8] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] < max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] > max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] < karyotype[data[i,1],4]) {
    data[i, 8] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (max(karyotype$End) * chr_width/2/3.543307 / 150 <= data[i,2] & data[i,3] <= karyotype[data[i,1],4]){
    data[i, 8] <- karyotype[data[i,1],8]
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] > karyotype[data[i,1],4] & data[i,3] <= karyotype[data[i,1],5]) {
    data[i, 8] <- karyotype[data[i,1],8]
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],4] <= data[i,2] & data[i,3] <= karyotype[data[i,1],5]){
    data[i, 8] <- karyotype[data[i,1],8] - ((data[i,2] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],4] & data[i,2] < karyotype[data[i,1],5] & data[i,3] > karyotype[data[i,1],5]) {
    data[i, 8] <- karyotype[data[i,1],8] - ((data[i,2] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],5] <= data[i,2] & data[i,3] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 8] <- karyotype[data[i,1],8]
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) & data[i,3] > karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) ){
    data[i, 8] <- karyotype[data[i,1],8]
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 8] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,2]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 9] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,2] * 150 / max(karyotype$End) * 3.543307
  }
}

names(data)[7] <- "x2"
names(data)[8] <- "y2"

for (i in 1:nrow(data)){
  if (data[i,3] <= max(karyotype$End) * chr_width/2/3.543307 / 150){
    data[i, 10] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] < max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] > max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] < karyotype[data[i,1],4]) {
    data[i, 10] <- karyotype[data[i,1],8]
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (max(karyotype$End) * chr_width/2/3.543307 / 150 <= data[i,2] & data[i,3] <= karyotype[data[i,1],4]){
    data[i, 10] <- karyotype[data[i,1],8]
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] > karyotype[data[i,1],4] & data[i,3] <= karyotype[data[i,1],5]) {
    data[i, 10] <- karyotype[data[i,1],8] - ((data[i,3] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],4] <= data[i,2] & data[i,3] <= karyotype[data[i,1],5]){
    data[i, 10] <- karyotype[data[i,1],8] - ((data[i,3] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],4] & data[i,2] < karyotype[data[i,1],5] & data[i,3] > karyotype[data[i,1],5]) {
    data[i, 10] <- karyotype[data[i,1],6]
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],5] <= data[i,2] & data[i,3] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 10] <- karyotype[data[i,1],8]
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) & data[i,3] > karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) ){
    data[i, 10] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 10] <- karyotype[data[i,1],8] - (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 11] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
}

names(data)[9] <- "x3"
names(data)[10] <- "y3"

for (i in 1:nrow(data)){
  if (data[i,3] <= max(karyotype$End) * chr_width/2/3.543307 / 150){
    data[i, 12] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - ((max(karyotype$End) * chr_width/2/3.543307 / 150 - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] < max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] > max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] < karyotype[data[i,1],4]) {
    data[i, 12] <- karyotype[data[i,1],6]
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (max(karyotype$End) * chr_width/2/3.543307 / 150 <= data[i,2] & data[i,3] <= karyotype[data[i,1],4]){
    data[i, 12] <- karyotype[data[i,1],6]
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] > karyotype[data[i,1],4] & data[i,3] <= karyotype[data[i,1],5]) {
    data[i, 12] <- karyotype[data[i,1],6] + ((data[i,3] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],4] <= data[i,2] & data[i,3] <= karyotype[data[i,1],5]){
    data[i, 12] <- karyotype[data[i,1],6] + ((data[i,3] - karyotype[data[i,1],4]) * (chr_width/3.543307* max(karyotype$End) / 150) / (karyotype[data[i,1],5] - karyotype[data[i,1],4])) * 150 / max(karyotype$End) * 3.543307
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],4] & data[i,2] < karyotype[data[i,1],5] & data[i,3] > karyotype[data[i,1],5]) {
    data[i, 12] <- karyotype[data[i,1],8]
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (karyotype[data[i,1],5] <= data[i,2] & data[i,3] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 12] <- karyotype[data[i,1],6]
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) & data[i,3] > karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) ){
    data[i, 12] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
  else if (data[i,2] >= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)){
    data[i, 12] <- karyotype[data[i,1],6] + (chr_width/2 - ((chr_width/2)^2 - (chr_width/2 - (karyotype[data[i,1],3] - data[i,3]) * 150 / max(karyotype$End) * 3.543307)^2)^0.5)
    data[i, 13] <- karyotype[data[i,1], 7] - chr_width/2 + data[i,3] * 150 / max(karyotype$End) * 3.543307
  }
}

names(data)[11] <- "x4"
names(data)[12] <- "y4"

for (i in 1:nrow(data)){
  if (data[i,3] <= max(karyotype$End) * chr_width/2/3.543307 / 150) {
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,6], ",", data[i,7], "\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] < max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] > max(karyotype$End) * chr_width/2/3.543307 / 150 & data[i,3] < karyotype[data[i,1],4]) {
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", karyotype[data[i,1],8], ",", karyotype[data[i,1],9], " L", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13], " L", karyotype[data[i,1],6], ",", karyotype[data[i,1],7], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,6], ",", data[i,7], "\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] > karyotype[data[i,1],4] & data[i,3] <= karyotype[data[i,1],5]){
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " L", karyotype[data[i,1],10], ",", karyotype[data[i,1],11], " L", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13]," L", karyotype[data[i,1],20], ",", karyotype[data[i,1],21], " Z" ,"\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] >= karyotype[data[i,1],4] & data[i,2] < karyotype[data[i,1],5] & data[i,3] > karyotype[data[i,1],5]){
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " L", karyotype[data[i,1],12], ",", karyotype[data[i,1],13], " L", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13]," L", karyotype[data[i,1],18], ",", karyotype[data[i,1],19], " Z" ,"\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] <= karyotype[data[i,1],4] & data[i,3] >= karyotype[data[i,1],5]){
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " L", karyotype[data[i,1],10], ",", karyotype[data[i,1],11], " L", karyotype[data[i,1],12], ",", karyotype[data[i,1],13]," L", data[i,12], ",", data[i,13], " L", data[i,10], ",", data[i,11]," L", karyotype[data[i,1],18], ",", karyotype[data[i,1],19], " L", karyotype[data[i,1],20], ",", karyotype[data[i,1],21], " Z" ,"\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] <= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150) & data[i,3] > karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)) {
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " L", karyotype[data[i,1],16], ",", karyotype[data[i,1],17], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", karyotype[data[i,1],14], ",", karyotype[data[i,1],15], " L", data[i,6], ",", data[i,7], "\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else if (data[i,2] >= karyotype[data[i,1],3] - (max(karyotype$End) * chr_width/2/3.543307 / 150)) {
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13], " A", chr_width/2, ",", chr_width/2, " 0 0,1 ", data[i,6], ",", data[i,7], "\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
  else{
    data[i,14] = paste("<path d=\"M", data[i,6], ",", data[i,7], " L", data[i,8], ",", data[i,9], " L", data[i,10], ",", data[i,11], " L", data[i,12], ",", data[i,13], " Z" ,"\" style=\"fill:", data[i, 5], "; stroke:", data[i, 5], "; stroke-width:0.5\"/>", sep = "")
  }
}

names(data)[14] <- "rect"

# legend for data1

data1_legend <- as.data.frame(1:1000)
names(data1_legend)[1] <- "Num"

# 颜色与上面data1的设置保持一致
color_rect <- as.data.frame(colorRampPalette(c("#8dd3c7", "#ffffb3", "#fb8072"))(max(data1_legend[,1])))

for (i in 1:nrow(data1_legend)){
  data1_legend[i,2] <- substring(color_rect[data1_legend[i,1], 1], 1, 7)
}
names(data1_legend)[2] <- "color"

for (i in 1:nrow(data1_legend)){
  data1_legend[i,3] <- (20 + 140) * 3.543307 + (i - 1) * ((20 * 3.543307) / 1000)
}
names(data1_legend)[3] <- "x1"

for (i in 1:nrow(data1_legend)){
  data1_legend[i,4] <- (25 + 10) * 3.543307
}
names(data1_legend)[4] <- "y1"

for (i in 1:nrow(data1_legend)){
  data1_legend[i,5] <- paste("<rect x=\"", data1_legend[i,3], "\" y=\"", data1_legend[i,4], "\" width=\"", (20 * 3.543307) / 1000,  "\" height=\"", 4 * 3.543307, "\" style=\"fill:", data1_legend[i,2], ";stroke:none\"/>", sep = "")
}
names(data1_legend)[5] <- "data1_legend"

legend1_text <- data.frame(paste("<text x=\"", data1_legend[1,3], "\" y=\"", data1_legend[1,4] + 8 * 3.543307 - 3, "\" font-size=\"12\" fill=\"black\" >Low</text>", sep = ""),
                           paste("<text x=\"", data1_legend[1000,3] - 25, "\" y=\"", data1_legend[1,4] + 8 * 3.543307 - 3, "\" font-size=\"12\" fill=\"black\" >High</text>", sep = "")
)


# read the data_2 files (染色体旁边的小图数据文件，六列，第一列名称，第二列形状（只支持box正方形、triangle三角形和circle圆），第四列染色体，第五、六列染色体上的位置，第七类颜色，注意同一类的type对应的颜色一致)
data_interval <- read.table("data_2.txt", sep = "\t", header = T, stringsAsFactors = F, colClasses = c("character", "character", "integer", "integer", "integer", "character"))

for (i in 1:nrow(data_interval)){
  data_interval[i, 7] <- karyotype[data_interval[i,3], 8] + chr_width / 3
  data_interval[i, 8] <- karyotype[data_interval[i,3], 9] - chr_width/2 + (data_interval[i,4] + data_interval[i,5]) / 2 * 150 / max(karyotype$End) * 3.543307
}
names(data_interval)[7] <- "x"
names(data_interval)[8] <- "y"

for (i in 1:nrow(data_interval)){
  if (data_interval[i, 2] == "triangle") {
    data_interval[i,9] = paste("<path d=\"M", data_interval[i,7] - chr_width / 6, ",", data_interval[i,8] + chr_width / 6, " L", data_interval[i,7] + chr_width / 6, ",", data_interval[i,8] + chr_width / 6, " L", data_interval[i,7], ",", data_interval[i,8] - chr_width / 6, " Z" ,"\" style=\"fill:#", data_interval[i, 6], ";stroke:none; opacity:0.6\"/>", sep = "")
  }
  else if (data_interval[i, 2] == "box") {
    data_interval[i,9] <- paste("<rect x=\"", data_interval[i,7] - chr_width / 6, "\" y=\"", data_interval[i,8] - chr_width / 6, "\" width=\"", chr_width / 3,  "\" height=\"", chr_width / 3, "\" style=\"fill:#", data_interval[i, 6], "; stroke:none; opacity:0.6\"/>", sep = "")
  }
  else if (data_interval[i, 2] == "circle") {
    data_interval[i,9] <- paste("<circle cx=\"", data_interval[i,7], "\" cy=\"", data_interval[i,8], "\" r=\"", chr_width / 6, "\" style=\"fill:#", data_interval[i, 6], "; stroke:none; opacity:0.6\"/>", sep = "")
  }
}
names(data_interval)[9] <- "shape"

# legend for data2
data2_legend <- data_interval[!duplicated(data_interval$Type), 1:6]
data2_legend <- data2_legend[order(data2_legend$Shape, data2_legend$color),]

for (i in 1:nrow(data2_legend)){
  data2_legend[i, 7] <- data1_legend[1, 3]
}
names(data2_legend)[7] <- "x"

for (i in 1:nrow(data2_legend)){
  data2_legend[i, 8] <- data1_legend[1, 4] + 4 + (12 + (i - 1) * 4 ) * 3.543307
}
names(data2_legend)[8] <- "y"

for (i in 1:nrow(data2_legend)){
  data2_legend[i, 9] <- data2_legend[i, 7] + 4
}
names(data2_legend)[9] <- "x1"

for (i in 1:nrow(data2_legend)){
  data2_legend[i, 10] <- data2_legend[i, 8] - 4 * 3.543307 / 2
}
names(data2_legend)[10] <- "y1"

for (i in 1:nrow(data2_legend)){
  if (data2_legend[i, 2] == "triangle") {
    data2_legend[i,11] = paste("<path d=\"M", data2_legend[i,9] - 4, ",", data2_legend[i,10] + 4, " L", data2_legend[i,9] + 4, ",", data2_legend[i,10] + 4, " L", data2_legend[i,9], ",", data2_legend[i,10] - 4, " Z" ,"\" style=\"fill:#", data2_legend[i, 6], "; stroke:none\"/>", sep = "")
  }
  else if (data2_legend[i, 2] == "box") {
    data2_legend[i,11] <- paste("<rect x=\"", data2_legend[i,9] - 4, "\" y=\"", data2_legend[i,10] - 4, "\" width=\"", 8,  "\" height=\"", 8, "\" style=\"fill:#", data2_legend[i, 6], ";stroke:none\"/>", sep = "")
  }
  else if (data2_legend[i, 2] == "circle") {
    data2_legend[i,11] <- paste("<circle cx=\"", data2_legend[i,9], "\" cy=\"", data2_legend[i,10], "\" r=\"", 4, "\" style=\"fill:#", data2_legend[i, 6], ";stroke:none\"/>", sep = "")
  }
}
names(data2_legend)[11] <- "shape"

for (i in 1:nrow(data2_legend)){
  data2_legend[i,12] <- paste("<text x=\"", data2_legend[i,7] + 15, "\" y=\"", data2_legend[i,8] - (4 * 3.543307 / 2 - 4), "\" font-size=\"12\" fill=\"black\" >", data2_legend[i,1],"</text>", sep = "")
}
names(data2_legend)[12] <- "name"


#write .svg file

first_line <- data.frame("<?xml version=\"1.0\" standalone=\"no\"?>",
                         "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                         "",
                         paste("<svg id=\"svg\" width=\"744.0945\" height=\"1052.362\">", "\t")
)

write.table(first_line[1, 1], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(first_line[1, 2], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 3], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 4], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 5], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

write.table(data$rect, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(karyotype$path, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(karyotype$text, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(data1_legend$data1_legend, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(legend1_text[1,1], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(legend1_text[1,2], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(data_interval$shape, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(data2_legend$shape, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(data2_legend$name, "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

last_line <- data.frame(paste("</svg>"))
write.table(last_line[1, 1], "chromsome.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

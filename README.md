# Idiogram-plotting

Download the R script and three input files in the same directory. Then, run the R script and you will get a .svg file under the same directory. You can use Inkscape to open the output file (.svg) and save as a PDF file.

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/Idiogram-plotting.jpg)

There are three input files you need to prepare if you want to plot your own figure:

karyotype.txt

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/karyotype.jpg)

data_1.txt

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/data_1.jpg)

data_2.txt

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/data_2.jpg)

Note, the Chr in data_1 and data_2 need to be filled with numbers. In this case, "X" in the karyotype file is changed to be 23 and "Y" in the karyotype file is changed to be 24 in both data file.

If the number of chromosome of you species is less than 24, you can changed the number "170" in the line 30.

10 chromosomes with the value "170" not changed

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/10_170.jpg)

10 chromosomes with the value "170" changed to be "150"

![image](https://github.com/TickingClock1992/Idiogram-plotting/blob/master/image/10_150.jpg)

So if you really want to do better, you had better try more than once.

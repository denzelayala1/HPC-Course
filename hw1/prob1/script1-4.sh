#!/bin/bash

grep " * residual (" ../HPC1_Model2.log > tmpfile.txt
awk '{print $5}' tmpfile.txt > prob1-3.out

countMax=$(wc -l tmpfile.txt | awk '{print $1}' )
echo "Total Number of lines is $countMax" | tee prob1-4.txt

awk '
BEGIN{countf=0.0; totSum=0.0;}

{
totSum+= $(NF-1);
countf+=$5
}

END{printf "The total residual count = %.3e\nThe total runtime = %.2f\n",countf, totSum}

' tmpfile.txt | tee -a prob1-4.txt

rm tmpfile.txt

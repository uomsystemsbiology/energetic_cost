#! /bin/sh

# Inform user
echo Creating $1_icd.m

#Create the component ICD file complete with headers.
echo "function icd = $1_icd" > $1_icd.m
echo "%% Component icd file ($1_icd.m)" >> $1_icd.m
echo "%% Generated by MTT at `date`" >> $1_icd.m

#Write out the variables 
gawk '{
       if (NF==2) {i++; print "icd."$1 "\t = \""$2"\";"}
     
     }
     END{
       if (i==0) print "icd = 0;"
        }' $1_icd.txt >> $1_icd.m
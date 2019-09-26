#!/usr/bin/env bash

wget https://github.com/mothur/mothur/releases/download/v.1.42.3/Mothur.linux_64.zip
mv Mothur.linux_64.zip code/
#Unzip into code/ directory as UoW may have bug in moving directories
unzip code/Mothur.linux_64.zip -d code/
#tidy up code/ folder
rm -R code/Mothur.linux_64.zip
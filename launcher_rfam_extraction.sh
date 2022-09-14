#!/bin/bash

for file in ./*.fa 
do 
	echo $file && cmscan --rfam --cut_ga --nohmmonly --tblout ${file/.fa/.rfam.tblout} --fmt 2 --clanin ../../disk1/Rfam/Rfam.14.1.clanin ../../disk1/Rfam/Rfam.cm $file > ${file/.fa/.cmscan}   
done


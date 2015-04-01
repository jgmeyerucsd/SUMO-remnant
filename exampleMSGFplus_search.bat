:: change directory
F:
cd \msgfplus.20140716\
:: the line below just gives the options, press any key to start the search
java -Xmx12000M -d64 -jar MSGFplus.jar
pause
:: each option in the lines below must be set according to the options printed from the above command depending on the specific type of data
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\exactive\2015\mg132\sbl_009779.mzXML -o output\exactive\2015\Mg132\mg132pl_4_light1.mzid -e 1 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod digly.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\exactive\2015\mg132\sbl_009779.mzXML -o output\exactive\2015\Mg132\mg132pl_4_heavy1.mzid -e 1 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod diglyh.txt
pause
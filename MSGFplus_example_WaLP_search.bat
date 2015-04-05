F:
cd \msgfplus.20140716\
:: the first 4 lines run searches with fixed heavy K, the last 4 lines run searches with fixed light K
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02118_Q.mzXML -o output\fusion\2014\december\qh3.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod diglyh.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02119_R.mzXML -o output\fusion\2014\december\rh3.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod diglyh.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02120_S.mzXML -o output\fusion\2014\december\sh3.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod diglyh.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02121_T.mzXML -o output\fusion\2014\december\th3.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod diglyh.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02118_Q.mzXML -o output\fusion\2014\december\q4.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod digly.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02119_R.mzXML -o output\fusion\2014\december\r4.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod digly.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02120_S.mzXML -o output\fusion\2014\december\s4.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod digly.txt
java -Xmx12000M -d64 -jar MSGFplus.jar -thread 8 -s C:\Inetpub\wwwroot\ISB\data\fusion\141205_MS2-CID-ETD_Meyer-J_UCSD\t02121_T.mzXML -o output\fusion\2014\december\t4.mzid -e 10 -ti 0,2 -tda 1 -inst 1 -ntt 0 -m 0 -t 10ppm -d database\110712_human.cc.fasta -mod digly.txt
pause

python timeRec.py "./ptre qrs.ptre 205L.csv 3 0" >logtecg1
python timeRec.py "./ptre qrs.ptre 221L.csv 3 0" >logtecg2
python timeRec.py "./ptre qrs.ptre 123L.csv 3 0" >logtecg3
python timeRec.py "./ptre ecgstl.ptre bflist205.txt 2 2" >logtspt1
python timeRec.py "./ptre ecgstl.ptre bflist221.txt 2 2" >logtspt2
python timeRec.py "./ptre ecgstl.ptre bflist123.txt 2 2" >logtspt3
python timeRec.py "./ptre param4ecg.ptre debugTest.csv 4 0 ecg2.label" >logtpi1
python timeRec.py "./ptre param5ecg.ptre debugTest.csv 5 0 ecg2.label" >logtpi2
python timeRec.py "./ptre param6ecg.ptre debugTest.csv 6 0 ecg2.label" >logtpi3

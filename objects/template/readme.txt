
To implement new analysis, do

1)atom-mkanalysis --skip-hep-data 1603.09222
2)create/edit info.txt
3)edit .info file (add keyword and luminosity)

To validate the analysis, do

1)cmake .
2)source setup.sh
3)generate events
4)make
5)atom -b script
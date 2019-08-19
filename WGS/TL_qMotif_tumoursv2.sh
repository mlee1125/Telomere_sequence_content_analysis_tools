#!/bin/sh

#for f in `ls *Tumour*qmotif_4TTAGGGu.ini.xml`
echo "id,includesN,unmappedN,includesT,unmappedT"

#for f in `ls *PrimaryTumour*qmotif_4TTAGGGu.ini.xml`
#for f in `ls *Tumour*qmotif_4TTAGGGu.ini.xml`
for f in `ls *tumour*qmotif_9TTAGGGc.qmotif.xml`
do
	id=`echo $f | sed -r "s/^[^_]*_([A-Z]*)([0-9]*)_.*$/\1_\2/"`
	id2=`echo $f | sed -r "s/^[^_]*_([A-Z]*)([0-9]*)_.*$/\1\2/"`
	includesT=`cat $f | grep "scaledIncludes" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`
	unmappedT=`cat $f | grep "scaledUnmapped" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`
	
	n=`ls *$id2*normal*qmotif_9TTAGGGc.qmotif.xml`
	includesN=`cat $n | grep "scaledIncludes" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`
	unmappedN=`cat $n | grep "scaledUnmapped" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`

	echo "$id,$includesN,$unmappedN,$includesT,$unmappedT"
done


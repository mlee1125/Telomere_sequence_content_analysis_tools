#!/bin/sh

echo "id,includes,unmapped"

for f in `ls *qmotif_9TTAGGGu.qmotif.xml`
do
	id=`echo $f`
	includes=`cat $f | grep "scaledIncludes" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`
	unmapped=`cat $f | grep "scaledUnmapped" | sed -r "s/.*count..([0-9]*)[^0-9]*/\1/"`
	
	echo "$id,$includes,$unmapped"
done


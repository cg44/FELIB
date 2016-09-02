#!/bin/tcsh
# To run and check all Level 0 Programs
#
set PROGRAMS="1p1 1p2 2p1 2p2 3p1 3p2 4p1 4p2 5p1 5p2"
set DIFFIL = "#differences"

/bin/rm -f $DIFFIL

switch ($1)
case clean:
		/bin/rm -f \#*
		breaksw
default:
	foreach i ($PROGRAMS)
        echo "*********************************************"
		echo "Running " seg$i
		set PROG=seg$i
		set DATA=dat$i
		set RES1=\#res$i
		set RES2=res$i
		echo "Data="$DATA" Results="$RES1
		$PROG <../data/$DATA >$RES1
        sed -e "/E[-+]/s/E/D/g" <  $RES1 > $RES1.tmp
		echo "Differencing results (old/new)"
		diff ../results/$RES2 $RES1.tmp >\#diff$i
        if(`ls -last \#diff$i | awk '{print $6}'` != 0) then
			echo "*** Differences found: see file #diff$i"
			echo "*** Differences found in "$PROG": see file #diff$i" >>$DIFFIL
		else
			echo "No differences found"
			/bin/rm -f \#diff$i
		endif
	endif
	end
	echo "*********************************************"
    echo "Contents of difference file"
	cat $DIFFIL
	echo "*********************************************"
endsw

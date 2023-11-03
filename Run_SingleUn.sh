( python Run.py --P WWW_Dim6 --UN $1 > WWW_Dim6_$1.debug 2>&1 ; echo "WWW_Dim6 $1 done" ) &  
( python Run.py --P WWZ_Dim6 --UN $1 > WWZ_Dim6_$1.debug 2>&1 ; echo "WWZ_Dim6 $1 done" ) & 
( python Run.py --P WZZ_Dim6 --UN $1 > WZZ_Dim6_$1.debug 2>&1 ; echo "WZZ_Dim6 $1 done" ) & 
( python Run.py --P ZZZ_Dim6 --UN $1 > ZZZ_Dim6_$1.debug 2>&1 ; echo "ZZZ_Dim6 $1 done" ) & 

( python Run.py --P WWW_Dim8 --UN $1 > WWW_Dim8_$1.debug 2>&1 ; echo "WWW_Dim8 $1 done" ) & 
( python Run.py --P WWZ_Dim8 --UN $1 > WWZ_Dim8_$1.debug 2>&1 ; echo "WWZ_Dim8 $1 done" ) & 
( python Run.py --P WZZ_Dim8 --UN $1 > WZZ_Dim8_$1.debug 2>&1 ; echo "WZZ_Dim8 $1 done" ) & 
( python Run.py --P ZZZ_Dim8 --UN $1 > ZZZ_Dim8_$1.debug 2>&1 ; echo "ZZZ_Dim8 $1 done" ) & 
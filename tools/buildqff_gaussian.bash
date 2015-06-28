#!/bin/bash

if [[ $1 == "" ]];then
   echo '-----------------------------------------------------'
   echo '          ERROR: NUMBER OF ATOMS UNDEFINED           ' 
   echo '-----------------------------------------------------'
   echo 'Usage:                                               '
   echo 'bash buildqff_gaussian <nat> <stencil_file> <restart>'
   echo '                                                     '
   echo '<nat> = # of atoms                                   '
   echo '<stencil_file> = File with stencil geoms (qumvia out)'
   echo '<restart> = 0 (no) or 1 (yes)                        '
   echo '                                                     '
   echo 'IMPORTANT:                                           '
   echo 'Gaussian formated checkpoint file of equilibrium     '
   echo 'geometry freq calculation must be renamed geo0.fchk  '
   echo 'and placed in the same directory this script runs.   '
   echo 'Formated check file can be obtained using formchk    '
   echo 'program included in the normal gaussian distribution.'
   echo 'Also if restart=1, the script will ask the user to   '
   echo 'define which step to restart from.                   '
   echo '-----------------------------------------------------'
   exit
fi

#-----------------------------------------------------------
# READING VARIABLES FROM COMMAND LINE

nat=$1              # Number of atoms
stencil_file=$2     # File with stencil geometries
dorestart=$3        # Is this a restart? yes=1, no=0.
#-----------------------------------------------------------

#-----------------------------------------------------------
# IF THIS IS A RESTART RUN THE USER IS ASKED TO DEFINE THE 
# INITIAL STEP. 
# IT SHOULD BE ANY NUMBER BETWEEN 1 AND 2*(3N-6)
if [ ${dorestart} -eq 0 ]; then
	inistep=1           # If yes, which step to start from?
else
	echo 'Please introduce the step to restart from \n'
        echo '>>> '
        read inistep
fi
#-----------------------------------------------------------

# DEFINING SOME SUBSIDIARY VARIABLES
nlne=$(($nat+1))    # Number of lines to read in geoms.qva
nvdf=$(($nat*3-6))  # Number of vibrational degrees of freedom.
ntot=$(($nvdf*2))   # Total number of hessians to be computed.

# ECHO INPUT 
echo "NUMBER OF ATOMS ${nat}"
echo "STENCIL FILE NAME: ${stencil_file}"
if [ $dorestart -eq 0 ];then
	echo "RESTART?: NO"
else
	echo "RESTART?: YES"
        echo "INITIAL STEP: $inistep"
fi

# IF THIS IS NOT A RESTART CLEAN UP FILES FROM PREVIOUS RUNS.
if [ ${dorestart} -eq 0 ]; then
#       DELETE FILES OF PREVIOUS RUNS IF PRESENT
	rm -f hess.gau geo*.inp
	rm -rf gaufiles
	mkdir gaufiles
	cp geo0.fchk ${stencil_file} gaufiles
fi

cd gaufiles


# NOW, LET'S GET SOME WORK DONE, BITCHES!

# BEGINNING MAIN LOOP
i=${inistep}
# 
while [[ i -le $ntot ]];do
# CREATING TEMPLATE FOR GAMESS INPUT FILE.
rm -f template.inp
cat>>template.inp<<eof
%NProcShared=4
%mem=7000MB
%Chk=geo${i}.chk
# MP2/cc-pvdz freq=savenormalmodes NoSymm

eof
   # Creating input file
   cat template.inp > geo${i}.inp
   b=$(($i-1))
   lne=1
   while [[ lne -le $nlne ]];do
      line=$(($b*$nlne+$lne))
      if [[ lne -eq 1 ]] ;then
         awk -v l=$line 'NR==l {print "LABEL" $0}' ${stencil_file} >> geo${i}.inp
         label=$(awk -v l=$line 'NR==l {print "LABEL" $0}' ${stencil_file})
         echo "" >> geo${i}.inp
         echo "0 1" >> geo${i}.inp
      else
         awk -v l=$line 'NR==l {print}' ${stencil_file} >> geo${i}.inp
      fi
      let lne=lne+1
   done
   echo '' >> geo${i}.inp

   g09 geo${i}.inp &> geo${i}.log
   formchk geo${i}.chk
   
   echo "$label" >> hess.gau
   awk 'BEGIN {switch=0} 
        {
            if (switch==1) {
               if (/Dipole Moment/) exit
               print
            }
            if (/Cartesian Force Constants/) switch=1
        }' geo${i}.fchk >> hess.gau
   echo "FINISHED $i/$ntot HESSIAN CALCULATIONS"
   
   let i=i+1
done

rm -f tmp tmp2
echo "LABEL  0  0    0" >> tmp
awk 'BEGIN {switch=0} 
     {
         if (switch==1) {
            if (/Dipole Moment/) exit
            print
         }
         if (/Cartesian Force Constants/) switch=1
     }' geo0.fchk >> tmp
cat tmp hess.gau > tmp2
mv tmp2 hess.gau
rm -f tmp
cp hess.gau ..

echo '###########################################'
echo '     FINISHED HESSIAN QFF CALCULATION      '
echo 'Results are stored in file hess.gau for    '
echo 'subsequent calculation of qff in QUMVIA.   '
echo '###########################################'

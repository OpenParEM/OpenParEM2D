#!/bin/bash

# no arguments

testCaseFile="regression_case_list.txt"
regressionResultsFile="regression_results.csv"
currentDirectory=$PWD

# start time
start_time=$SECONDS

#--------------------------------------------------------------------------
# input processing
#--------------------------------------------------------------------------

# check for arguments
if [ $# -gt 0 ]
then
   echo "ERROR: process.sh has command line arguments"
   exit 1
fi

#--------------------------------------------------------------------------
# data processing
#--------------------------------------------------------------------------

# see if the list of projects to process exists
if [[ ! -f $testCaseFile ]]
then
   echo "ERROR: Test case list file \""$testCaseFile"\" does not exist."
   exit 1
fi

# save the old data just in case
if [[ -f $regressionResultsFile ]]
then
   echo "#Moving "$regressionResultsFile" to "$regressionResultsFile".archive"
   echo ""
   mv $regressionResultsFile $regressionResultsFile".archive"
fi

# loop for each line in testCaseFile and process
while read -r line
do

   # skip on comment
   if [[ ${line::1} == "#" ]]
   then
      continue
   fi

   # skip on blank line
   if [[ $line == "" ]]
   then
      continue
   fi

   # split on " "
   IFS=' ' read -ra arr1 <<< "$line"
   numProc=${arr1[1]}

   # split on "/"
   IFS='/' read -ra arr <<< ${arr1[0]}

   # assemble the path and project
   projectPath=${arr[0]}"/"
   i=1;
   ilimit=$(( ${#arr[@]}-1 ))
   while [ $i -lt $ilimit ]
   do
      projectPath+=${arr[$i]}"/"
      i=$(( $i + 1 ))
   done
   projectFile=${arr[$i]}
   projectName="${projectFile%.*}"

   # add some labeling to help find things if necessary
   echo "#Executing "${arr1[0]}
   echo "#"$line >> $regressionResultsFile

   # change to the project directory
   cd $projectPath

   # process
   process.sh $projectFile $numProc < /dev/null

   # accumulate the results
   cat $projectName"_test_test_results.csv" >> $currentDirectory"/"$regressionResultsFile

   # reset back to the working directory
   cd $currentDirectory
   echo ""

done < $testCaseFile

# finish up with elapsed time
end_time=$SECONDS
elapsed=$(( end_time - start_time ))
echo "Regression elapsed time (s): "$elapsed




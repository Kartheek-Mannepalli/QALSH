#!/bin/bash

make
rm -rf ./src/*.o 2> /dev/null
rm -rf ./src/*.d 2> /dev/null

# Parameters of QALSH:
#    -alg  (integer)   options of algorithms (0 - 3)
#    -d    (integer)   dimensionality of the dataset
#    -n    (integer)   cardinality of the dataset
#    -qn   (integer)   number of queries
#    -B    (integer)   page size
#    -c    (real)      approximation ratio (c > 1)
#    -ds   (string)    file path of the dataset
#    -qs   (string)    file path of the query set
#    -ts   (string)    file path of the ground truth set
#    -of   (string)    output folder to store info of qalsh
#
# The options of algorithms (-alg) are:
#    0 - Ground-Truth
#        Parameters: -alg 0 -n -qn -d -ds -qs -ts
#
#    1 - Indexing
#        Parameters: -alg 1 -n -d -B -c -ds -of
#
#    2 - QALSH
#        Parameters: -alg 2 -qn -d -qs -ts -of
#
#    3 - Linear Scan;
#        Parameters: -alg 3 -n -qn -d -B -qs -ts -of
#
# NOTE: Each parameter is required to be separated by one space


names=(Sift)
card=(1000000)
dim=(128)

i=0
while [ $i -lt 1 ]
do

        echo "##############################"
        echo "#           QALSH            #"
        echo "##############################"

        algorithm=QALSH
        dsName=${names[i]}

        epochTime=`date +%s`

        ############ Parameters ############
        dsPath=~/Desktop/LSH/Datasets/${dsName}/${dsName}.ds
        resultsPath=~/Desktop/LSH/Results/${algorithm}/${dsName}
        indexPath=~/Desktop/LSH/Indexes/${algorithm}/${dsName}/
        dataBinPath=~/Desktop/LSH/Datasets/${dsName}/${dsName}.bin

        nPoints=${card[i]}
        nAttributes=${dim[i]}

        isQinDS=1
        pageSize=4096

	      qName=Sift_7789
        nQueries=1
        qsPath=~/Desktop/LSH/Queries/${dsName}/${qName}.qs
        gtPath=~/Desktop/LSH/Ground_truth/${dsName}/${qName}

        ####################################

        mkdir -p ${indexPath}/indexes/
        mkdir -p ${resultsPath}/linear/
        mkdir -p ${gtPath}/


        echo "\n########## Indexing ##########"
        ./qalsh -alg 1 -n ${nPoints} -d ${nAttributes} -B ${pageSize} -c 2.0 -ds ${dsPath} -of ${resultsPath} \
        -indexPath ${indexPath} -dataBinPath ${dataBinPath} \
        |& tee ${indexPath}/indexing.out


        echo "\n### Ground Truth Creation ###"
        ./qalsh -alg 0 -n ${nPoints} -qn ${nQueries} -d ${nAttributes} -ds ${dsPath} -qs ${qsPath} -ts ${gtPath}/gt.gt2 \
         -indexPath ${indexPath} -isQinDS ${isQinDS} \
        |& tee ${gtPath}/ground_truth.out


        echo "\n########### QALSH ###########"
        ./qalsh -alg 2 -qn ${nQueries} -d ${nAttributes} -qs ${qsPath} -ts ${gtPath}/gt.gt2 -of ${resultsPath} \
        -indexPath ${indexPath} -dataBinPath ${dataBinPath} -isQinDS ${isQinDS} \
        |& tee ${resultsPath}/lsh.out

#	echo "########### Linear Search ###########"
#	./qalsh -alg 3 -n ${nPoints} -qn ${nQueries} -d ${nAttributes} -B ${pageSize} -ds ${dsPath} -qs ${qsPath} \
#		-ts ${gtPath}/gt.gt2 \
#		-of ${resultsPath}/linear \
#	     	|& tee ${resultsPath}/linear/lsh.out

	echo "Password0" | sudo -S bash clear_cache.sh

    i=$[$i+1]
done


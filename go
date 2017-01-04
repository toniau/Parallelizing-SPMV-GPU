#
# Input file declarations
#
matrixFile="../proj4.testcases/turon_m.mtx"
vectorFile="../proj4.testcases/turon_m_vector.txt"

matrixFile="../proj4.testcases/watson_2.mtx"
vectorFile="../proj4.testcases/watson_2_vector.txt"

matrixFile="../proj4.testcases/mac_econ_fwd500.mtx"
vectorFile="../proj4.testcases/mac_econ_fwd500_vector.txt"

matrixFile="../proj4.testcases/scircuit.mtx"
vectorFile="../proj4.testcases/scircuit_vector.txt"

matrixFile="../proj4.testcases/cop20k_A.mtx"
vectorFile="../proj4.testcases/cop20k_A_vector.txt"



##############################################################
#
#
# Run multiple atomic adds, each with a different block size
#
threadsPerBlock="128"
threadsPerBlock="1 2 4 8 16 32 64 128 256 384 512 1024"

algorithm="atom"

echo
echo -e "----------  Atomic Add algorithm  ------------"
echo
for t in $threadsPerBlock
do
     rm -f ./*.txt
echo spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./helper/verify $matrixFile  ./output.txt
     echo 

done




##############################################################
#
#
#
# Run multiple segmented scans, each with a different block size
#
threadsPerBlock="128"
threadsPerBlock="2 4 8 16 32 64 128 256 384 512 1024"

algorithm="scan"

echo
echo -e "----------  Segmented Scan algorithm  ------------"
echo
for t in $threadsPerBlock
do
     rm -f ./*.txt
echo spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./helper/verify $matrixFile  ./output.txt
     echo 

done




##############################################################
#
#
#
# Run multiple optimized segmented scans, each with a different block size
#
threadsPerBlock="128"
threadsPerBlock="2 4 8 16 32 64 128 256 384 512 1024"

algorithm="optscan"

echo
echo -e "----------  Optimized Segmented Scan algorithm  ------------"
echo
for t in $threadsPerBlock
do
     rm -f ./*.txt
echo spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./spmv -mat $matrixFile -ivec $vectorFile -alg $algorithm -blocksize $t
     ./helper/verify $matrixFile  ./output.txt
     echo 

done


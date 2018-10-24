#for degree in 4 5
#do
#	echo "processing degree = $degree for naive protocol\n"
#	./main-naive data/test_${degree}Degree10AS 10 3 10 > result_${degree}Degree10AS_naive
#	echo "processing degree = $degree for sdni protocol\n"
#	./main-sdni data/test_${degree}Degree10AS 10 3 10 > result_${degree}Degree10AS_sdni
#done

for load in 0.1 0.25 0.5 0.65 0.8 1 1.2 1.4 1.6 2 3 5 8 10 13 15
do
	echo "processing load = $load for naive protocol"
	./main-naive data/test_6Degree10AS 10 3 $load > result_load${load}_naive
	echo "processing load = $load for sdni protocol"
	./main-sdni data/test_6Degree10AS 10 3 $load > result_load${load}_sdni
done

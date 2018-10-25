for load in 0.1 0.25 0.5 0.65 0.8 1 1.2 1.4 1.6 2 3 5 8 13 15
do
	echo "processing load = $load for maxthr algorithm"
	./main-benchmark data/test_5Degree10AS_biodirection 10 $load > result_load${load}_benchmark_biodirection
done

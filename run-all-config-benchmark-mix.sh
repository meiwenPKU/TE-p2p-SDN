for load in 0.65 0.8 1 1.2 1.4 1.6 2 3 5 8 10 13 15
do
	echo "processing load = $load for mix algorithm"
	./main-benchmark-mix data/test_5Degree10AS_biodirection 10 $load > result_load${load}_benchmark_mix_biodirection
done

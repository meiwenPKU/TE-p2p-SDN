for load in 0.5 0.65 0.8 1 1.2 1.4 1.6 2 3 5 8 10 13 15
do
	echo "processing load = $load for google algorithm"
	./main-naive-google data/test_5Degree10AS_biodirection 10 3 $load > result_load${load}_naive_google_biodirection
done
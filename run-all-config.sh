# for degree in 3 4 5 6 7 8 9
# do
# 	echo "processing degree = $degree for naive protocol\n"
# 	./main-naive data/test_${degree}Degree10AS_biodirection 10 3 10 > result_${degree}Degree10AS_naive_biodirection
# 	echo "processing degree = $degree for sdni protocol\n"
# 	./main-sdni data/test_${degree}Degree10AS_biodirection 10 3 10 > result_${degree}Degree10AS_sdni_biodirection
# done

# for load in 0.1 0.25 0.5 0.65 0.8 1 1.2 1.4 1.6 2 3 5 8 10 13 15
# do
# 	echo "processing load = $load for naive protocol"
# 	./main-naive data/test_6Degree10AS 10 3 $load > result_load${load}_naive
# 	echo "processing load = $load for sdni protocol"
# 	./main-sdni data/test_6Degree10AS 10 3 $load > result_load${load}_sdni
# done

for kpath in 4 5 6 7
do
	echo "processing kpath = $kpath for naive protocol\n"
	./main-naive data/test_5Degree10AS_biodirection 10 $kpath 10 > result_kpath${kpath}_naive_biodirection
	echo "processing kpath = $kpath for sdni protocol\n"
	./main-sdni data/test_5Degree10AS_biodirection 10 $kpath 10 > result_kpath${kpath}_sdni_biodirection
done

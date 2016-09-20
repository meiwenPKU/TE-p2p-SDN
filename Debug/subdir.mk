################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DijkstraShortestPathAlg.cpp \
../FordFulkersonAlg.cpp \
../GlobalVariable.cpp \
../Graph.cpp \
../MainP.cpp \
../TopoTable.cpp \
../YenTopKShortestPathsAlg.cpp 

OBJS += \
./DijkstraShortestPathAlg.o \
./FordFulkersonAlg.o \
./GlobalVariable.o \
./Graph.o \
./MainP.o \
./TopoTable.o \
./YenTopKShortestPathsAlg.o 

CPP_DEPS += \
./DijkstraShortestPathAlg.d \
./FordFulkersonAlg.d \
./GlobalVariable.d \
./Graph.d \
./MainP.d \
./TopoTable.d \
./YenTopKShortestPathsAlg.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



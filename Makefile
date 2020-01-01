.PHONY: debug

driver: saye_utils.o saye_algorithm.o driver.o
	g++ -g -std=c++11 driver.o saye_utils.o saye_algorithm.o -o driver.x

surftest: saye_utils.o saye_algorithm.o surftest.o
	g++ -g -std=c++11 surftest.o saye_utils.o saye_algorithm.o -o surftest.x

costtest: saye_utils.o saye_algorithm.o costtest.o grad_descent.o
	g++ -g -std=c++11 costtest.o saye_utils.o saye_algorithm.o grad_descent.o -o costtest.x

debug: saye_utils.o saye_algorithm.o driver.o
	g++ -g -Wall -std=c++11 driver.o saye_utils.o saye_algorithm.o -o driver.x

saye_utils.o: saye_utils.cpp saye_utils.h
	g++ -g -std=c++11 -c saye_utils.cpp

saye_algorithm.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h base_functions.h
	g++ -g -std=c++11 -c saye_algorithm.cpp

grad_descent.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h grad_descent.cpp grad_descent.h
	g++ -g -std=c++11 -c grad_descent.cpp

driver.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h driver.cpp
	g++ -g -std=c++11 -c driver.cpp

surftest.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h surftest.cpp
	g++ -g -std=c++11 -c surftest.cpp

costtest.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h costtest.cpp grad_descent.cpp grad_descent.h
	g++ -g -std=c++11 -c costtest.cpp

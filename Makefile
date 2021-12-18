CC = g++ 
CFLAGS = -g -Wall -std=c++11

.PHONY: debug

driver: saye_utils.o saye_algorithm.o level_set_plot.o rbf_network.o grad_descent.o driver.o
	g++ -g -std=c++11 driver.o saye_utils.o saye_algorithm.o rbf_network.o grad_descent.o level_set_plot.o -larmadillo -I/usr/include/python2.7 -lpython2.7 -o driver.x

total_test_quad: saye_utils.o saye_algorithm.o level_set_plot.o rbf_network.o grad_descent.o total_test_quad.o
	$(CC) -std=c++11 total_test_quad.o grad_descent.o saye_utils.o saye_algorithm.o rbf_network.o -larmadillo level_set_plot.o -I/usr/include/python2.7 -lpython2.7 -o total_test_quad.x

total_test_para: saye_utils.o saye_algorithm.o level_set_plot.o rbf_network.o grad_descent.o total_test_para.o
	$(CC) -std=c++11 total_test_para.o grad_descent.o saye_utils.o saye_algorithm.o rbf_network.o -larmadillo level_set_plot.o -I/usr/include/python2.7 -lpython2.7 -o total_test_para.x

driver.o: saye_algorithm.cpp saye_algorithm.h saye_utils.cpp saye_utils.h driver.cpp
	g++ -g -std=c++11 -c driver.cpp

total_test_quad.o: total_test_quad.cpp
	$(CC) $(CFLAGS) -c total_test_quad.cpp

total_test_para.o: total_test_para.cpp
	$(CC) $(CFLAGS) -c total_test_para.cpp

grad_descent.o: grad_descent.cpp grad_descent.h
	$(CC) $(CFLAGS) -c grad_descent.cpp

level_set_plot.o: level_set_plot.cpp level_set_plot.h
	$(CC) $(CFLAGS) -c level_set_plot.cpp -I/usr/include/python2.7 -lpython2.7

rbf_network.o: rbf_network.cpp rbf_network.h
	$(CC) $(CFLAGS) -c rbf_network.cpp

saye_utils.o: saye_utils.cpp saye_utils.h
	$(CC) $(CFLAGS) -c saye_utils.cpp

saye_algorithm.o: saye_algorithm.cpp saye_algorithm.h
	$(CC) $(CFLAGS) -c saye_algorithm.cpp

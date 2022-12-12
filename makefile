make cn:
	g++ main.cpp src/cn.cpp -std=c++11 -I include -larmadillo  -o main
	./main

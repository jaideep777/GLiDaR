#ifndef _TIMER_FUNCTION_LINUX_H_
#define _TIMER_FUNCTION_LINUX_H_

#include <chrono>
#include <utility>
#include <iostream>

template<typename F, typename... Args>
void funcTime(F func, Args&&... args){
	std::chrono::_V2::system_clock::time_point t1, t2;
	t1 = std::chrono::high_resolution_clock::now();
    func(std::forward<Args>(args)...);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "[" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms]\n";
}

class MeasureTime {
	private:
    std::chrono::time_point<std::chrono::high_resolution_clock>  _start;

    public:
	MeasureTime() : _start(std::chrono::high_resolution_clock::now()) {}

    ~MeasureTime() {
        auto stop = std::chrono::high_resolution_clock::now();
		std::cout << "[" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - _start).count() << " ms]\n";
    }
};

#endif

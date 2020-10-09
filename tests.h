#pragma once

#include "flow_algorithms.h"
#include "my_random.h"

#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

class LogDuration {
public:
    explicit LogDuration(const std::string& msg = "")
        : message(msg + ": ")
        , start(steady_clock::now())
    {
    }

    ~LogDuration() {
        auto finish = steady_clock::now();
        auto dur = finish - start;
        std::cerr << message
            << duration_cast<milliseconds>(dur).count()
            << " ms" << std::endl;
    }
private:
    std::string message;
    steady_clock::time_point start;
};

#define UNIQ_ID_IMPL(lineno) _a_local_var_##lineno
#define UNIQ_ID(lineno) UNIQ_ID_IMPL(lineno)

#define LOG_DURATION(message) \
    LogDuration UNIQ_ID(__LINE__){message};

template<class T, class U>
void AssertEqual(const T& t, const U& u, const string& hint = {}) {
    if (!(t == u)) {
        ostringstream os;
        os << "Assertion failed: " << t << " != " << u;
        if (!hint.empty()) {
            os << " hint: " << hint;
        }
        throw runtime_error(os.str());
    }
}

inline void Assert(bool b, const string& hint) {
    AssertEqual(b, true, hint);
}

class TestRunner {
public:
    template <class TestFunc>
    void RunTest(TestFunc func, const string& test_name) {
        try {
            func();
            cerr << test_name << " OK" << endl;
        } catch (exception& e) {
            ++fail_count;
            cerr << test_name << " fail: " << e.what() << endl;
        } catch (...) {
            ++fail_count;
            cerr << "Unknown exception caught" << endl;
        }
    }

    ~TestRunner() {
        if (fail_count > 0) {
            cerr << fail_count << " unit tests failed. Terminate" << endl;
            exit(1);
        }
    }

private:
    int fail_count = 0;
};

#define ASSERT_EQUAL(x, y) {            \
    ostringstream os;                     \
    os << #x << " != " << #y << ", "      \
        << __FILE__ << ":" << __LINE__;     \
    AssertEqual(x, y, os.str());          \
}

#define ASSERT(x) {                     \
    ostringstream os;                     \
    os << #x << " is false, "             \
        << __FILE__ << ":" << __LINE__;     \
    Assert(x, os.str());                  \
}

#define RUN_TEST(tr, func) \
    tr.RunTest(func, #func)

template<typename T>
T randomBetween(T l, T r) {
	return randomDouble01() * (r - l) + l;
}

template<typename T>
EdgesContainer<T> getRandomEdgesContainer(int n, int m, T min_capacity, T max_capacity) {
	EdgesContainer<T> result(n);
	result.setSource(0);
	result.setSink(n - 1);

	double edge_probability = min(1.0, static_cast<double>(m) / (n * (n - 1) / 2));
	for (int from = 0; from < n; ++from) {
		for (int to = 0; to < n; ++ to) {
			if (from == to) {
				continue;
			}
			if (randomDouble01() < edge_probability) {
				T capacity = randomBetween<T>(min_capacity, max_capacity);
				result.addEdge(from, to, capacity);
			}
		}
	}

	return result;
}

template<typename T, typename A>
void runTimeTest(EdgesContainer<T>& edges_container, A getFlowAlgorithm, const string& message) {
	{	LOG_DURATION(message)
		getFlowAlgorithm(edges_container);
	}
}

template<typename T, typename A>
void runTimeTest(int n, int m, T min_capacity, T max_capacity,
		A getFlowAlgorithm, const string& message) {
	EdgesContainer<T> edges_container = getRandomEdgesContainer<T>(n, m, min_capacity, max_capacity);
	{	LOG_DURATION(message)
		getFlowAlgorithm(edges_container);
	}
}

template<typename T>
class TestsContainer {
public:
	int getTestsNumber() const { return tests.size(); }
	void addTest(const EdgesContainer<T>& edges_container) { tests.push_back(edges_container); }
	EdgesContainer<T> getEdgesContainer(int test_number) const { 	return tests[test_number]; }
private:
	vector<EdgesContainer<T>> tests;
};


#include "flow_algorithms.h"
#include "tests.h"

#include <bits/stdc++.h>

using namespace std;

static const int tests_number = 20;
static const int min_n = 50;
static const int delta_n = 50;
static const long long min_capacity = 2;
static const long long max_capacity = 1000;

int getEdgesNumberForTest(int n) {
	return n * 5;
}

template<typename T>
int getEdgesNumber(EdgesContainer<T>& edges_container) {
	int n = edges_container.getVertexNumber();
	int answer = 0;
	for (int from = 0; from < n; ++from) {
		for (typename EdgesContainer<T>::EdgeIterator it = edges_container.Begin(from);
				it.isValid(); ++it) {
			++answer;
		}
	}
	return answer;
}

string ToString(int n) {
	string answer;
	while(n > 0) {
		char last_digit = (n % 10) + '0';
		string last_digit_str;
		last_digit_str += last_digit;
		n /= 10;
		answer = last_digit + answer;
	}
	return answer;
}

template<typename T>
string getTestInfo(EdgesContainer<T>& edges_container) {
    int n = edges_container.getVertexNumber();
    int m = getEdgesNumber<T>(edges_container);
    return ToString(n) + " " + ToString(m);
}

template<typename A, typename T>
void runTimeTests(const TestsContainer<T>& tests_container, A getFlowAlgorithm) {
	for (int i = 0; i < tests_container.getTestsNumber(); ++i) {
		EdgesContainer<T> edges_container = tests_container.getEdgesContainer(i);
		string message = getTestInfo<T>(edges_container);
		runTimeTest(edges_container, getFlowAlgorithm, message);
	}
}

template<typename T>
TestsContainer<T> getTestsContainer() {
	TestsContainer<T> tests_container;
	for (int test_number = 0; test_number < tests_number; ++test_number) {
		int n = min_n + test_number * delta_n;
		int m = getEdgesNumberForTest(n);
		tests_container.addTest(getRandomEdgesContainer(n, m, min_capacity, max_capacity));
	}
	return tests_container;
}

template<typename T>
void runTests() {
	TestsContainer<T> tests_container = getTestsContainer<T>();
	cerr << "Dinitz time test:" << endl;
	runTimeTests(tests_container, getMaxFlowDinitz<T>);
	cerr << "Push & relabel time test:" << endl;
	runTimeTests(tests_container, getMaxFlowPushRelabel<T>);
	cerr << "Malhotra-Kumar-Maheshwari time test:" << endl;
	runTimeTests(tests_container, getMaxFlowMalhotraKumarMaheshwari<T>);
}

int main() {
	runTests<long long>();
	return 0;
}

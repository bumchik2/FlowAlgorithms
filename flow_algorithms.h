#pragma once

#include "flow_algorithms.h"

#include <iostream>
#include <vector>
#include <stdexcept>
#include <list>
#include <queue>
#include <algorithm>
#include <string>

using std::min;
using std::max;
using std::queue;
using std::list;
using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::runtime_error;
using std::istream;
using std::ostream;
using std::string;

const long long INF_long_long = 100'000'000'000;
const long long INF_long_long2 = 1'000'000'000'000'000;
const int INF_int = 1'000'000'000;

template<typename T>
struct Edge {
public:
    int from;
    int to;
    T capacity;
    T flow = T(0);

    bool deleted = false;
    // this field is used in Malhotra-Kumar-Maheshwari algorithm.
    // it would be more correct to make an external vector<bool> for
    // storing such information, but that would make code more complicated

    T remainingCapacity() const;
    bool isFilled() const;
    void print() const;
};

template<typename T>
T Edge<T>::remainingCapacity() const {
    return capacity - flow;
}

template<typename T>
bool Edge<T>::isFilled() const {
    return remainingCapacity() <= 0;
}

template<typename T>
void Edge<T>::print() const {
    cout << "from: " << from << "; to: " << to << "; c: " <<
            capacity <<"; f: " << flow << endl;
}

template<typename T>
class EdgesContainer {
public:
    EdgesContainer(int n): n_(n), ends_(vector<int>(n, -1)) {}

    void setSource(int s) { s_ = s; }
    int getSource() const { return s_; }
    void setSink(int t) { t_ = t; }
    int getSink() const { return t_; }
    int getVertexNumber() const { return n_; }

    void addEdge(int from, int to, T capacity, T back_capacity = 0);
    void pushFlow(int edge_number, T flow);

    class EdgeIterator {
    private:
        EdgesContainer* edges_container_ptr_ = nullptr;
        int edge_number_ = -1;
    public:
        EdgeIterator() = default;
        EdgeIterator(EdgesContainer* edges_container_ptr, int edge_number):
            edges_container_ptr_(edges_container_ptr), edge_number_(edge_number) {}
        EdgeIterator(const EdgeIterator&) = default;
        EdgeIterator& operator = (const EdgeIterator&) = default;

        EdgeIterator& operator ++();
        bool isValid();

        void checkValid(const string& hint);
        Edge<T> getEdge();
        Edge<T> getBackwardEdge();
        void pushFlow(T flow);
        void pushBackwardFlow(T flow);
        void setDeleted(bool deleted);

        template<typename U>
        friend class EdgesContainer;
    };

    EdgeIterator Begin(int vertex_number);

    T getCurrentFlow();
    void print() const;

private:
    int n_;
    int s_ = -1;
    int t_ = -1;

    vector<Edge<T>> edges_;
    vector<int> ends_;
    vector<int> prevs_;

    int getReversedEdgeNumber_(int edge_number);
    void addSingleDirectionEdge_(int from, int to, T capacity);
    void pushSingleDirectionFlow_(int edge_number, T flow);
};

template<typename T>
typename EdgesContainer<T>::EdgeIterator EdgesContainer<T>::Begin(int vertex_number) {
    return EdgeIterator(this, ends_[vertex_number]);
}

template<typename T>
void EdgesContainer<T>::addEdge(int from, int to, T capacity, T back_capacity) {
    addSingleDirectionEdge_(from, to, capacity);
    addSingleDirectionEdge_(to, from, back_capacity);
}

template<typename T>
void EdgesContainer<T>::pushFlow(int edge_number, T flow) {
    pushSingleDirectionFlow_(edge_number, flow);
    pushSingleDirectionFlow_(getReversedEdgeNumber_(edge_number), -flow);
}

template<typename T>
T EdgesContainer<T>::getCurrentFlow() {
    T answer = T(0);
    for (EdgeIterator it = Begin(s_); it.isValid(); ++it) {
        Edge<T> current_edge = it.getEdge();
        answer += current_edge.flow;
    }
    return answer;
}

template<typename T>
void EdgesContainer<T>::print() const {
    cout << "Edges container info:" << endl;
    cout << "n: " << n_ << "; s: " << s_ << "; t: " << t_ << endl;
    cout << "edges: " << endl;
    for (unsigned int i = 0; i < edges_.size(); ++i) {
        edges_[i].print();
    }
}

template<typename T>
bool EdgesContainer<T>::EdgeIterator::isValid() {
    return edges_container_ptr_ != nullptr && edge_number_ != -1;
}

template<typename T>
void EdgesContainer<T>::EdgeIterator::checkValid(const std::string& hint) {
    if (!isValid()) {
        throw runtime_error(hint);
    }
}

template<typename T>
typename EdgesContainer<T>::EdgeIterator& EdgesContainer<T>::EdgeIterator::operator ++() {
    checkValid("trying to increment invalid iterator");
    edge_number_ = edges_container_ptr_->prevs_[edge_number_];
    return *this;
}

template<typename T>
Edge<T> EdgesContainer<T>::EdgeIterator::getEdge() {
    checkValid("trying to call getEdge from invalid iterator");
    return edges_container_ptr_->edges_[edge_number_];
}

template<typename T>
Edge<T> EdgesContainer<T>::EdgeIterator::getBackwardEdge() {
    checkValid("trying to call getBackwardEdge from invalid iterator");
    int backward_edge_number = edges_container_ptr_->getReversedEdgeNumber_(edge_number_);
    return edges_container_ptr_->edges_[backward_edge_number];
}

template<typename T>
void EdgesContainer<T>::EdgeIterator::pushFlow(T flow) {
    checkValid("trying to push flow in invalid edge");
    edges_container_ptr_->pushFlow(edge_number_, flow);
}

template<typename T>
void EdgesContainer<T>::EdgeIterator::pushBackwardFlow(T flow) {
    checkValid("trying to push backward flow in invalid edge");
    int backward_edge_number = edges_container_ptr_->getReversedEdgeNumber_(edge_number_);
    edges_container_ptr_->pushFlow(backward_edge_number, flow);
}

template<typename T>
void EdgesContainer<T>::EdgeIterator::setDeleted(bool deleted) {
    checkValid("trying to delete invalid edge");
    edges_container_ptr_->edges_[edge_number_].deleted = deleted;
    int backward_edge_number = edges_container_ptr_->getReversedEdgeNumber_(edge_number_);
    edges_container_ptr_->edges_[backward_edge_number].deleted = deleted;
}

template<typename T>
int EdgesContainer<T>::getReversedEdgeNumber_(int edge_number) {
    return edge_number ^ 1;
}

template<typename T>
void EdgesContainer<T>::addSingleDirectionEdge_(int from, int to, T capacity) {
    edges_.push_back(Edge<T>{from, to, capacity, 0});
    prevs_.push_back(ends_[from]);
    ends_[from] = edges_.size() - 1;
}

template<typename T>
void EdgesContainer<T>::pushSingleDirectionFlow_(int edge_number, T flow) {
    edges_[edge_number].flow += flow;
}

/// EdgesContainer part is over

template<typename T>
void preflowBFS(EdgesContainer<T>& edges_container, vector<int>& time_in) {
    int n = edges_container.getVertexNumber();
    int s = edges_container.getSource();

    queue<int> q;
    q.push(s);

    vector<bool> used(n, false);
    used[s] = true;

    time_in = vector<int>(n, INF_int);
    time_in[s] = 0;

    while(!q.empty()) {
        int from = q.front();
        for (typename EdgesContainer<T>::EdgeIterator it =
                edges_container.Begin(from); it.isValid(); ++it) {
            Edge<T> current_edge = it.getEdge();
            int to = current_edge.to;
            if (used[to] || current_edge.isFilled()) {
                continue;
            }
            q.push(to);
            time_in[to] = time_in[from] + 1;
            used[to] = true;
        }
        q.pop();
    }
}

template<typename T>
bool dinitzDFS(EdgesContainer<T>& edges_container, const vector<int>& time_in,
        vector<bool>& used, int vertex_number, T local_min, T& global_min,
        vector<typename EdgesContainer<T>::EdgeIterator>& iterators) {
    used[vertex_number] = true;

    int t = edges_container.getSink();
    if (vertex_number == t) {
        return true;
    }

    for (typename EdgesContainer<T>::EdgeIterator it = iterators[vertex_number];
            it.isValid(); ++it) {
        Edge<T> current_edge = it.getEdge();
        int to = current_edge.to;
        if (time_in[to] != time_in[vertex_number] + 1 || used[to] ||
                    current_edge.isFilled()) {
            ++iterators[vertex_number];
            continue;
        }

        T new_local_min = min(local_min, current_edge.remainingCapacity());
        bool managed_to_find_path = dinitzDFS(edges_container, time_in,
                used, to, new_local_min, global_min, iterators);
        if (managed_to_find_path) {
            global_min = min(global_min, new_local_min);
            it.pushFlow(global_min);
            return true;
        }
    }

    return false;
}

template<typename T>
T getMaxFlowDinitz(EdgesContainer<T>& edges_container) {
    int n = edges_container.getVertexNumber();
    int s = edges_container.getSource();
    int t = edges_container.getSink();

    while(true) {
        vector<int> time_in;
        preflowBFS(edges_container, time_in);
        if (time_in[t] == INF_int) {
            break;
        }

        vector<typename EdgesContainer<T>::EdgeIterator> iterators(n);
        for (int i = 0; i < n; ++i) {
            iterators[i] = edges_container.Begin(i);
        }

        while(true) {
            T global_min = T(INF_long_long);
            vector<bool> used (n, false);
            bool managed_to_find_path = dinitzDFS(edges_container, time_in,
                    used, s, INF_long_long, global_min, iterators);
            if (!managed_to_find_path) {
                break;
            }
        }
    }

    return edges_container.getCurrentFlow();
}

template<typename T>
void addP(int from, int to, T delta, vector<T>& pin, vector<T>& pout, vector<T>& p) {
    pout[from] += delta;
    p[from] = min(pin[from], pout[from]);
    if (p[from] < 0) {
        throw runtime_error("from wtf");
    }
    pin[to] += delta;
    p[to] = min(pin[to], pout[to]);
    if (p[to] < 0) {
        throw runtime_error("to wtf");
    }
}

template<typename T>
void removeP(int from, int to, T delta, vector<T>& pin, vector<T>& pout, vector<T>& p) {
    addP(from, to, -delta, pin, pout, p);
}

template<typename T>
void markUseless(EdgesContainer<T>& edges_container, int vertex_number, vector<bool>& useless,
        vector<T>& pin, vector<T>& pout, vector<T>& p, const vector<int>& time_in) {
    if (p[vertex_number] != T(0)) {
        throw runtime_error("only vertex with p == 0 can become useless");
    }
    if (useless[vertex_number]) {
        return;
    }
    useless[vertex_number] = true;

    for (typename EdgesContainer<T>::EdgeIterator it = edges_container.Begin(vertex_number);
            it.isValid(); ++it) {
        Edge<T> forward_edge = it.getEdge();
        int to = forward_edge.to;
        if (useless[to]  || forward_edge.deleted || abs(time_in[to] - time_in[vertex_number]) != 1) {
            continue;
        }
        Edge<T> edge = (time_in[to] > time_in[vertex_number]) ? it.getEdge() : it.getBackwardEdge();

        it.setDeleted(true);
        removeP(edge.from, edge.to, edge.remainingCapacity(), pin, pout, p);
        if (p[to] == 0) {
            markUseless(edges_container, to, useless, pin, pout, p, time_in);
        }
    }
}

template<typename T>
void pushSingleDirection(EdgesContainer<T>& edges_container, int vertex_number,
        T flow, vector<T>& pin, vector<T>& pout, vector<T>& p, vector<bool>& useless, const vector<int>& time_in,
        vector<typename EdgesContainer<T>::EdgeIterator>& iterators, bool pushing_forward) {
    if (vertex_number == edges_container.getSource() || vertex_number == edges_container.getSink()) {
        // we don't need to push from source or sink
        return;
    }

    for (typename EdgesContainer<T>::EdgeIterator it = iterators[vertex_number];
            it.isValid(); ++it) {
        Edge<T> edge = pushing_forward ? it.getEdge() : it.getBackwardEdge();
        int to = pushing_forward ? edge.to : edge.from;

        bool we_can_try_to_push = true;
        we_can_try_to_push &= !useless[to];
        we_can_try_to_push &= !edge.deleted;
        we_can_try_to_push &= !edge.isFilled();
        we_can_try_to_push &= (time_in[edge.to] == time_in[edge.from] + 1);
        if (!we_can_try_to_push) {
            ++iterators[vertex_number];
            continue;
        }

        T new_flow_to_push = min(edge.remainingCapacity(), flow);
        if (edge.remainingCapacity() <= flow) {
            ++iterators[vertex_number];
        }

        if (pushing_forward) {
            it.pushFlow(new_flow_to_push);
        } else {
            it.pushBackwardFlow(new_flow_to_push);
        }

        removeP(edge.from, edge.to, new_flow_to_push, pin, pout, p);

        if (new_flow_to_push == 0) {
            continue;
        }

        pushSingleDirection(edges_container, to, new_flow_to_push,
                pin, pout, p, useless, time_in, iterators, pushing_forward);
        flow -= new_flow_to_push;
        if (flow == T(0)) {
            break;
        }
    }
    return;
}

template<typename T>
void pushBothDirections(EdgesContainer<T>& edges_container, int vertex_number,
        T flow, vector<T>& pin, vector<T>& pout, vector<T>& p, vector<bool>& useless,
        const vector<int>& time_in, vector<typename EdgesContainer<T>::EdgeIterator>& forward_iterators,
        vector<typename EdgesContainer<T>::EdgeIterator>& backward_iterators) {
    pushSingleDirection(edges_container, vertex_number, flow, pin, pout,
            p, useless, time_in, forward_iterators, true); // pushing forward
    pushSingleDirection(edges_container, vertex_number, flow, pin, pout,
            p, useless, time_in, backward_iterators, false); // pushing backward
}

template<typename T>
void initializeMalhotraKumarMaheshwari(EdgesContainer<T>& edges_container,
        vector<T>& pin, vector<T>& pout, vector<T>& p, vector<bool>& useless, vector<int>& time_in,
        vector<typename EdgesContainer<T>::EdgeIterator>& forward_iterators,
        vector<typename EdgesContainer<T>::EdgeIterator>& backward_iterators,
        int n, int s, int t) {
    // initializing pin, pout, p
    for (int from = 0; from < n; ++from) {
        for (typename EdgesContainer<T>::EdgeIterator it = edges_container.Begin(from);
                it.isValid(); ++it) {
            it.setDeleted(false);
            Edge<T> current_edge = it.getEdge();
            int to = current_edge.to;
            if (!current_edge.isFilled() && time_in[to] == time_in[from] + 1) {
                addP(from, to, current_edge.remainingCapacity(), pin, pout, p);
            }
        }
    }
    pin[s] = pin[t] = pout[s] = pout[t] = p[s] = p[t] = T(INF_long_long2);

    // initializing useless vertexes
    for (int i = 0; i < n; ++i) {
        if (p[i] == 0) {
            markUseless(edges_container, i, useless, pin, pout, p, time_in);
        }
    }

    // initializing iterators
    for (int i = 0; i < n; ++i) {
        forward_iterators[i] = backward_iterators[i] = edges_container.Begin(i);
    }
}

template<typename T>
T getMaxFlowMalhotraKumarMaheshwari(EdgesContainer<T>& edges_container) {
    int n = edges_container.getVertexNumber();
    int s = edges_container.getSource();
    int t = edges_container.getSink();

    // if we don't fill direct s-->t edges, in some cases algortihm won't work correctly
    for (typename EdgesContainer<T>::EdgeIterator it = edges_container.Begin(s);
    		it.isValid(); ++it) {
    	Edge<T> edge = it.getEdge();
    	if (edge.to == t) {
        	it.pushFlow(edge.remainingCapacity());
    	}
    }

    while(true) {
        vector<int> time_in;
        preflowBFS(edges_container, time_in);
        if (time_in[t] == INF_int) {
            break;
        }

        for (int i = 0; i < n; ++i) {
            if (time_in[i] >= time_in[t] && i != t) {
                time_in[i] = INF_int; // we don't need these vertexes
            }
        }

        vector<T> pin (n, T(0));
        vector<T> pout (n, T(0));
        vector<T> p (n, T(0));
        vector<bool> useless(n, false);

        vector<typename EdgesContainer<T>::EdgeIterator> forward_iterators(n);
        vector<typename EdgesContainer<T>::EdgeIterator> backward_iterators(n);

        initializeMalhotraKumarMaheshwari(edges_container, pin, pout, p,
                useless, time_in, forward_iterators, backward_iterators, n, s, t);

        while(true) {
            T min_p = T(INF_long_long);
            int min_p_vertex = -1;
            for (int i = 0; i < n; ++i) {
                if (!useless[i] && p[i] < min_p && i != s && i != t) {
                    min_p_vertex = i;
                    min_p = p[i];
                }
            }

            if (min_p_vertex != -1) {
                pushBothDirections(edges_container, min_p_vertex, min_p,
                        pin, pout, p, useless, time_in, forward_iterators, backward_iterators);
                markUseless(edges_container, min_p_vertex, useless, pin, pout, p, time_in);
            } else {
                break;
            }
        }
    }

    return edges_container.getCurrentFlow();
}

/// Blocking flow algorithms part is over

template<typename T>
bool canPush(const vector<T> excess, const vector<int>& height, const Edge<T>& edge) {
    return excess[edge.from] > 0 && !edge.isFilled() &&
            height[edge.from] == height[edge.to] + 1;
}

template<typename T>
void push(EdgesContainer<T>& edges_container, vector<T>& excess,
        typename EdgesContainer<T>::EdgeIterator it) {
    Edge<T> edge = it.getEdge();
    T flow = min(excess[edge.from], edge.remainingCapacity());
    excess[edge.from] -= flow;
    excess[edge.to] += flow;
    it.pushFlow(flow);
}

template<typename T>
void relabel(EdgesContainer<T>& edges_container, vector<int>& height, int vertex_number) {
    int new_height = INF_int;
    for (typename EdgesContainer<T>::EdgeIterator it =
            edges_container.Begin(vertex_number); it.isValid(); ++it) {
        Edge<T> current_edge = it.getEdge();
        if (current_edge.remainingCapacity() > 0) {
            new_height = min(new_height, height[current_edge.to] + 1);
        }
    }
    height[vertex_number] = new_height;
}

template<typename T>
bool discharge(EdgesContainer<T>& edges_container, vector<T>& excess, vector<int>& height,
        vector<typename EdgesContainer<T>::EdgeIterator>& iterators, int vertex_number) {
    bool height_was_changed = false;

    while(excess[vertex_number] > 0) {
        if (!iterators[vertex_number].isValid()) {
            relabel(edges_container, height, vertex_number);
            height_was_changed = true;
            iterators[vertex_number] = edges_container.Begin(vertex_number);
        } else {
            if (canPush(excess, height, iterators[vertex_number].getEdge())) {
                push(edges_container, excess, iterators[vertex_number]);
            } else {
                ++iterators[vertex_number];
            }
        }
    }

    return height_was_changed;
}

template<typename T>
T getMaxFlowPushRelabel(EdgesContainer<T>& edges_container) {
    int n = edges_container.getVertexNumber();
    int s = edges_container.getSource();
    int t = edges_container.getSink();

    vector<int> height(n, 0);
    height[s] = n;

    vector<T> excess(n, T(0));
    for (typename EdgesContainer<T>::EdgeIterator it = edges_container.Begin(s);
            it.isValid(); ++it) {
        Edge<T> current_edge = it.getEdge();
        it.pushFlow(current_edge.capacity);
        excess[current_edge.to] += current_edge.capacity;
    }

    list<int> vertex_list;
    for (int i = 0; i < n; ++i) {
        if (i != s && i != t) {
            vertex_list.push_back(i);
        }
    }

    vector<typename EdgesContainer<T>::EdgeIterator> iterators(n);
    for (int i = 0; i < n; ++i) {
        iterators[i] = edges_container.Begin(i);
    }

    while(true) {
        bool finish_reached = true;
        for (list<int>::iterator it = vertex_list.begin(); it !=
                vertex_list.end(); ++it) {
            int vertex_number = *it;
            if (excess[vertex_number] > 0) {
                bool height_was_changed = discharge(edges_container,
                        excess, height, iterators, vertex_number);
                if (height_was_changed) {
                    vertex_list.erase(it);
                    vertex_list.push_front(vertex_number);
                    finish_reached = false;
                    break;
                }
            }
        }
        if (finish_reached) {
            break;
        }
    }

    return edges_container.getCurrentFlow();
}

/// Push relabel algorithm part is over

const long long delta_value = 1001;

EdgesContainer<long long> getEdgesContainer(int n, const vector<long long>& value,
        const vector<vector<int>>& dependencies) {
    EdgesContainer<long long> edges_container(n + 2);
    int s = n;
    int t = n + 1;
    edges_container.setSource(n);
    edges_container.setSink(n + 1);
    for (int i = 0; i < n; ++i) {
        edges_container.addEdge(i, t, delta_value - value[i]);
        edges_container.addEdge(s, i, delta_value);
    }
    for (int i = 0; i < n; ++i) {
        for (unsigned j = 0; j < dependencies[i].size(); ++j) {
            edges_container.addEdge(i, dependencies[i][j], INF_long_long);
        }
    }
    return edges_container;
}

void readInput(int& n, vector<long long>& value,
        vector<vector<int>>& dependencies, istream& in) {
    in >> n;

    value = vector<long long>(n);
    for (int i = 0; i < n; ++i) {
        in >> value[i];
    }

    dependencies = vector<vector<int>>(n);
    for (int i = 0; i < n; ++i) {
        int sz;
        in >> sz;
        dependencies[i] = vector<int>(sz);
        for (int j = 0; j < sz; ++j) {
            in >> dependencies[i][j];
            --dependencies[i][j];
        }
    }
}

template<typename A>
void solveProblem(istream& in, ostream& out, A algorithm) {
    int n;
    vector<long long> value;
    vector<vector<int>> dependencies;
    readInput(n, value, dependencies, in);

    EdgesContainer<long long> edges_container = getEdgesContainer(n, value, dependencies);
    long long flow = algorithm(edges_container);
    out << - (flow - delta_value * n);
}

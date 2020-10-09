#include <random>

double randomDouble01() {
    double random_value = rand();
    double max_value = 32768;
    return random_value / max_value;
}

int randomIntBetween(int l, int r) {
    return randomDouble01() * (r - l) + l + 1;
}

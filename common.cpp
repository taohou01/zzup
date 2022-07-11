#include "common.h"

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const DpcEvent &evt) {
    const char* evt_name[] = 
        {"l_min", "r_min", "l_max", "r_max", "min", "max", "inc", "dec", "oppo"};

    os << std::setprecision(std::numeric_limits<Decimal>::digits10) << evt.d  << " " 
        << evt.t << " " << evt_name[evt.type] << " " << evt.e1;

    if (evt.type > LOCAL_MAX) {
        os << " " << evt.e2;
    }

    return os;
}

Integer combChoose(Integer n, Integer k) {
    if (k == 0) { return 1; }
    return (n * combChoose(n - 1, k - 1)) / k;
}

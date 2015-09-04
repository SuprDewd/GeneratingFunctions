#include <bits/stdc++.h>
#include "gmpxx.h"
using namespace std;
typedef long long ll;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;
template <class T> int size(const T &x) { return x.size(); }
const int INF = 2147483647;
#define rep(i,a,b) for (__typeof(a) i=(a); i<(b); ++i)
#define iter(it,c) for (__typeof((c).begin()) it = (c).begin(); it != (c).end(); ++it)

mpq_class tompq(double x) { return mpq_class(x * 1e9) / mpz_class("1000000000"); }
mpq_class mpq_log(mpq_class z) { return tompq(log(z.get_d())); }
mpq_class mpq_exp(mpq_class z) { return tompq(exp(z.get_d())); }
mpq_class mpq_sin(mpq_class z) { return tompq(sin(z.get_d())); }
mpq_class mpq_cos(mpq_class z) { return tompq(cos(z.get_d())); }

#define DEFAULT_CAP 15

struct power_series {
    int cap;
    vector<mpq_class> a;
    power_series(int cap=DEFAULT_CAP);
    power_series(int cap, mpq_class c);
    void clean();
    power_series with_cap(int cap) const;
    mpq_class head() const;
    power_series tail() const;
    power_series pow(mpz_class p) const;
    power_series integral() const;
    power_series derivative() const;
    power_series log() const;
    power_series exp() const;
    power_series operator +(const power_series &other) const;
    power_series operator -(const power_series &other) const;
    power_series operator -() const;
    power_series operator *(const mpq_class &c) const;
    power_series operator *(const power_series &other) const;
    power_series operator /(const power_series &other) const;
};
ostream& operator <<(ostream &outs, const power_series &p) {
    // if (size(p.a) == 0) {
    //     outs << "0";
    // } else {
    //     int i = 0;
    //     while (i < size(p.a)) {
    //         if (i != 0 || p.a[i] < 0) {
    //             if (p.a[i] < 0) {
    //                 outs << "-";
    //             } else {
    //                 outs << "+";
    //             }
    //         }
    //         outs << abs(p.a[i]);
    //         if (i > 0) {
    //             outs << "x";
    //             if (i > 1) {
    //                 outs << "^" << i;
    //             }
    //         }
    //         i++;
    //     }
    // }
    outs << "{";
    int i = 0;
    while (i < p.cap) {
        if (i > 0) {
            outs << ",";
        }
        if (i < size(p.a)) {
            outs << p.a[i];
        } else {
            outs << "0";
        }
        i++;
    }
    outs << "}";
    return outs;
}

power_series X(DEFAULT_CAP);

power_series::power_series(int cap) {
    this->cap = max(cap,0);
}
power_series::power_series(int cap, mpq_class c) {
    this->cap = max(cap,0);
    a.push_back(c);
}
void power_series::clean() {
    while (size(a) > cap) {
        a.pop_back();
    }
    while (!a.empty() && a.back() == 0) {
        a.pop_back();
    }
}
power_series power_series::with_cap(int cap) const {
    power_series res(max(cap,0));
    res.a = a;
    res.clean();
    return res;
}
mpq_class power_series::head() const {
    if (size(a) >= 1) {
        return a[0];
    } else {
        return 0;
    }
}
power_series power_series::tail() const {
    power_series res(cap-1);
    int i = 1;
    while (i < size(a)) {
        res.a.push_back(a[i]);
        i++;
    }
    return res;
}
power_series power_series::pow(mpz_class p) const {
    power_series a = *this;
    power_series res(cap, 1);
    while (p > 0) {
        if ((p & 1) != 0) {
            res = res * a;
        }
        a = a * a;
        p >>= 1;
    }
    return res;
}
power_series power_series::integral() const {
    power_series res(cap+1);
    res.a.push_back(0);
    int i = 0;
    while (i < size(a)) {
        res.a.push_back(a[i] / (i+1));
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::derivative() const {
    power_series res(cap-1);
    int i = 1;
    while (i < size(a)) {
        res.a.push_back(a[i] * i);
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::log() const {
    power_series res = (this->derivative() / (*this)).integral();
    res.a[0] = mpq_log(a[0]);
    return res.with_cap(cap);
}
power_series power_series::exp() const {
    if (size(a) == 0) {
        return power_series(cap,1);
    }
    power_series res = (with_cap(cap-1).exp() * derivative()).integral();
    if (size(res.a) == 0) res.a.push_back(0);
    res.a[0] = mpq_exp(a[0]);
    res.clean();
    return res;
}
power_series power_series::operator +(const power_series &other) const {
    power_series res(max(cap,other.cap));
    int i = 0;
    while (i < size(a) || i < size(other.a)) {
        mpq_class here = 0;
        if (i < size(a)) here += a[i];
        if (i < size(other.a)) here += other.a[i];
        res.a.push_back(here);
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::operator -(const power_series &other) const {
    power_series res(max(cap,other.cap));
    int i = 0;
    while (i < size(a) || i < size(other.a)) {
        mpq_class here = 0;
        if (i < size(a)) here += a[i];
        if (i < size(other.a)) here -= other.a[i];
        res.a.push_back(here);
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::operator -() const {
    power_series res(cap);
    int i = 0;
    while (i < size(a)) {
        res.a.push_back(-a[i]);
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::operator *(const mpq_class &c) const {
    power_series res(cap);
    int i = 0;
    while (i < size(a)) {
        res.a.push_back(a[i] * c);
        i++;
    }
    res.clean();
    return res;
}
power_series power_series::operator *(const power_series &other) const {
    if (size(a) == 0) return power_series(max(cap,other.cap));
    if (size(a) == 2 && a[0] == 0 && a[1] == 1) {
        power_series res(max(cap,other.cap));
        res.a.push_back(0);
        int i = 0;
        while (i < size(other.a)) {
            res.a.push_back(other.a[i]);
            i++;
        }
        res.clean();
        return res;
    }
    return power_series(max(cap,other.cap), head() * other.head()) + X.with_cap(max(cap,other.cap)) * (other.tail() * head() + tail() * other.with_cap(other.cap - 1));
}
power_series power_series::operator /(const power_series &other) const {
    if (size(a) == 0) return power_series();
    return power_series(max(cap,other.cap), head() / other.head()) + X.with_cap(max(cap,other.cap)) * ((tail() - other.tail() * (head() / other.head())).with_cap(cap-1) / other);
}

int main() {
    X.a.push_back(0);
    X.a.push_back(1);

    // power_series p;
    // p.a.push_back(1);
    // p.a.push_back(1);
    // cout << p.pow(10000000) << endl;

    // cout << (power_series(DEFAULT_CAP, 1) / (power_series(DEFAULT_CAP, 1) - X)).pow(1000) << endl;

    cout << (power_series(DEFAULT_CAP, 2) / (power_series(DEFAULT_CAP, 1) - X)).log().exp() << endl;

    // power_series res = p;
    // rep(i,0,100) {
    //     cout << res << endl;
    //     res = res * p;
    // }
    // rep(i,0,100) {
    //     res = res / p;
    //     cout << res << endl;
    // }
    // cout << power_series(DEFAULT_CAP, 1) / (power_series(DEFAULT_CAP, 1) - X) << endl;
    return 0;
}


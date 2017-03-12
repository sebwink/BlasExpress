#include <cassert>

#include <iostream>
#include <vector>
#include <type_traits>

template <typename E>
class VecExpr {

public:

	auto operator[](size_t i) {
		return static_cast<const E&>(*this)[i]; 
	}

	size_t size() const {
		return static_cast<const E&>(*this).size();
	}

	E& operator ()() {
		return static_cast<E&>(*this);
	}

	const E& operator()() const {
		return static_cast<const E&>(*this);
	}
};

template <typename ValueType>
class Vec : public VecExpr<Vec<ValueType>> {
	
	std::vector<ValueType> elems;

public:

	using vtype = ValueType;

	ValueType operator[](size_t i) const {
		return elems[i];
	}

	ValueType& operator[](size_t i) {
		return elems[i];
	}

	size_t size() const {
		return elems.size();
	}

	template <typename T>
	friend std::ostream& operator<<(std::ostream&, const Vec<T>&);

	Vec(size_t n) : elems(n) { }

	Vec(std::vector<ValueType> v) : elems(v) { }

	template <typename E>
	Vec(const VecExpr<E>& vec) : elems(vec.size()) {
		for (size_t i = 0; i != vec.size(); ++i) {
			elems[i] = vec[i];
		}
	}

};

template <typename ValueType>
std::ostream& operator<<(std::ostream& os, const Vec<ValueType>& v) {
	for (const ValueType& elem : v.elems)
		os << elem;
	return os;
}

template <typename A, typename B>
class VecSum : public VecExpr<VecSum<A,B>> {

	
	const A& _a;
	const B& _b;

public:

	VecSum(const A& a, const B& b) : _a(a), _b(b) {
		assert(a.size() == b.size());
	}

	auto operator[](size_t i) const { 
		return _a[i] + _b[i]; 
	}

	auto operator[](size_t i) {
		return _a[i] + _b[i];
	}

	size_t size() const {
		return _a.size();
	} 

};

template <typename A, typename B>
VecSum<A,B> const operator+(const A& a, const B& b) {
	return VecSum<A,B>(a,b);
}

int main() {
	std::vector<double> stdv {1.0,1.0, 1.0};
	Vec<double> v(stdv);
	Vec<double> u(stdv);
	Vec<double> w(stdv);
	Vec<double> x = u + v + w;
}

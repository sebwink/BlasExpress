#include <iostream>

using std::cout;
using std::endl;

#include <boost/proto/proto.hpp>

using boost::proto::_;
using boost::proto::or_;
using boost::proto::when;
using boost::proto::_value;
using boost::proto::_left;
using boost::proto::_right;
using boost::proto::otherwise;
using boost::proto::_default;
using boost::proto::terminal;
using boost::proto::plus;
using boost::proto::multiplies;

struct eval : boost::proto::or_< when< terminal<_>
                                     , _value>
                               , otherwise< _default<eval> >
							   >
{ };

struct epa : boost::proto::or_< terminal<_>
                              , plus<epa,epa>
							  , multiplies<epa,epa>
							  >
{ };

template <typename AST> struct expr_;

struct domain_ : boost::proto::domain< boost::proto::generator<expr_>, epa > { };
	
template <typename AST>
struct expr_ : boost::proto::extends<AST, expr_<AST>, domain_> {
	using extendee = boost::proto::extends<AST, expr_<AST>, domain_>;

	expr_(const AST& ast = AST()) : extendee(ast) { }

	BOOST_PROTO_EXTENDS_USING_ASSIGN(expr_)

	using result_type = double;
	
	result_type operator()() const {
		eval callee;
		return callee(*this);
	}

};

template <typename T>
struct variable : expr_<typename terminal<T>::type> {
	
	variable(const T& v = T{}) {
		boost::proto::value(*this) = v; 
	}

	template <typename X>
	variable operator=(const expr_<X>&  x) {
		boost::proto::value(*this) = x();
		return *this;
	}
};


struct Int {
	variable<int> repr;

    Int(){}
	Int(int i) {
		repr = variable<int>(i);
	}
	template<typename T>
	Int(expr_<T> xrepr) {
		repr = xrepr;
	}	

	template<typename T>
	Int& operator=(const expr_<T>& xrepr) {
		repr = xrepr;
		return *this;
	}

	friend std::ostream& operator<<(std::ostream&, Int);
};

std::ostream& operator<<(std::ostream& out, Int i) {
	out << i.repr();
	return out;
}

Int operator+(Int i, Int j) {
	return Int(i.repr + j.repr);
}

Int operator*(Int i, Int j) {
	return Int(i.repr * j.repr);
}

Int operator+(Int i, int j) {
	return Int(i.repr + j);
}

int main() {
	Int x(2);
	Int y(x+x*x+3);
	cout << y << endl;
}	
